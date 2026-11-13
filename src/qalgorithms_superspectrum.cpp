#include "qalgorithms_superspectrum.h"
#include "qalgorithms_superspectrum_tracks.h"
#include "qalgorithms_input_output.h"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>

namespace qAlgorithms
{
    namespace
    {
        constexpr size_t scansPerWindow = 5;
        using GridIndex = int64_t;
        constexpr GridIndex maxGridPoints = 1000000;

        struct GriddedRun
        {
            GridIndex first;
            std::vector<float> logs;
        };

        struct ProfileRun
        {
            size_t scan, first, end;
            double lower, upper;
        };

        struct ProfileScan
        {
            uint32_t index;
            float retentionTime;
            std::vector<double> mz, intensity;
            std::vector<ProfileRun> runs;
            std::map<std::pair<size_t, double>, GriddedRun> grid;
        };

        int massRegion(double mz)
        {
            int exponent;
            const double fraction = std::frexp(mz, &exponent);
            return 4 * exponent + static_cast<int>(8 * fraction) - 4;
        }

        double regionLower(int region)
        {
            const int band = (region % 4 + 4) % 4;
            return std::ldexp(.5 + band / 8.0, (region - band) / 4);
        }

        struct GridSettings
        {
            SuperspectrumSettings fallback{};
            std::map<int, SuperspectrumSettings> regions;

            const SuperspectrumSettings &at(double mz) const
            {
                if (regions.empty())
                    return fallback;
                const auto next = regions.find(massRegion(mz));
                if (next != regions.end())
                    return next->second;
                return fallback;
            }
        };

        ProfileScan readProfile(const XML_File &file, uint32_t index, float rt)
        {
            ProfileScan scan{index, rt, {}, {}, {}, {}};
            if (get_spectrum(&file, &scan.mz, &scan.intensity, scan.index) != 0 ||
                scan.mz.size() != scan.intensity.size() || (!scan.mz.empty() && scan.mz.front() <= 0) ||
                !std::is_sorted(scan.mz.begin(), scan.mz.end()) ||
                std::adjacent_find(scan.mz.begin(), scan.mz.end()) != scan.mz.end())
                throw std::runtime_error("Invalid profile data for superspectrum scan " + std::to_string(scan.index));
            return scan;
        }

        struct ProfileMeasurements
        {
            std::vector<double> spacing, widths;
        };

        double quantile(const std::vector<double> &sorted, double fraction)
        {
            const double index = fraction * (sorted.size() - 1);
            const auto lower = static_cast<size_t>(index);
            const auto upper = static_cast<size_t>(std::ceil(index));
            return sorted[lower] + (index - lower) * (sorted[upper] - sorted[lower]);
        }

        SuperspectrumSettings estimateSpacing(std::vector<double> &spacing)
        {
            std::sort(spacing.begin(), spacing.end());
            const double median = quantile(spacing, .5);
            const double step = std::max(median / 2, std::numeric_limits<double>::denorm_min());
            const double maxGap = std::max(3 * median, 1.5 * quantile(spacing, .9));
            return {step, std::min(maxGap, std::numeric_limits<double>::max()), maxscale_global};
        }

        size_t estimateHalfWidth(std::vector<double> &widths, double step)
        {
            if (widths.size() < 12)
                return maxscale_global;
            std::sort(widths.begin(), widths.end());
            const double halfWidth = std::ceil(1.5 * quantile(widths, .9) / step);
            return static_cast<size_t>(std::clamp(halfWidth, 2.0, static_cast<double>(maxscale_global)));
        }

        struct GridEvent
        {
            GridIndex index;
            GridIndex end;
            const ProfileRun *run;
            bool begin;
        };

        struct Candidate
        {
            SharedShapeRegCoeffs fit;
            double apex, logRise, score;
        };

        struct WindowCandidate
        {
            Candidate candidate;
            GridIndex center;
            double step;
            std::vector<size_t> scans;
        };

        double mass(GridIndex index, double step)
        {
            return static_cast<double>(index) * step;
        }

        GridIndex gridIndex(double mz, double step, bool lowerBound)
        {
            const double quotient = mz / step;
            // Consecutive indices must remain distinguishable as doubles.
            if (!(quotient < 0x1p52))
                throw std::invalid_argument("Superspectrum m/z step is too small for the measured range");
            auto index = static_cast<GridIndex>(lowerBound ? std::ceil(quotient) : std::floor(quotient));
            if (lowerBound)
            {
                if (mass(index, step) < mz)
                    ++index;
                else if (index > 0 && mass(index - 1, step) >= mz)
                    --index;
            }
            else if (mass(index, step) > mz)
                --index;
            else if (mass(index + 1, step) <= mz)
                ++index;
            return index;
        }

        void findProfileRuns(ProfileScan &profile, const GridSettings &settings)
        {
            size_t first = 0;
            while (first < profile.mz.size())
            {
                if (profile.intensity[first] <= 0)
                {
                    ++first;
                    continue;
                }
                size_t end = first + 1;
                while (end < profile.mz.size() && profile.intensity[end] > 0 &&
                       profile.mz[end] - profile.mz[end - 1] <=
                           settings.at(profile.mz[end - 1] + (profile.mz[end] - profile.mz[end - 1]) / 2).maxGap)
                    ++end;
                if (end - first >= 2)
                    profile.runs.push_back({0, first, end, profile.mz[first], profile.mz[end - 1]});
                first = end;
            }
        }

        GridSettings estimateGrid(const XML_File &file, const std::vector<uint32_t> &indices, bool silent)
        {
            constexpr size_t sampleLimit = 24;
            constexpr size_t minimumRegionalSpacings = 24;
            std::map<int, ProfileMeasurements> regional;
            ProfileMeasurements pooled;
            std::vector<ProfileScan> sample;
            const size_t count = std::min(sampleLimit, indices.size());
            auto observe = [&](ProfileScan scan)
            {
                const auto &mz = scan.mz;
                for (size_t i = 1; i < mz.size(); ++i)
                {
                    const double distance = mz[i] - mz[i - 1];
                    double local = distance;
                    if (i > 1)
                        local = std::min(local, mz[i - 1] - mz[i - 2]);
                    if (i + 1 < mz.size())
                        local = std::min(local, mz[i + 1] - mz[i]);
                    if (distance <= 2 * local)
                    {
                        regional[massRegion(mz[i - 1] + distance / 2)].spacing.push_back(distance);
                        pooled.spacing.push_back(distance);
                    }
                }
                sample.push_back(std::move(scan));
            };
            for (size_t i = 0; i < count; ++i)
            {
                const size_t intervals = std::max<size_t>(1, count - 1);
                const size_t last = indices.size() - 1;
                const size_t position = i * (last / intervals) + i * (last % intervals) / intervals;
                observe(readProfile(file, indices[position], 0));
            }
            if (pooled.spacing.empty() && count < indices.size())
            {
                for (uint32_t index : indices)
                {
                    sample.clear();
                    observe(readProfile(file, index, 0));
                    if (!pooled.spacing.empty())
                        break;
                }
            }
            if (pooled.spacing.empty())
                return {};

            GridSettings result{estimateSpacing(pooled.spacing), {}};
            for (auto &[region, data] : regional)
                if (data.spacing.size() >= minimumRegionalSpacings)
                    result.regions.emplace(region, estimateSpacing(data.spacing));
            for (auto &scan : sample)
            {
                findProfileRuns(scan, result);
                for (const auto &run : scan.runs)
                    for (size_t i = run.first + 1; i + 1 < run.end; ++i)
                    {
                        const auto &y = scan.intensity;
                        if (!(y[i] > y[i - 1] && y[i] > y[i + 1]))
                            continue;
                        const double half = y[i] / 2;
                        size_t left = i, right = i;
                        while (left > run.first && y[left] > half && y[left - 1] < y[left])
                            --left;
                        while (right + 1 < run.end && y[right] > half && y[right + 1] < y[right])
                            ++right;
                        if (left == i || right == i || y[left] > half || y[right] > half)
                            continue;
                        const auto &mz = scan.mz;
                        const double lower = mz[left] + (mz[left + 1] - mz[left]) *
                            ((half - y[left]) / (y[left + 1] - y[left]));
                        const double upper = mz[right] - (mz[right] - mz[right - 1]) *
                            ((half - y[right]) / (y[right - 1] - y[right]));
                        const double width = upper - lower;
                        if (width > 0)
                        {
                            regional[massRegion(mz[i])].widths.push_back(width);
                            pooled.widths.push_back(width);
                        }
                    }
            }
            result.fallback.maxHalfWidth = estimateHalfWidth(pooled.widths, result.fallback.gridStep);
            for (auto &[region, settings] : result.regions)
                settings.maxHalfWidth = estimateHalfWidth(regional.at(region).widths, settings.gridStep);
            if (!silent)
            {
                auto report = [](const SuperspectrumSettings &settings)
                {
                    printf("step=%.9g, gap=%.9g, half-width=%zu\n",
                        settings.gridStep, settings.maxGap, settings.maxHalfWidth);
                };
                printf("    Superspectrum pooled settings: ");
                report(result.fallback);
                for (const auto &[region, settings] : result.regions)
                {
                    printf("    Superspectrum m/z [%.9g, %.9g): ", regionLower(region), regionLower(region + 1));
                    report(settings);
                }
            }
            return result;
        }

        const float *gridTrace(ProfileScan &scan, const ProfileRun &run, GridIndex start, double step)
        {
            const auto key = std::make_pair(run.first, step);
            const auto cached = scan.grid.find(key);
            if (cached != scan.grid.end())
                return cached->second.logs.data() + (start - cached->second.first);
            const auto first = gridIndex(run.lower, step, true);
            const auto last = gridIndex(run.upper, step, false);
            if (last - first >= maxGridPoints)
                throw std::invalid_argument("Superspectrum grid run exceeds one million points; increase GRID_STEP");
            const auto length = static_cast<size_t>(last - first + 1);
            std::vector<float> logs;
            logs.reserve(length);
            size_t right = run.first + 1;
            for (size_t i = 0; i < length; ++i)
            {
                const double mz = mass(first + static_cast<GridIndex>(i), step);
                while (right + 1 < run.end && scan.mz[right] < mz)
                    ++right;
                const double fraction = (mz - scan.mz[right - 1]) / (scan.mz[right] - scan.mz[right - 1]);
                const double low = scan.intensity[right - 1];
                const double high = scan.intensity[right];
                const double intensity = fraction <= .5 ? std::fma(fraction, high - low, low) :
                    std::fma(1 - fraction, low - high, high);
                logs.push_back(static_cast<float>(std::log(intensity)));
            }
            const auto inserted = scan.grid.emplace(key, GriddedRun{first, std::move(logs)}).first;
            return inserted->second.logs.data() + (start - first);
        }

        std::vector<Candidate> findCandidates(const std::vector<const float *> &traces, size_t length,
            GridIndex gridStart, const SuperspectrumSettings &settings, const GridSettings &grid)
        {
            const size_t maxHalfWidth = settings.maxHalfWidth;
            std::vector<double> mean(length, 0);
            for (const auto *trace : traces)
                for (size_t i = 0; i < length; ++i)
                    mean[i] += static_cast<double>(trace[i]) / traces.size();

            // Center the trace mean before forming prefix sums for the null fit.
            std::vector<long double> sums(length + 1), squares(length + 1);
            for (size_t i = 0; i < length; ++i)
            {
                const long double value = mean[i] - mean.front();
                sums[i + 1] = sums[i] + value;
                squares[i + 1] = squares[i] + value * value;
            }
            struct Moments { long double second = 0, third = 0, fourth = 0; };
            std::vector<Moments> moments(maxHalfWidth + 1);
            for (size_t i = 1; i <= maxHalfWidth; ++i)
            {
                const long double x = i;
                moments[i] = {moments[i - 1].second + x * x,
                              moments[i - 1].third + x * x * x,
                              moments[i - 1].fourth + x * x * x * x};
            }
            std::vector<std::optional<Candidate>> best(length);
            forEachCoefficientsMultiTrace(traces, length, [&](const SharedShapeRegCoeffs &fit)
            {
                if (!(fit.b2 < 0 && fit.b3 < 0))
                    return;
                const double curvature = fit.b1 < 0 ? fit.b2 : fit.b3;
                const double apex = -fit.b1 / (2 * curvature);
                const double scale = fit.scale;
                if (!(std::abs(apex) < scale - 1))
                    return;
                const double logRise = -curvature * apex * apex -
                    std::max(-fit.b1 * scale + fit.b2 * scale * scale,
                              fit.b1 * scale + fit.b3 * scale * scale);
                if (!(logRise >= std::log(2.0)))
                    return;

                // OLS explained sum of squares on the mean trace. The
                // intercept drops out when the fitted values are centered.
                const size_t first = fit.x0 - fit.scale;
                const size_t end = fit.x0 + fit.scale + 1;
                const long double n = end - first;
                const long double sum = sums[end] - sums[first];
                const long double baseline = squares[end] - squares[first] - sum * sum / n;
                const auto &m = moments[fit.scale];
                const long double b1 = fit.b1, b2 = fit.b2, b3 = fit.b3;
                const long double fittedSum = (b2 + b3) * m.second;
                const long double explained = 2 * b1 * b1 * m.second + 2 * b1 * (b3 - b2) * m.third
                    + (b2 * b2 + b3 * b3) * m.fourth - fittedSum * fittedSum / n;
                if (!(baseline > 0 && explained > baseline / 2 && explained <= baseline * 1.0001L))
                    return;
                if (!grid.regions.empty())
                {
                    const double center = mass(gridStart + static_cast<GridIndex>(fit.x0), settings.gridStep);
                    const double apexMz = std::fma(apex, settings.gridStep, center);
                    if (&grid.at(apexMz) != &settings)
                        return;
                }
                auto &candidate = best[static_cast<size_t>(std::llround(static_cast<double>(fit.x0) + apex))];
                const double score = static_cast<double>(explained);
                if (!candidate || score > candidate->score)
                    candidate = Candidate{fit, apex, logRise, score};
            }, maxHalfWidth);
            std::vector<Candidate> candidates;
            for (auto &candidate : best)
                if (candidate)
                    candidates.push_back(std::move(*candidate));
            return candidates;
        }

        void fitBlock(std::vector<WindowCandidate> &candidates, std::deque<ProfileScan> &scans,
            const std::array<const ProfileRun *, scansPerWindow> &active, GridIndex first, GridIndex end,
            const SuperspectrumSettings &settings, const GridSettings &grid)
        {
            if (end - first < 5)
                return;
            std::vector<size_t> scanIds;
            for (size_t scan = 0; scan < active.size(); ++scan)
                if (active[scan])
                    scanIds.push_back(scan);
            if (scanIds.size() < 2)
                return;
            const auto length = static_cast<size_t>(end - first);
            std::vector<const float *> logs;
            for (size_t scan : scanIds)
                logs.push_back(gridTrace(scans[scan], *active[scan], first, settings.gridStep));
            for (auto &candidate : findCandidates(logs, length, first, settings, grid))
            {
                const auto center = first + static_cast<GridIndex>(candidate.fit.x0);
                const double apexOffset = candidate.fit.b1 * candidate.apex / 2;
                const bool validHeight = std::all_of(candidate.fit.b0.begin(), candidate.fit.b0.end(),
                    [apexOffset](double intercept)
                    {
                        const double height = std::exp(intercept + apexOffset);
                        return height > 0 && std::isfinite(height);
                    });
                if (validHeight)
                    candidates.push_back({std::move(candidate), center, settings.gridStep, scanIds});
            }
        }

        void selectCandidates(std::vector<WindowCandidate> &selected, std::vector<WindowCandidate> &candidates)
        {
            std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b)
            {
                // Compare explained variation integrated over m/z so a
                // denser lattice has no advantage from sample count alone.
                const long double scoreA = static_cast<long double>(a.candidate.score) * a.scans.size() * a.step;
                const long double scoreB = static_cast<long double>(b.candidate.score) * b.scans.size() * b.step;
                if (scoreA != scoreB)
                    return scoreA > scoreB;
                const double centerA = mass(a.center, a.step);
                const double centerB = mass(b.center, b.step);
                if (centerA != centerB)
                    return centerA < centerB;
                if (a.candidate.fit.scale != b.candidate.fit.scale)
                    return a.candidate.fit.scale < b.candidate.fit.scale;
                return a.scans < b.scans;
            });
            std::map<double, double> occupied;
            for (auto &window : candidates)
            {
                const auto &candidate = window.candidate;
                const auto &fit = candidate.fit;
                const double step = window.step;
                const double first = mass(window.center - fit.scale, step);
                const double last = mass(window.center + fit.scale, step);
                const auto next = occupied.lower_bound(first);
                if ((next != occupied.end() && next->first <= last) ||
                    (next != occupied.begin() && std::prev(next)->second >= first))
                    continue;
                occupied.emplace_hint(next, first, last);
                selected.push_back(std::move(window));
            }
        }

        void writeWindow(std::ostream &csv, SuperspectrumTracks &tracks,
            const std::deque<ProfileScan> &scans, const std::vector<WindowCandidate> &selected)
        {
            std::vector<SuperspectrumObservation> observations;
            observations.reserve(selected.size());
            for (const auto &window : selected)
            {
                const auto &candidate = window.candidate;
                const auto &fit = candidate.fit;
                const double step = window.step;
                const double center = mass(window.center, step);
                const double apexMz = std::fma(candidate.apex, step, center);
                SuperspectrumObservation observation{apexMz, mass(window.center - fit.scale, step),
                    mass(window.center + fit.scale, step), step, {}};
                for (size_t i = 0; i < window.scans.size(); ++i)
                {
                    const size_t position = window.scans[i];
                    const auto &scan = scans[position];
                    const double height = std::exp(fit.b0[i] + fit.b1 * candidate.apex / 2);
                    const size_t middle = scansPerWindow / 2;
                    const size_t distance = position < middle ? middle - position : position - middle;
                    observation.heights.push_back({scan.index, scan.retentionTime, height, distance});
                }
                observations.push_back(std::move(observation));
            }
            const auto ids = tracks.append(observations, scans.front().index);
            for (size_t i = 0; i < selected.size(); ++i)
            {
                const auto &window = selected[i];
                const auto &candidate = window.candidate;
                const auto &fit = candidate.fit;
                const double step = window.step;
                const double center = mass(window.center, step);
                for (const auto &height : observations[i].heights)
                {
                    csv << scans.front().index << ',' << scans.back().index << ',' << i + 1 << ','
                        << height.scan << ',' << height.retentionTime << ',' << observations[i].mz << ',' << height.height << ','
                        << center << ',' << fit.scale << ',' << fit.b1 << ',' << fit.b2 << ',' << fit.b3
                        << ',' << candidate.logRise << ',' << step << ',' << ids[i] << '\n';
                }
            }
        }

        void fitComponent(std::vector<WindowCandidate> &candidates, std::deque<ProfileScan> &scans,
            const std::vector<ProfileRun> &runs, size_t first, size_t end,
            const SuperspectrumSettings &settings, const GridSettings &grid,
            double lower, double upper)
        {
            std::vector<GridEvent> events;
            for (size_t i = first; i < end; ++i)
            {
                if (runs[i].upper < lower || runs[i].lower > upper)
                    continue;
                const auto begin = gridIndex(std::max(runs[i].lower, lower), settings.gridStep, true);
                const auto last = gridIndex(std::min(runs[i].upper, upper), settings.gridStep, false);
                if (begin <= last)
                {
                    events.push_back({begin, last + 1, &runs[i], true});
                    events.push_back({last + 1, last + 1, &runs[i], false});
                }
            }
            std::sort(events.begin(), events.end(), [](const auto &a, const auto &b)
            {
                if (a.index != b.index)
                    return a.index < b.index;
                return a.begin < b.begin;
            });
            std::array<const ProfileRun *, scansPerWindow> active{};
            std::array<GridIndex, scansPerWindow> endAt{};
            for (size_t i = 0; i < events.size();)
            {
                const GridIndex index = events[i].index;
                std::array<bool, scansPerWindow> starts{};
                do
                {
                    const auto &event = events[i++];
                    const size_t scan = event.run->scan;
                    active[scan] = event.begin ? event.run : nullptr;
                    starts[scan] = event.begin;
                    if (event.begin)
                        endAt[scan] = event.end;
                } while (i < events.size() && events[i].index == index);

                // Each common interval starts at one run boundary and
                // ends at another. A short extra trace must not truncate
                // a fit supported across the full peak by other scans.
                if (std::none_of(starts.begin(), starts.end(), [](bool start) { return start; }))
                    continue;
                std::array<GridIndex, scansPerWindow> seen{};
                size_t count = 0;
                for (size_t edge = 0; edge < active.size(); ++edge)
                {
                    const auto boundary = endAt[edge];
                    if (!active[edge] || std::find(seen.begin(), seen.begin() + count, boundary) != seen.begin() + count)
                        continue;
                    seen[count++] = boundary;
                    std::array<const ProfileRun *, scansPerWindow> covering{};
                    bool startsHere = false;
                    for (size_t scan = 0; scan < active.size(); ++scan)
                        if (active[scan] && endAt[scan] >= boundary)
                        {
                            covering[scan] = active[scan];
                            startsHere |= starts[scan];
                        }
                    if (startsHere)
                        fitBlock(candidates, scans, covering, index, boundary, settings, grid);
                }
            }
        }

        void fitWindow(std::ostream &csv, SuperspectrumTracks &tracks, std::deque<ProfileScan> &scans,
            const GridSettings &grid)
        {
            std::vector<ProfileRun> runs;
            for (size_t scan = 0; scan < scans.size(); ++scan)
                for (auto run : scans[scan].runs)
                {
                    run.scan = scan;
                    runs.push_back(run);
                }
            std::stable_sort(runs.begin(), runs.end(), [](const auto &a, const auto &b) { return a.lower < b.lower; });
            std::vector<WindowCandidate> selected;
            for (size_t first = 0; first < runs.size();)
            {
                size_t end = first + 1;
                double upper = runs[first].upper;
                while (end < runs.size() && runs[end].lower <= upper)
                {
                    upper = std::max(upper, runs[end].upper);
                    ++end;
                }
                const int firstRegion = massRegion(runs[first].lower);
                const int lastRegion = massRegion(upper);
                std::vector<WindowCandidate> candidates;
                size_t regionalGrids = 0;
                auto region = grid.regions.lower_bound(firstRegion);
                for (; region != grid.regions.end() && region->first <= lastRegion; ++region)
                {
                    const auto &settings = region->second;
                    const double halo = 2 * settings.maxHalfWidth * settings.gridStep;
                    fitComponent(candidates, scans, runs, first, end, settings, grid,
                        regionLower(region->first) - halo, regionLower(region->first + 1) + halo);
                    ++regionalGrids;
                }
                if (regionalGrids < static_cast<size_t>(lastRegion - firstRegion + 1))
                    fitComponent(candidates, scans, runs, first, end, grid.fallback, grid, runs[first].lower, upper);
                selectCandidates(selected, candidates);
                first = end;
            }
            writeWindow(csv, tracks, scans, selected);
        }

        struct TrackChromatograms
        {
            std::vector<EIC> bins;
            std::vector<size_t> trackIds;
            std::vector<ChromatogramMass> masses;

            void append(const std::vector<uint32_t> &indices, const std::vector<float> &retentionTimes,
                size_t id, double mz,
                std::optional<double> spread, const std::map<uint32_t, SuperspectrumTrackHeight> &heights)
            {
                auto scan = indices.begin();
                for (auto first = heights.begin(); first != heights.end();)
                {
                    scan = std::lower_bound(scan, indices.end(), first->first);
                    if (scan == indices.end() || *scan != first->first)
                        throw std::runtime_error("Superspectrum track references an unknown profile scan");
                    auto end = std::next(first);
                    auto nextScan = std::next(scan);
                    size_t measured = 1;
                    while (end != heights.end() && nextScan != indices.end())
                    {
                        if (end->first != *nextScan)
                        {
                            if (std::next(nextScan) == indices.end() || end->first != *std::next(nextScan))
                                break;
                            ++nextScan;
                        }
                        ++end;
                        ++nextScan;
                        ++measured;
                    }
                    const auto length = static_cast<size_t>(nextScan - scan);
                    if (measured >= MINLENGTH)
                    {
                        if (length > UINT16_MAX)
                            throw std::runtime_error("Superspectrum chromatogram exceeds the qPeaks sample limit");
                        if (!(mz >= std::numeric_limits<float>::denorm_min() && mz <= std::numeric_limits<float>::max()))
                            throw std::runtime_error("Superspectrum track mass exceeds the qPeaks m/z range");
                        EIC bin;
                        bin.scanNumbers.clear();
                        bin.RT.clear();
                        bin.ints_area.clear();
                        bin.DQSC.clear();
                        bin.df.clear();
                        bin.mz.clear();
                        auto point = first;
                        uint16_t count = 0;
                        for (auto index = scan; index != nextScan; ++index)
                        {
                            const float rt = retentionTimes.at(index - indices.begin());
                            double height;
                            if (*index == point->first)
                            {
                                height = point->second.measurement.height;
                                ++point;
                                ++count;
                            }
                            else
                            {
                                const auto &left = std::prev(point)->second.measurement;
                                const auto &right = point->second.measurement;
                                const double fraction = (double(rt) - left.retentionTime) /
                                    (double(right.retentionTime) - left.retentionTime);
                                height = std::exp((1 - fraction) * std::log(left.height) +
                                    fraction * std::log(right.height));
                            }
                            if (height > std::numeric_limits<float>::max())
                                throw std::runtime_error("Superspectrum height exceeds the qPeaks intensity range");
                            bin.scanNumbers.push_back(*index);
                            bin.RT.push_back(rt);
                            // findFeatures uses ints_area as the chromatogram signal.
                            bin.ints_area.push_back(static_cast<float>(height));
                            bin.df.push_back(count);
                        }
                        bins.push_back(std::move(bin));
                        trackIds.push_back(id);
                        masses.push_back({static_cast<float>(mz), spread ? static_cast<float>(*spread) :
                            std::numeric_limits<float>::quiet_NaN()});
                    }
                    first = end;
                    scan = nextScan;
                }
            }
        };
    }

    void processSuperspectra(const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<float> &retentionTimes, const std::filesystem::path &outputDir,
        const std::string &filename,
        const SuperspectrumSettings *fixed, bool noOverwrite, bool silent)
    {
        const std::string stem = filename + "_superspectrum";
        const std::array<std::string, 4> names{stem + "_centroids.csv", stem + "_tracks.csv",
            stem + "_features.csv", stem + "_unresolved_groups.csv"};
        auto countOutputs = [&]()
        {
            return std::count_if(names.begin(), names.end(), [&](const auto &name)
            {
                const auto path = outputDir / name;
                if (!std::filesystem::exists(path))
                    return false;
                if (!std::filesystem::is_regular_file(path))
                    throw std::runtime_error("Superspectrum output paths must be regular files");
                return true;
            });
        };
        if (noOverwrite)
        {
            const auto existing = countOutputs();
            if (static_cast<size_t>(existing) == names.size())
                return;
            if (existing)
                throw std::runtime_error("Superspectrum output is incomplete; rerun without --skip-existing");
        }
        if ((fixed && !(fixed->gridStep > 0 && std::isfinite(fixed->gridStep) && fixed->maxGap > 0 &&
              std::isfinite(fixed->maxGap) && fixed->maxHalfWidth >= 2 &&
              fixed->maxHalfWidth <= maxscale_global)) || retentionTimes.size() != indices.size())
            throw std::invalid_argument("Invalid superspectrum settings or scan metadata");
        const auto grid = fixed ? GridSettings{*fixed, {}} : estimateGrid(file, indices, silent);

        auto staging = outputDir / names.front();
        staging += ".partial";
        if (!std::filesystem::create_directory(staging))
            throw std::runtime_error("Superspectrum partial output already exists");
        std::array<bool, names.size()> published{}, backedUp{};
        try
        {
            std::ofstream csv(staging / names[0]);
            std::ofstream trackCsv(staging / names[1]);
            if (!csv || !trackCsv)
                throw std::runtime_error("Could not open superspectrum output");
            csv << "firstScan,lastScan,centroid,scan,retentionTimeSeconds,mz,height,gridCenterMz,halfWidth,b1,b2,b3,logRise,gridStepMz,track\n";
            csv << std::setprecision(std::numeric_limits<double>::max_digits10);
            trackCsv << std::setprecision(std::numeric_limits<double>::max_digits10);
            TrackChromatograms chromatograms;
            SuperspectrumTracks tracks(trackCsv, scansPerWindow,
                [&](size_t id, double mz, std::optional<double> spread, const auto &heights)
                {
                    chromatograms.append(indices, retentionTimes, id, mz, spread, heights);
                });
            std::deque<ProfileScan> scans;
            for (size_t i = 0; grid.fallback.gridStep > 0 && i < indices.size(); ++i)
            {
                auto scan = readProfile(file, indices[i], retentionTimes[i]);
                findProfileRuns(scan, grid);
                scans.push_back(std::move(scan));
                if (scans.size() == scansPerWindow)
                {
                    fitWindow(csv, tracks, scans, grid);
                    scans.pop_front();
                }
            }
            tracks.finish();
            csv.close();
            trackCsv.close();
            if (!csv || !trackCsv)
                throw std::runtime_error("Could not complete superspectrum output");
            std::vector<FeaturePeak> features;
            findFeatures(&chromatograms.bins, &retentionTimes, &features, &chromatograms.masses);
            if (!printFeatureList(&features, staging, stem, &chromatograms.bins, true, false, nullptr,
                    &chromatograms.trackIds))
                throw std::runtime_error("Could not write superspectrum feature output");
            if (countOutputs() && noOverwrite)
                throw std::runtime_error("Superspectrum output appeared while processing");
            for (size_t i = 0; i < names.size(); ++i)
            {
                const auto target = outputDir / names[i];
                if (std::filesystem::exists(target))
                {
                    std::filesystem::rename(target, staging / ("previous-" + names[i]));
                    backedUp[i] = true;
                }
                std::filesystem::rename(staging / names[i], target);
                published[i] = true;
            }
            if (!silent)
                printf("    Superspectrum feature detection: %zu chromatograms, %zu features\n",
                    chromatograms.bins.size(), features.size());
        }
        catch (...)
        {
            std::error_code error;
            for (size_t i = 0; i < names.size(); ++i)
            {
                if (published[i])
                    std::filesystem::remove(outputDir / names[i], error);
                if (backedUp[i])
                    std::filesystem::rename(staging / ("previous-" + names[i]), outputDir / names[i], error);
                std::filesystem::remove(staging / names[i], error);
            }
            std::filesystem::remove(staging, error);
            throw;
        }
        std::error_code error;
        for (const auto &name : names)
            std::filesystem::remove(staging / ("previous-" + name), error);
        std::filesystem::remove(staging, error);
    }
}
