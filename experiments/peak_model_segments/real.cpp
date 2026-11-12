#include "qalgorithms_peak_search.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <utility>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

static void compareModels(size_t scan, const std::vector<float> &raw,
    const SegmentedHyperbolicSearchResult &hyperbolic, std::ostream &spans, std::ostream &peaks)
{
    std::vector<float> sampleIndices(raw.size());
    std::iota(sampleIndices.begin(), sampleIndices.end(), 0.0f);
    std::vector<std::pair<double, size_t>> apices;
    for (size_t index : hyperbolic.diagnosticIndices)
    {
        const auto &entry = hyperbolic.search.windows[index];
        const double apex = static_cast<double>(entry.window.center) + entry.inference.properties.front().apexU.value
            / hyperbolic.search.fixedStep;
        apices.emplace_back(apex, index);
    }
    std::sort(apices.begin(), apices.end());
    for (const auto &span : hyperbolic.segmentation.spans)
    {
        if (span.status != IntensitySpanStatus::eligible)
        {
            continue;
        }
        const auto first = std::lower_bound(apices.begin(), apices.end(),
            std::pair<double, size_t>{static_cast<double>(span.start), 0});
        const auto last = std::lower_bound(first, apices.end(),
            std::pair<double, size_t>{static_cast<double>(span.end), 0});
        const auto groups = std::count_if(hyperbolic.diagnosticGroups.begin(), hyperbolic.diagnosticGroups.end(),
            [&](const auto &group) { return group.start >= span.start && group.end < span.end; });
        spans << scan << ',' << span.start << ',' << span.end << ',' << last - first << ',' << groups << ',';
        PeakSearchRequest request;
        request.intensity = raw.data() + span.start;
        request.coordinates = sampleIndices.data() + span.start;
        request.length = span.end - span.start;
        try
        {
            const auto quadratic = std::get<PiecewiseQuadraticSearchResult>(findPeaks(request));
            spans << quadratic.status << ',' << quadratic.peaks.size() << ",\n";
            for (const auto &peak : quadratic.peaks)
            {
                peaks << scan << ',' << span.start << ',' << span.end << ',' << peak.position << ',';
                if (first != last)
                {
                    const auto nearest = std::min_element(first, last, [&](const auto &a, const auto &b)
                    {
                        return std::abs(a.first - peak.position) < std::abs(b.first - peak.position);
                    });
                    const auto &window = hyperbolic.search.windows[nearest->second].window;
                    peaks << nearest->first << ',' << nearest->first - peak.position << ',' << window.center << ',' << window.scale;
                }
                else
                {
                    peaks << ",,,";
                }
                peaks << '\n';
            }
        }
        catch (const std::runtime_error &error)
        {
            // The legacy selector can reject a whole span when its candidate
            // constraints cannot be satisfied. Keep it in the denominator.
            spans << ",,\"" << error.what() << "\"\n";
        }
    }
}

int main(int argc, char **argv)
{
    try
    {
        const bool supportedOnly = argc == 8 && std::string_view(argv[7]) == "--supported-only";
        if (argc != 4 && argc != 6 && argc != 7 && !supportedOnly)
        {
            throw std::invalid_argument(
                "Usage: real_peak_segments INPUT.mzML SCANS.csv WINDOWS.csv [SPANS.csv PEAKS.csv [APICES.csv [--supported-only]]]");
        }
        const std::filesystem::path input(argv[1]);
        XML_File file(input.c_str(), mzML);
        const auto scans = filter_spectra(&file, true, Polarities::positive, true);
        if (file.defective || scans.empty())
        {
            throw std::runtime_error("Missing profile MS1 data");
        }
        std::ofstream summary(argv[2]), windows(argv[3]);
        std::ofstream comparisonSpans, comparisonPeaks, apexGroups;
        if (argc >= 6)
        {
            comparisonSpans.open(argv[4]);
            comparisonPeaks.open(argv[5]);
            comparisonSpans << "scan,start,end,hyperbolic_windows,hyperbolic_groups,quadratic_status,quadratic_peaks,error\n";
            comparisonPeaks << "scan,start,end,quadratic_apex,hyperbolic_apex,delta,center,scale\n";
            comparisonPeaks.precision(17);
        }
        if (argc >= 7)
        {
            apexGroups.open(argv[6]);
            apexGroups << "scan,start,end,maximum_first,maximum_last,windows\n";
        }
        summary << "scan,points,zeros,positive_below_threshold,segments,eligible_segments,short_segments,"
                   "eligible_points,coordinate_breaks,declared,fitted,available,finite_coefficients,integrable,"
                   "segment_us,logs_us,running_us,search_us\n";
        windows << "scan,center,scale,available,integrable\n";
        constexpr double threshold = 2.7182818284590452;
        constexpr double step = .35;
        constexpr size_t maxScale = 8;
        const HyperbolicModel::Geometry geometry(maxScale, step);
        const size_t count = std::min<size_t>(80, scans.size());
        using Clock = std::chrono::steady_clock;
        auto microseconds = [](auto elapsed) { return std::chrono::duration<double, std::micro>(elapsed).count(); };
        for (size_t k = 0; k < count; ++k)
        {
            const size_t scan = scans[(2 * k + 1) * scans.size() / (2 * count)];
            std::vector<float> mz, raw;
            if (get_spectrum(&file, &mz, &raw, scan) != 0 || mz.size() != raw.size() || raw.size() < MINLENGTH)
            {
                throw std::runtime_error("Invalid selected profile");
            }
            std::vector<double> coordinates, gaps;
            for (float value : mz)
            {
                coordinates.push_back(std::sqrt(static_cast<double>(value)));
                if (coordinates.size() > 1)
                {
                    gaps.push_back(coordinates.back() - coordinates[coordinates.size() - 2]);
                }
            }
            const auto middle = gaps.begin() + gaps.size() / 3;
            std::nth_element(gaps.begin(), middle, gaps.end());
            const double coordinateStep = *middle;
            const std::vector<AlignedIntensityTrace> traces = {{raw.data(), coordinates.data(), raw.size()}};
            const auto plan = makeUniformHyperbolicSearchPlan(raw.size(), step, 2, std::min(maxScale, (raw.size() - 1) / 2));

            const auto start = Clock::now();
            const auto segments = segmentIntensities(traces, coordinateStep, .45, threshold, MINLENGTH);
            const auto afterSegments = Clock::now();
            size_t eligible = 0, shortSegments = 0, eligiblePoints = 0;
            std::vector<std::vector<float>> logs;
            for (const auto &span : segments.spans)
            {
                if (span.status == IntensitySpanStatus::tooShort)
                {
                    ++shortSegments;
                }
                if (span.status == IntensitySpanStatus::eligible)
                {
                    ++eligible;
                    eligiblePoints += span.end - span.start;
                    auto &values = logs.emplace_back();
                    values.reserve(span.end - span.start);
                    for (size_t i = span.start; i < span.end; ++i)
                    {
                        values.push_back(static_cast<float>(std::log(static_cast<double>(raw[i]))));
                    }
                }
            }
            const auto afterLogs = Clock::now();
            size_t finite = 0, integrable = 0;
            for (const auto &values : logs)
            {
                HyperbolicModel::runningFit({values.data()}, values.size(), geometry, [&](const auto &fit)
                {
                    const auto &p = fit.locals.front();
                    const bool valid = fit.identifiable && std::isfinite(fit.sharedSlope)
                        && std::isfinite(p.intercept) && std::isfinite(p.curvature) && std::isfinite(p.skew);
                    finite += valid;
                    integrable += valid && detail::hyperbolic::validIntegrablePeak(p);
                });
            }
            const auto afterRunning = Clock::now();
            PeakSearchRequest request;
            request.model = PeakModelKind::Hyperbolic;
            request.intensity = raw.data();
            request.length = raw.size();
            request.hyperbolicPlan = &plan;
            const PeakSearchRequest::Segmentation segmentation{coordinates.data(), coordinateStep, .45, threshold,
                supportedOnly ? SegmentedInferenceDetail::supportedWindows : SegmentedInferenceDetail::allWindows};
            request.segmentation = &segmentation;
            const auto result = std::get<SegmentedHyperbolicSearchResult>(findPeaks(request));
            const auto afterSearch = Clock::now();
            if (argc >= 6)
            {
                compareModels(scan, raw, result, comparisonSpans, comparisonPeaks);
            }
            if (argc >= 7)
            {
                for (size_t i = 0; i < result.apexGroups.size(); ++i)
                {
                    const auto &bounds = result.diagnosticGroups[i];
                    const auto &group = result.apexGroups[i];
                    for (const auto &maximum : group.maxima)
                    {
                        apexGroups << scan << ',' << bounds.start << ',' << bounds.end << ',' << maximum.first << ','
                                   << maximum.last << ',' << maximum.windowIndices.size() << '\n';
                    }
                    apexGroups << scan << ',' << bounds.start << ',' << bounds.end << ",,,"
                               << group.unassignedWindowIndices.size() << '\n';
                }
            }
            size_t fitted = 0, available = 0;
            for (size_t i = 0; i < result.search.windows.size(); ++i)
            {
                if (result.windowStatus[i] != SegmentedWindowStatus::evaluated)
                {
                    continue;
                }
                ++fitted;
                const auto &entry = result.search.windows[i];
                available += entry.inference.status == InferenceStatus::available;
                const bool shape = !entry.inference.properties.empty() && entry.inference.properties.front().integrable;
                if (entry.window.scale == 4 || entry.window.scale == 8)
                {
                    windows << scan << ',' << entry.window.center << ',' << entry.window.scale << ','
                            << (entry.inference.status == InferenceStatus::available) << ',' << shape << '\n';
                }
            }
            const auto zeros = std::count(raw.begin(), raw.end(), 0.0f);
            const auto weak = std::count_if(raw.begin(), raw.end(), [](float value) { return value > 0 && value <= threshold; });
            summary << scan << ',' << raw.size() << ',' << zeros << ',' << weak << ',' << eligible + shortSegments
                    << ',' << eligible << ',' << shortSegments << ',' << eligiblePoints << ','
                    << segments.coordinateBreaks.size() << ',' << plan.familySize() << ',' << fitted << ','
                    << available << ',' << finite << ',' << integrable << ','
                    << microseconds(afterSegments - start) << ',' << microseconds(afterLogs - afterSegments) << ','
                    << microseconds(afterRunning - afterLogs) << ',' << microseconds(afterSearch - afterRunning) << '\n';
        }
        summary.close();
        windows.close();
        if (argc >= 6)
        {
            comparisonSpans.close();
            comparisonPeaks.close();
        }
        if (argc >= 7)
        {
            apexGroups.close();
        }
        if (!summary || !windows || !comparisonSpans || !comparisonPeaks || !apexGroups)
        {
            throw std::runtime_error("Incomplete segment evidence");
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
