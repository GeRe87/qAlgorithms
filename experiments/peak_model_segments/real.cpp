#include "qalgorithms_peak_search.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
        {
            throw std::invalid_argument("Usage: real_peak_segments INPUT.mzML SCANS.csv WINDOWS.csv");
        }
        const std::filesystem::path input(argv[1]);
        XML_File file(input.c_str(), mzML);
        const auto scans = filter_spectra(&file, true, Polarities::positive, true);
        if (file.defective || scans.empty())
        {
            throw std::runtime_error("Missing profile MS1 data");
        }
        std::ofstream summary(argv[2]), windows(argv[3]);
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
            const auto result = plan.evaluateSegmented(traces, coordinateStep, .45);
            const auto afterSearch = Clock::now();
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
        if (!summary || !windows)
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
