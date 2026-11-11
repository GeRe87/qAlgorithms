#include "qalgorithms_intensity_segments.h"

#include <cmath>
#include <stdexcept>

namespace qAlgorithms
{
    IntensitySegments segmentIntensities(const std::vector<AlignedIntensityTrace> &traces,
        double coordinateStep, double relativeTolerance, double threshold, size_t minimumLength)
    {
        if (traces.empty() || !std::isfinite(coordinateStep) || coordinateStep <= 0
            || !(relativeTolerance >= 0 && relativeTolerance < 1)
            || !std::isfinite(threshold) || threshold < 0 || minimumLength == 0)
        {
            throw std::invalid_argument("Invalid intensity segmentation request");
        }
        const size_t length = traces.front().length;
        for (const auto &trace : traces)
        {
            if (trace.length != length || (length && (!trace.intensity || !trace.coordinates)))
            {
                throw std::invalid_argument("Trace lengths or arrays do not match");
            }
        }

        IntensitySegments result;
        const auto *coordinates = traces.front().coordinates;
        size_t start = 0;
        auto previous = IntensitySpanStatus::eligible;
        auto append = [&](size_t end)
        {
            if (end > start)
            {
                const auto status = previous == IntensitySpanStatus::eligible && end - start < minimumLength
                    ? IntensitySpanStatus::tooShort : previous;
                result.spans.push_back({start, end, status});
            }
        };
        for (size_t i = 0; i < length; ++i)
        {
            const double x = coordinates[i];
            if (!std::isfinite(x) || (i && !(x > coordinates[i - 1])))
            {
                throw std::invalid_argument("Expected finite increasing coordinates");
            }
            auto status = IntensitySpanStatus::eligible;
            for (const auto &trace : traces)
            {
                if (trace.coordinates[i] != x)
                {
                    throw std::invalid_argument("Trace coordinates are not aligned");
                }
                const double value = trace.intensity[i];
                const auto sample = !std::isfinite(value) ? IntensitySpanStatus::nonfinite
                    : value < 0 ? IntensitySpanStatus::negative
                    : value <= threshold ? IntensitySpanStatus::belowThreshold : IntensitySpanStatus::eligible;
                if (sample > status)
                {
                    status = sample;
                }
            }
            const bool coordinateBreak = i && std::abs((x - coordinates[i - 1]) / coordinateStep - 1) > relativeTolerance;
            if (coordinateBreak)
            {
                result.coordinateBreaks.push_back(i);
            }
            if (status != previous || (coordinateBreak && status == IntensitySpanStatus::eligible))
            {
                append(i);
                start = i;
                previous = status;
            }
        }
        append(length);
        return result;
    }
}
