#ifndef QALGORITHMS_INTENSITY_SEGMENTS_H
#define QALGORITHMS_INTENSITY_SEGMENTS_H

#include <cstddef>
#include <vector>

namespace qAlgorithms
{
    struct AlignedIntensityTrace
    {
        const float *intensity;
        const double *coordinates;
        size_t length;
    };

    enum class IntensitySpanStatus { eligible, tooShort, belowThreshold, negative, nonfinite };

    struct IntensitySpan
    {
        size_t start, end; // half-open indices in every original trace
        IntensitySpanStatus status;
    };

    struct IntensitySegments
    {
        std::vector<IntensitySpan> spans;
        std::vector<size_t> coordinateBreaks; // index on the right of each incompatible edge
    };

    // Coordinates must be finite, strictly increasing and identical across traces.
    // Split at incompatible coordinate steps and intersect intensity eligibility.
    IntensitySegments segmentIntensities(const std::vector<AlignedIntensityTrace> &traces,
        double coordinateStep, double relativeTolerance, double threshold, size_t minimumLength);
}

#endif
