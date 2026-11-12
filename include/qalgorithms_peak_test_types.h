#ifndef QALGORITHMS_PEAK_TEST_TYPES_H
#define QALGORITHMS_PEAK_TEST_TYPES_H

#include <cstddef>
#include <cstdint>
#include <limits>

namespace qAlgorithms
{
    namespace peak_model
    {
        enum class PeakValidationStatus { notEvaluated, failed, passed, unavailable };

        struct PeakBackgroundOptions
        {
            double alpha = .05;
            size_t rotations = 999;
            uint64_t seed = 0;
            size_t minimumSamplesPerSide = 2;
        };

        struct PeakBackgroundComparison
        {
            double backgroundRSS = std::numeric_limits<double>::quiet_NaN();
            double explainedFraction = std::numeric_limits<double>::quiet_NaN();
            double pValue = std::numeric_limits<double>::quiet_NaN();
        };

        struct PeakBackgroundValidation
        {
            PeakValidationStatus status = PeakValidationStatus::unavailable;
            bool supportedPeakShape = false;
            PeakBackgroundComparison constant, linear;
        };
    }
}

#endif
