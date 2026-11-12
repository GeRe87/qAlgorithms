#ifndef QALGORITHMS_DETAIL_PEAK_BACKGROUND_H
#define QALGORITHMS_DETAIL_PEAK_BACKGROUND_H

#include "qalgorithms_peak_test_types.h"

#include <functional>
#include <vector>

namespace qAlgorithms
{
    namespace peak_model
    {
        // Replay is synchronous and uses the same fitting procedure as the
        // observed template. On success it fills one finite, nonnegative
        // shape value per raw sample.
        // The grid matches the input length, as does a supported template.
        using PeakBackgroundReplay = std::function<bool(const std::vector<double> &, std::vector<double> &)>;

        PeakBackgroundValidation runPeakBackground(const std::vector<double> &values,
            const std::vector<double> &centeredGrid, const std::vector<double> &observedShape,
            bool supported, const PeakBackgroundReplay &replay, const PeakBackgroundOptions &options);
    }
}

#endif
