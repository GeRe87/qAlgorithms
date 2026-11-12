#ifndef QALGORITHMS_DETAIL_PEAK_INFERENCE_H
#define QALGORITHMS_DETAIL_PEAK_INFERENCE_H

#include "qalgorithms_peak_inference.h"

namespace qAlgorithms
{
    namespace peak_model
    {
        // Internal path for a kernel fit and matching, already validated traces.
        // The fit and both pointer arrays use the same trace-relative indices.
        HyperbolicInference inferFittedHyperbolicWindow(
            const std::vector<const float *> &rawTraces, const std::vector<const float *> &logs,
            size_t length, const HyperbolicModel::Fit &fit, const HyperbolicModel::Geometry &cache,
            const HyperbolicInferenceOptions &options, bool supportedOnly);
    }
}

#endif
