#ifndef QALGORITHMS_DETAIL_PEAK_PROPERTIES_H
#define QALGORITHMS_DETAIL_PEAK_PROPERTIES_H

#include "qalgorithms_peak_validation.h"
#include <algorithm>
#include <cmath>

namespace qAlgorithms
{
    namespace peak_model
    {
        // Fixed-window properties of the original log-OLS fit. A numerical
        // failure or nonsmooth shape leaves the joint estimate unavailable.
        std::optional<CandidateFitOutput::Properties> describeLogFit(const RegCoeffs &fit, double mse);
        std::optional<CandidateFitOutput::Properties> describeLogFit(const HyperbolicModel::Fit &fit,
            const HyperbolicModel::Geometry &geometry, double mse);

        inline bool validPositionAreaCovariance(const std::array<double, 4> &covariance)
        {
            return std::all_of(covariance.begin(), covariance.end(), [](double value) { return std::isfinite(value); })
                && covariance[0] >= 0 && covariance[3] >= 0 && covariance[1] == covariance[2]
                && std::abs(covariance[1]) <= std::sqrt(covariance[0]) * std::sqrt(covariance[3]);
        }
    }
}

#endif
