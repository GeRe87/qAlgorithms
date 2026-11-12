#include "qalgorithms_peak_model_abi.h"

#include <cmath>
#include <cstdlib>
#include <new>

namespace
{
    struct Gaussian
    {
        double intercept, slope, curvature, edge;
        double interceptFactor, slopeFactor, curvatureFactor, interceptCurvature;
    };

    struct FitLifetime
    {
        bool live = false;
        ~FitLifetime()
        {
            if (live)
            {
                std::abort();
            }
        }
    } lifetime;

    int32_t fitWindow(const double *logs, uint64_t length, double step, void **output)
    {
        if (lifetime.live)
        {
            std::abort();
        }
        const double half = static_cast<double>(length / 2);
        double sum = 0, xx = 0, xxxx = 0, xy = 0, xxy = 0;
        for (uint64_t i = 0; i < length; ++i)
        {
            const double x = (static_cast<double>(i) - half) * step;
            const double squared = x * x;
            sum += logs[i];
            xx += squared;
            xxxx += squared * squared;
            xy += x * logs[i];
            xxy += squared * logs[i];
        }
        const double n = static_cast<double>(length);
        const double curvature = (xxy - xx * sum / n) / (xxxx - xx * xx / n);
        const double determinant = n * xxxx - xx * xx;
        *output = new (std::nothrow) Gaussian{(sum - curvature * xx) / n, xy / xx, curvature, half * step,
            xxxx / determinant, 1 / xx, n / determinant, -xx / determinant};
        lifetime.live = *output != nullptr;
        if (lifetime.live && logs[0] < 0)
        {
            return QALGORITHMS_PEAK_FAILED; // exercises destruction even on fit failure
        }
        return *output ? QALGORITHMS_PEAK_PASSED : QALGORITHMS_PEAK_UNAVAILABLE;
    }

    void destroy(void *fit)
    {
        delete static_cast<Gaussian *>(fit);
        lifetime.live = false;
    }

    int32_t shape(const void *fit)
    {
        return static_cast<const Gaussian *>(fit)->curvature < 0 ? QALGORITHMS_PEAK_PASSED : QALGORITHMS_PEAK_FAILED;
    }

    int32_t apex(const void *fit, double *position)
    {
        const auto &g = *static_cast<const Gaussian *>(fit);
        *position = -g.slope / (2 * g.curvature);
        return QALGORITHMS_PEAK_PASSED;
    }

    int32_t predict(const void *fit, double x, double *value)
    {
        const auto &g = *static_cast<const Gaussian *>(fit);
        *value = g.intercept + (g.slope + g.curvature * x) * x;
        return QALGORITHMS_PEAK_PASSED;
    }

    // Integration fixture only. These margins exercise the callback
    // protocol; they are not a statistical calibration of this model.
    int32_t formConfidence(const void *fit, const QAlgorithmsPeakNoiseV1 *noise)
    {
        const auto &g = *static_cast<const Gaussian *>(fit);
        return -g.curvature * g.edge * g.edge > 5 * std::sqrt(noise->log_mse)
            ? QALGORITHMS_PEAK_PASSED : QALGORITHMS_PEAK_FAILED;
    }

    double windowFraction(const Gaussian &g, double left, double right)
    {
        const double position = -g.slope / (2 * g.curvature);
        const double root = std::sqrt(-g.curvature);
        return .5 * (std::erf((right - position) * root) - std::erf((left - position) * root));
    }

    int32_t area(const void *fit, const QAlgorithmsPeakNoiseV1 *noise, double minimum, QAlgorithmsPeakAreaV1 *out)
    {
        const auto &g = *static_cast<const Gaussian *>(fit);
        const double logTotal = g.intercept - g.slope * g.slope / (4 * g.curvature)
            + .5 * std::log(std::acos(-1.0) / -g.curvature);
        *out = {1, logTotal + std::log(windowFraction(g, -g.edge, g.edge)) - 5 * std::sqrt(noise->log_mse),
            -g.edge, g.edge};
        return out->log_lower_bound > std::log(minimum) ? QALGORITHMS_PEAK_PASSED : QALGORITHMS_PEAK_FAILED;
    }

    int32_t coverage(const void *fit, const QAlgorithmsPeakNoiseV1 *noise,
        double left, double right, double minimum, QAlgorithmsPeakCoverageV1 *out)
    {
        const auto &g = *static_cast<const Gaussian *>(fit);
        const double fraction = windowFraction(g, left, right);
        const double margin = std::sqrt(noise->log_mse);
        *out = {fraction, margin, fraction > 5 * margin ? fraction - 5 * margin : 0};
        return out->lower_confidence_bound >= minimum ? QALGORITHMS_PEAK_PASSED : QALGORITHMS_PEAK_FAILED;
    }
}

#ifdef _WIN32
#define MODULE_EXPORT __declspec(dllexport)
#else
#define MODULE_EXPORT __attribute__((visibility("default")))
#endif

extern "C" MODULE_EXPORT const QAlgorithmsPeakModelV1 *qalgorithms_peak_model_v1()
{
#ifdef INCOMPATIBLE_VERSION
    constexpr uint32_t version = 99;
#else
    constexpr uint32_t version = QALGORITHMS_PEAK_MODEL_ABI_V1;
#endif
    static const QAlgorithmsPeakModelV1 model{version, sizeof(QAlgorithmsPeakModelV1), "external-gaussian-ols", 3,
        fitWindow, destroy, shape, apex, predict, formConfidence, area, coverage};
    return &model;
}

#ifndef VALIDATION_ONLY
extern "C" MODULE_EXPORT int32_t qalgorithms_peak_model_properties_v1(
    const void *fit, const QAlgorithmsPeakNoiseV1 *noise, QAlgorithmsPeakPropertiesV1 *out)
{
    const auto &g = *static_cast<const Gaussian *>(fit);
    apex(fit, &out->apex);
    out->area = std::exp(g.intercept - g.slope * g.slope / (4 * g.curvature))
        * std::sqrt(std::acos(-1.0) / -g.curvature);
    const double positionSlope = -1 / (2 * g.curvature);
    const double positionCurve = g.slope / (2 * g.curvature * g.curvature);
    const double areaSlope = out->area * out->apex;
    const double areaCurve = out->area * (out->apex * out->apex - .5 / g.curvature);
    out->covariance[0] = noise->log_mse * (positionSlope * positionSlope * g.slopeFactor
        + positionCurve * positionCurve * g.curvatureFactor);
    out->covariance[1] = out->covariance[2] = noise->log_mse * (positionSlope * areaSlope * g.slopeFactor
        + positionCurve * (areaCurve * g.curvatureFactor + out->area * g.interceptCurvature));
    out->covariance[3] = noise->log_mse * (out->area * out->area * g.interceptFactor + areaSlope * areaSlope * g.slopeFactor
        + areaCurve * areaCurve * g.curvatureFactor + 2 * out->area * areaCurve * g.interceptCurvature);
    return QALGORITHMS_PEAK_PASSED;
}
#endif
