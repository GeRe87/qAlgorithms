#include "detail/qalgorithms_single_trace_hyperbolic.h"
#include "detail/qalgorithms_peak_inference.h"
#include "qalgorithms_utils.h"

#include <cmath>
#include <limits>
#include <tuple>

namespace qAlgorithms::detail
{
    namespace
    {
        constexpr double stepU = hyperbolicSingleTraceStep;

        const hyperbolic::GeometryCache &geometry()
        {
            static const hyperbolic::GeometryCache cache(
                std::min(maxscale_global, static_cast<size_t>(30 / stepU)), stepU);
            return cache;
        }

        hyperbolic::LocalFit localFit(const RegCoeffs &c)
        {
            return {c.b0, c.b2, c.b3};
        }
    }

    void SingleTraceHyperbolic::runningFit(const float *logs, size_t length, std::vector<Coefficients> *output)
    {
        hyperbolic::forEachRunning({logs}, length, geometry(), [&](const hyperbolic::FitView &fit)
        {
            const auto &local = fit.locals.front();
            if (fit.identifiable && std::isfinite(fit.sharedSlope) && std::isfinite(local.intercept) &&
                hyperbolic::validIntegrablePeak(local))
            {
                output->push_back({{local.intercept, fit.sharedSlope, local.curvature, local.skew,
                    static_cast<uint16_t>(fit.scale), static_cast<uint16_t>(fit.center)}, local.intercept});
            }
        });
        // The shared conflict resolver examines smaller windows first.
        std::sort(output->begin(), output->end(), [](const auto &a, const auto &b)
        {
            return std::tie(a.values.scale, a.values.x0) < std::tie(b.values.scale, b.values.x0);
        });
    }

    Span_i32 SingleTraceHyperbolic::regressionSpan(const Coefficients *c)
    {
        return {int32_t(center(*c) - halfWidth(*c)), int32_t(2 * halfWidth(*c) + 1)};
    }

    double SingleTraceHyperbolic::logValue(const Coefficients *c, double offset)
    {
        return hyperbolic::predict(localFit(c->values), c->values.b1, stepU * offset);
    }

    double SingleTraceHyperbolic::apex(const Coefficients *c)
    {
        return hyperbolic::apexU(localFit(c->values), c->values.b1) / stepU;
    }

    bool SingleTraceHyperbolic::apexHasSupport(float position, const Coefficients *c)
    {
        return std::abs(position) < double(halfWidth(*c)) - 1;
    }

    bool SingleTraceHyperbolic::insignificantCurvature(const Coefficients *c, double mse, size_t df)
    {
        const auto &g = geometry().at(halfWidth(*c));
        const double lowerCurvature = c->values.b2 - std::abs(c->values.b3);
        const double se = std::sqrt(mse * (g.varCurvatureFactor() + g.varLocalSkewFactor(1)));
        return !(lowerCurvature > se * t_critical_95(df));
    }

    double SingleTraceHyperbolic::refitAmplitude(const float *raw, Span_i32 span, float *prediction, Coefficients *c)
    {
        const auto fitted = hyperbolic::refitRawAmplitude(raw, size_t(span.endIdx()) + 1,
            center(*c), geometry().at(halfWidth(*c)), stepU, localFit(c->values), c->values.b1, 1);
        c->values.b0 = fitted.logAmplitude;
        for (int32_t i = span.startIdx; i <= span.endIdx(); ++i)
        {
            const double value = std::exp(logValue(c, double(i) - center(*c)));
            prediction[i] = static_cast<float>(std::min(value, double(std::numeric_limits<float>::max())));
        }
        return fitted.logAmplitude;
    }

    bool SingleTraceHyperbolic::measure(Result &result, const float *raw, const float *logs,
        size_t length, double, double step)
    {
        const auto &c = result.coeffs;
        const hyperbolic::FitView fit{center(c), halfWidth(c), c.values.b1, true,
            {{c.logIntercept, c.values.b2, c.values.b3}}};
        const auto inference = peak_model::inferFittedHyperbolicWindow({raw}, {logs}, length,
            fit, geometry(), {}, true);
        if (inference.status != peak_model::InferenceStatus::available)
            return false;
        const auto &properties = inference.properties.front();
        const double unitsPerU = step / stepU;
        const double area = properties.areaU.value * unitsPerU;
        const double areaUncertainty = std::sqrt(inference.covariance[8]) * unitsPerU;
        const double positionUncertainty = std::sqrt(inference.covariance[0]) / stepU;
        const double heightUncertainty = std::sqrt(inference.covariance[4]);
        constexpr double smallest = std::numeric_limits<float>::denorm_min();
        constexpr double largest = std::numeric_limits<float>::max();
        const double height = properties.height.value;
        const double fwhm = width(&c, height, step);
        if (!(height >= smallest && height <= largest && area >= smallest && area <= largest &&
              fwhm >= smallest && fwhm <= largest && areaUncertainty <= largest &&
              positionUncertainty <= largest && positionUncertainty * step <= largest &&
              heightUncertainty <= largest))
            return false;
        result.area = static_cast<float>(area);
        result.area_unc = static_cast<float>(areaUncertainty);
        result.position_unc = static_cast<float>(positionUncertainty);
        result.height_unc = static_cast<float>(heightUncertainty);
        return true;
    }

    double SingleTraceHyperbolic::width(const Coefficients *c, double, double step)
    {
        const double peak = apex(c);
        const double halfHeight = logValue(c, peak) - std::log(2.0);
        auto flank = [&](double direction)
        {
            double lower = 0, upper = 1;
            while (logValue(c, peak + direction * upper) > halfHeight && upper < 4096)
                upper *= 2;
            if (!(logValue(c, peak + direction * upper) <= halfHeight))
                return std::numeric_limits<double>::quiet_NaN();
            for (int i = 0; i < 45; ++i)
            {
                const double mid = (lower + upper) / 2;
                if (logValue(c, peak + direction * mid) > halfHeight)
                    lower = mid;
                else
                    upper = mid;
            }
            return (lower + upper) / 2;
        };
        return (flank(-1) + flank(1)) * step;
    }
}
