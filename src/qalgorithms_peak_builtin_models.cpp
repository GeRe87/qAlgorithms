#include "detail/qalgorithms_peak_background.h"
#include "qalgorithms_peak_validation.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <type_traits>

// Built-in single-window adapters for the public enum entry points.
// This file connects model fits to the validation and background routines.
namespace qAlgorithms
{
    namespace peak_model
    {
        namespace
        {
            // Replay the log-domain shape fit on every rotation. The shape's
            // intercept cancels because raw amplitude is estimated separately.
            class PeakTemplate
            {
            public:
                PeakTemplate(PeakModelKind model, size_t length, double step, size_t minimum)
                    : model_(model), center_(length / 2), minimum_(minimum), step_(step), logs_(length)
                {
                    if (model_ == PeakModelKind::Hyperbolic)
                    {
                        try
                        {
                            geometry_.emplace(center_, step_);
                        }
                        catch (const std::invalid_argument &)
                        {
                            // Unsupported or degenerate hyperbolic geometry.
                        }
                    }
                }

                bool identifiable() const
                {
                    return model_ == PeakModelKind::PiecewiseQuadratic
                        || (geometry_ && geometry_->at(center_).identifiable);
                }

                template<typename Visit>
                bool fitWindow(const std::vector<double> &values, Visit visit)
                {
                    for (size_t i = 0; i < values.size(); ++i)
                    {
                        if (!(values[i] > 0 && std::isfinite(values[i])))
                        {
                            return false;
                        }
                        logs_[i] = static_cast<float>(std::log(values[i]));
                    }
                    if (model_ == PeakModelKind::PiecewiseQuadratic)
                    {
                        quadratic_.clear();
                        Model<PeakModelKind::PiecewiseQuadratic>::runningFit(logs_.data(), logs_.size(), &quadratic_);
                        visit(quadratic_.back()); // full central window, last in scale order
                    }
                    else
                    {
                        HyperbolicModel::runningFit({logs_.data()}, logs_.size(), *geometry_, [&](const auto &fit)
                        {
                            if (fit.center == center_ && fit.scale == center_)
                            {
                                visit(fit);
                            }
                        });
                    }
                    return true;
                }

                bool evaluate(const std::vector<double> &values, const std::vector<double> &x, std::vector<double> &shape)
                {
                    bool supported = false;
                    fitWindow(values, [&](const auto &fit) { supported = makeShape(fit, x, shape); });
                    return supported;
                }

                bool makeShape(const RegCoeffs &fit, const std::vector<double> &x, std::vector<double> &shape) const
                {
                    const auto check = validatePeakGeometry(fit, minimum_);
                    // The template uses the entire tested window. A
                    // mixed flank must not turn upward within that window.
                    const double edge = static_cast<double>(center_);
                    if (check.apexSupport != PeakValidationStatus::passed
                        || (fit.b2 > 0 && fit.b1 < 2 * edge * fit.b2)
                        || (fit.b3 > 0 && -fit.b1 < 2 * edge * fit.b3))
                    {
                        return false;
                    }
                    for (size_t i = 0; i < x.size(); ++i)
                    {
                        const double u = x[i];
                        shape[i] = (fit.b1 + u * (u < 0 ? fit.b2 : fit.b3)) * u;
                    }
                    return normalize(shape);
                }

                bool makeShape(const HyperbolicModel::Fit &fit, const std::vector<double> &x, std::vector<double> &shape) const
                {
                    auto coefficients = fit.locals.front();
                    const auto check = validatePeakGeometry(coefficients, fit.sharedSlope, center_, step_, minimum_);
                    if (!fit.identifiable || check.apexSupport != PeakValidationStatus::passed)
                    {
                        return false;
                    }
                    coefficients.intercept = 0;
                    for (size_t i = 0; i < x.size(); ++i)
                    {
                        shape[i] = detail::hyperbolic::predict(coefficients, fit.sharedSlope, x[i] * step_);
                    }
                    return normalize(shape);
                }

                const HyperbolicModel::Geometry &geometry() const { return *geometry_; }
                const std::vector<float> &logs() const { return logs_; }

            private:
                static bool normalize(std::vector<double> &shape)
                {
                    if (!std::all_of(shape.begin(), shape.end(), [](double value) { return std::isfinite(value); }))
                    {
                        return false;
                    }
                    const double maximum = *std::max_element(shape.begin(), shape.end());
                    for (double &value : shape)
                    {
                        value = std::exp(value - maximum);
                    }
                    return true;
                }

                PeakModelKind model_;
                size_t center_, minimum_;
                double step_;
                std::vector<float> logs_;
                std::vector<RegCoeffs> quadratic_;
                std::optional<HyperbolicModel::Geometry> geometry_;
            };

            void validateRequest(PeakModelKind model, const float *intensities,
                size_t length, double fixedStep, const PeakBackgroundOptions &options)
            {
                if (intensities == nullptr || length < 5 || length % 2 == 0 || options.rotations == 0
                    || !(options.alpha > 0 && options.alpha < 1) || !(fixedStep > 0 && std::isfinite(fixedStep))
                    || (model != PeakModelKind::PiecewiseQuadratic && model != PeakModelKind::Hyperbolic)
                    || (model == PeakModelKind::PiecewiseQuadratic && (fixedStep != 1 || length / 2 > maxscale_global)))
                {
                    throw std::invalid_argument("Invalid raw peak validation request");
                }
                if (!std::all_of(intensities, intensities + length, [](float value) { return std::isfinite(value); }))
                {
                    throw std::invalid_argument("Raw peak validation requires finite intensities");
                }
            }

            std::vector<double> centeredGrid(size_t length)
            {
                std::vector<double> grid(length);
                std::iota(grid.begin(), grid.end(), -static_cast<double>(length / 2));
                return grid;
            }

            template<typename Predict>
            PeakFormNoise estimateLogNoise(const std::vector<float> &logs, const std::vector<double> &x,
                double alpha, Predict predict)
            {
                double rss = 0;
                for (size_t i = 0; i < logs.size(); ++i)
                {
                    const double residual = logs[i] - predict(x[i]);
                    rss += residual * residual;
                }
                const size_t df = logs.size() - 4; // both single-trace fits have four parameters
                return {rss / static_cast<double>(df), df, alpha};
            }
        }

        PeakBackgroundValidation validatePeakBackground(PeakModelKind model, const float *intensities,
            size_t length, double fixedStep, const PeakBackgroundOptions &options)
        {
            validateRequest(model, intensities, length, fixedStep, options);
            PeakTemplate peak(model, length, fixedStep, options.minimumSamplesPerSide);
            if (!peak.identifiable())
            {
                return {};
            }
            const std::vector<double> values(intensities, intensities + length), x = centeredGrid(length);
            std::vector<double> shape(length);
            const bool supported = peak.evaluate(values, x, shape);
            return runPeakBackground(values, x, shape, supported, [&](const auto &draw, auto &templateValues)
            {
                return peak.evaluate(draw, x, templateValues);
            }, options);
        }

        FullPeakCandidateValidation validatePeakCandidate(PeakModelKind model, const float *intensities,
            size_t length, double fixedStep, const FullPeakCandidateOptions &options)
        {
            validateRequest(model, intensities, length, fixedStep, options.background);
            if (!(options.minimumArea > 0 && std::isfinite(options.minimumArea))
                || !(options.minimumCoverage > 0 && options.minimumCoverage < 1)
                || !(options.logAlpha > 0 && options.logAlpha < .5))
            {
                throw std::invalid_argument("Invalid area, coverage or alpha for full peak validation");
            }
            FullPeakCandidateValidation result;
            result.stages.geometry.form = PeakValidationStatus::unavailable;
            const size_t minimum = options.background.minimumSamplesPerSide;
            PeakTemplate peak(model, length, fixedStep, minimum);
            if (!peak.identifiable())
            {
                return result;
            }
            const std::vector<double> values(intensities, intensities + length), x = centeredGrid(length);
            std::vector<double> shape;
            bool supported = false;
            peak.fitWindow(values, [&](const auto &fit)
            {
                using Fit = std::decay_t<decltype(fit)>;
                if constexpr (std::is_same_v<Fit, RegCoeffs>)
                {
                    using Quadratic = Model<PeakModelKind::PiecewiseQuadratic>;
                    result.stages.geometry = validatePeakGeometry(fit, minimum);
                    if (result.stages.geometry.apexSupport == PeakValidationStatus::passed)
                    {
                        result.logNoise = estimateLogNoise(peak.logs(), x, options.logAlpha,
                            [&](double sample) { return Quadratic::logValue(&fit, sample); });
                        result.stages = validatePeakCandidate(fit, result.logNoise, minimum,
                            options.minimumCoverage, options.minimumArea);
                    }
                }
                else
                {
                    const auto &c = fit.locals.front();
                    result.stages.geometry = validatePeakGeometry(c, fit.sharedSlope, fit.scale, fixedStep, minimum);
                    if (result.stages.geometry.apexSupport == PeakValidationStatus::passed)
                    {
                        result.logNoise = estimateLogNoise(peak.logs(), x, options.logAlpha, [&](double sample)
                        {
                            return detail::hyperbolic::predict(c, fit.sharedSlope, sample * fixedStep);
                        });
                        result.stages = validatePeakCandidate(fit, 0, peak.geometry(), result.logNoise, minimum,
                            options.minimumCoverage, options.minimumArea);
                    }
                }
                if (result.stages.windowCoverage.status == PeakValidationStatus::passed)
                {
                    shape.resize(length);
                    supported = peak.makeShape(fit, x, shape);
                }
            });
            for (PeakValidationStatus stage : {result.stages.geometry.form, result.stages.geometry.apexSupport,
                result.stages.formConfidence, result.stages.areaEvidence.status, result.stages.windowCoverage.status})
            {
                if (stage != PeakValidationStatus::passed)
                {
                    result.status = stage;
                    return result;
                }
            }
            result.background = runPeakBackground(values, x, shape, supported, [&](const auto &draw, auto &templateValues)
            {
                return peak.evaluate(draw, x, templateValues);
            }, options.background);
            result.status = result.background.status;
            return result;
        }
    }
}
