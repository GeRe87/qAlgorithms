#include "qalgorithms_peak_validation.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

namespace qAlgorithms
{
    namespace peak_model
    {
        namespace
        {
            struct Background
            {
                double offset = 0, slope = 0;
            };

            double dot(const std::vector<double> &a, const std::vector<double> &b)
            {
                return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
            }

            Background fitBackground(const std::vector<double> &values, const std::vector<double> &x, double xNorm)
            {
                return {std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size()),
                    xNorm > 0 ? dot(values, x) / xNorm : 0};
            }

            void subtractBackground(std::vector<double> &values, const std::vector<double> &x, Background fit)
            {
                for (size_t i = 0; i < values.size(); ++i)
                {
                    values[i] = (values[i] - fit.offset) - fit.slope * x[i];
                }
            }

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

                bool evaluate(const std::vector<double> &values, const std::vector<double> &x, std::vector<double> &shape)
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
                        const auto &fit = quadratic_.back(); // full central window, last in scale order
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
                    }
                    else
                    {
                        bool supported = false;
                        HyperbolicModel::runningFit({logs_.data()}, logs_.size(), *geometry_, [&](const auto &fit)
                        {
                            if (fit.center != center_ || fit.scale != center_ || !fit.identifiable)
                            {
                                return;
                            }
                            auto coefficients = fit.locals.front();
                            const auto check = validatePeakGeometry(coefficients, fit.sharedSlope, center_, step_, minimum_);
                            if (check.apexSupport != PeakValidationStatus::passed)
                            {
                                return;
                            }
                            coefficients.intercept = 0;
                            for (size_t i = 0; i < x.size(); ++i)
                            {
                                shape[i] = detail::hyperbolic::predict(coefficients, fit.sharedSlope, x[i] * step_);
                            }
                            supported = true;
                        });
                        if (!supported)
                        {
                            return false;
                        }
                    }
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

            private:
                PeakModelKind model_;
                size_t center_, minimum_;
                double step_;
                std::vector<float> logs_;
                std::vector<RegCoeffs> quadratic_;
                std::optional<HyperbolicModel::Geometry> geometry_;
            };

            // Projecting both vectors onto the background residual space gives
            // the least-squares improvement from adding a nonnegative template.
            double explainedFraction(const std::vector<double> &residual, std::vector<double> &shape,
                const std::vector<double> &x, double xNorm, double rss)
            {
                subtractBackground(shape, x, fitBackground(shape, x, xNorm));
                const double norm = dot(shape, shape), cross = dot(shape, residual);
                if (!(rss > 0 && norm > 0 && cross > 0))
                {
                    return 0;
                }
                const double correlation = cross / std::sqrt(rss) / std::sqrt(norm);
                return std::min(1.0, correlation * correlation);
            }

            PeakBackgroundComparison compareBackground(const std::vector<double> &values, const std::vector<double> &x,
                double xNorm, const std::vector<double> &observedShape, bool supported, PeakTemplate &peak,
                size_t rotations, std::mt19937_64 &engine)
            {
                const auto background = fitBackground(values, x, xNorm);
                std::vector<double> residual(values);
                subtractBackground(residual, x, background);
                const double rss = dot(residual, residual);
                auto shape = observedShape;
                const double observed = supported ? explainedFraction(residual, shape, x, xNorm, rss) : 0;
                if (observed == 0)
                {
                    return {rss, 0, 1}; // all simulated statistics are at least zero
                }
                std::normal_distribution<double> normal;
                std::vector<double> simulated(values.size());
                size_t exceedances = 0;
                for (size_t draw = 0; draw < rotations; ++draw)
                {
                    // Gaussian residual direction is uniform on the sphere
                    // conditional on the background projection and RSS.
                    for (double &value : residual)
                    {
                        value = normal(engine);
                    }
                    subtractBackground(residual, x, fitBackground(residual, x, xNorm));
                    const double norm = dot(residual, residual);
                    if (!(norm > 0 && std::isfinite(norm)))
                    {
                        return {rss, observed};
                    }
                    const double scale = std::sqrt(rss / norm);
                    for (size_t i = 0; i < values.size(); ++i)
                    {
                        simulated[i] = background.offset + background.slope * x[i] + scale * residual[i];
                    }
                    if (peak.evaluate(simulated, x, shape))
                    {
                        const auto replayBackground = fitBackground(simulated, x, xNorm);
                        residual = simulated;
                        subtractBackground(residual, x, replayBackground);
                        exceedances += explainedFraction(residual, shape, x, xNorm, dot(residual, residual)) >= observed;
                    }
                }
                const double pValue = (1 + static_cast<double>(exceedances)) / (1 + static_cast<double>(rotations));
                return {rss, observed, pValue};
            }
        }

        PeakBackgroundValidation validatePeakBackground(PeakModelKind model, const float *intensities,
            size_t length, double fixedStep, const PeakBackgroundOptions &options)
        {
            if (intensities == nullptr || length < 5 || length % 2 == 0 || options.rotations == 0
                || !(options.alpha > 0 && options.alpha < 1) || !(fixedStep > 0 && std::isfinite(fixedStep))
                || (model != PeakModelKind::PiecewiseQuadratic && model != PeakModelKind::Hyperbolic)
                || (model == PeakModelKind::PiecewiseQuadratic && (fixedStep != 1 || length / 2 > maxscale_global)))
            {
                throw std::invalid_argument("Invalid raw peak background test request");
            }
            std::vector<double> values(intensities, intensities + length), x(length), shape(length);
            for (size_t i = 0; i < length; ++i)
            {
                if (!std::isfinite(values[i]))
                {
                    throw std::invalid_argument("Raw peak background test requires finite intensities");
                }
                x[i] = static_cast<double>(i) - static_cast<double>(length / 2);
            }
            PeakTemplate peak(model, length, fixedStep, options.minimumSamplesPerSide);
            PeakBackgroundValidation result;
            if (!peak.identifiable())
            {
                return result;
            }
            result.supportedPeakShape = peak.evaluate(values, x, shape);
            std::mt19937_64 engine(options.seed);
            result.constant = compareBackground(values, x, 0, shape, result.supportedPeakShape, peak, options.rotations, engine);
            result.linear = compareBackground(values, x, dot(x, x), shape, result.supportedPeakShape, peak, options.rotations, engine);
            if (std::isfinite(result.constant.pValue) && std::isfinite(result.linear.pValue))
            {
                result.status = std::max(result.constant.pValue, result.linear.pValue) <= options.alpha
                    ? PeakValidationStatus::passed : PeakValidationStatus::failed;
            }
            return result;
        }
    }
}
