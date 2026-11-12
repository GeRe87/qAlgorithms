#include "qalgorithms_peak_validation.h"
#include "cephes.h"
#include "libcerf_reduced.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace qAlgorithms
{
    namespace peak_model
    {
        namespace
        {
            PeakGeometryValidation checkGeometry(bool shape, double apex, size_t scale,
                double step, size_t minimumSamplesPerSide)
            {
                PeakGeometryValidation result;
                result.form = shape ? PeakValidationStatus::passed : PeakValidationStatus::failed;
                if (!shape)
                {
                    return result;
                }
                result.apex = apex;
                result.apexSupport = PeakValidationStatus::unavailable;
                if (scale == 0 || !(step > 0 && std::isfinite(step)))
                {
                    return result;
                }
                const double edge = static_cast<double>(scale) * step;
                if (!std::isfinite(edge))
                {
                    return result;
                }
                const double side = minimumSamplesPerSide == 0 ? 0
                    : (static_cast<double>(scale) - static_cast<double>(minimumSamplesPerSide - 1)) * step;
                const bool supported = std::isfinite(apex) && (minimumSamplesPerSide == 0
                    ? std::abs(apex) <= edge : std::abs(apex) < side);
                result.apexSupport = supported ? PeakValidationStatus::passed : PeakValidationStatus::failed;
                return result;
            }

            bool validNoise(double mse, size_t df, double alpha)
            {
                return mse > 0 && std::isfinite(mse) && df > 0 && alpha > 0 && alpha < 1;
            }

            bool matchesHyperbolicDf(const HyperbolicModel::Fit &fit, size_t df)
            {
                const size_t max = std::numeric_limits<size_t>::max();
                if (fit.locals.empty() || fit.scale < 2 || fit.scale - 1 > max / 2)
                {
                    return false;
                }
                const size_t localDf = 2 * (fit.scale - 1);
                return fit.locals.size() <= max / localDf && df == fit.locals.size() * localDf - 1;
            }

            double positiveContrastP(double estimate, double variance, size_t df)
            {
                if (!std::isfinite(estimate) || !(variance > 0 && std::isfinite(variance)))
                {
                    return std::numeric_limits<double>::quiet_NaN();
                }
                const double deviation = std::sqrt(variance);
                const double magnitude = std::abs(estimate);
                if (df == 1)
                {
                    const double tail = std::atan2(deviation, magnitude) / std::acos(-1.0);
                    return estimate >= 0 ? tail : 1 - tail;
                }
                const double degrees = static_cast<double>(df);
                const double scaledDeviation = std::sqrt(degrees) * deviation;
                const double ratio = scaledDeviation / std::hypot(scaledDeviation, magnitude);
                const double x = ratio * ratio;
                const double tail = x == 0 ? 0 : x == 1 ? .5 : .5 * cephes::incbeta(degrees / 2, .5, x);
                if (!(tail >= 0 && tail <= .5))
                {
                    return std::numeric_limits<double>::quiet_NaN();
                }
                return estimate >= 0 ? tail : 1 - tail;
            }

            constexpr double unavailable = std::numeric_limits<double>::quiet_NaN();

            double studentUpperQuantile(size_t df, double alpha)
            {
                if (!(alpha > 0 && alpha < .5))
                {
                    return unavailable;
                }
                // Invert the same one-sided Student tail used by form
                // confidence. Keep the upper bracket for a conservative bound.
                double left = 0, right = 1;
                for (;;)
                {
                    const double p = positiveContrastP(right, 1, df);
                    if (!std::isfinite(p))
                    {
                        return unavailable;
                    }
                    if (p <= alpha)
                    {
                        break;
                    }
                    left = right;
                    right *= 2;
                    if (!std::isfinite(right))
                    {
                        return unavailable;
                    }
                }
                for (int i = 0; i < 48; ++i)
                {
                    const double mid = left + (right - left) / 2;
                    const double p = positiveContrastP(mid, 1, df);
                    if (!std::isfinite(p))
                    {
                        return unavailable;
                    }
                    if (p > alpha)
                    {
                        left = mid;
                    }
                    else
                    {
                        right = mid;
                    }
                }
                return right;
            }

            struct CoefficientBounds
            {
                double lower, upper;
            };

            template<typename Variance>
            std::optional<std::array<CoefficientBounds, 4>> confidenceBox(
                const std::array<double, 4> &coefficients, const PeakFormNoise &noise, Variance variance)
            {
                if (!validNoise(noise.logMse, noise.residualDf, noise.alpha))
                {
                    return std::nullopt;
                }
                // Four two-sided intervals. Bonferroni needs no assumption
                // about independence of the fitted coefficients.
                const double critical = studentUpperQuantile(noise.residualDf, noise.alpha / 8);
                std::array<CoefficientBounds, 4> bounds{};
                for (size_t i = 0; i < bounds.size(); ++i)
                {
                    const double var = noise.logMse * variance(i);
                    if (!(var > 0 && std::isfinite(var)))
                    {
                        return std::nullopt;
                    }
                    const double radius = critical * std::sqrt(var);
                    bounds[i] = {coefficients[i] - radius, coefficients[i] + radius};
                    if (!std::isfinite(bounds[i].lower) || !std::isfinite(bounds[i].upper))
                    {
                        return std::nullopt;
                    }
                }
                return bounds;
            }

            PeakAreaEvidence uncertifiedArea()
            {
                PeakAreaEvidence result;
                result.status = PeakValidationStatus::failed;
                result.logLowerBound = -std::numeric_limits<double>::infinity();
                return result;
            }

            template<typename MeanLog>
            PeakAreaEvidence testArea(double left, double right, double minimum, MeanLog meanLog)
            {
                PeakAreaEvidence result;
                if (!(left < 0 && right > 0 && std::isfinite(left) && std::isfinite(right)))
                {
                    return result;
                }
                // Jensen's inequality on each bin bounds the integral from
                // below. Splitting at zero keeps basis signs fixed in each bin.
                // Sum in log space to handle unrepresentable raw amplitudes.
                constexpr size_t bins = 32;
                std::array<double, 2 * bins> terms{};
                size_t i = 0;
                for (const auto &interval : {std::array<double, 2>{{left, 0}}, std::array<double, 2>{{0, right}}})
                {
                    const double width = (interval[1] - interval[0]) / bins;
                    for (size_t k = 0; k < bins; ++k)
                    {
                        const double midpoint = interval[0] + (static_cast<double>(k) + .5) * width;
                        terms[i++] = std::log(width) + meanLog(midpoint, width / 2);
                    }
                }
                if (!std::all_of(terms.begin(), terms.end(), [](double term) { return std::isfinite(term); }))
                {
                    return result;
                }
                const double maximum = *std::max_element(terms.begin(), terms.end());
                double sum = 0;
                for (double term : terms)
                {
                    sum += std::exp(term - maximum);
                }
                result.logLowerBound = maximum + std::log(sum);
                if (!std::isfinite(result.logLowerBound))
                {
                    return {};
                }
                result.left = left;
                result.right = right;
                result.shapeCertified = true;
                result.status = result.logLowerBound > std::log(minimum)
                    ? PeakValidationStatus::passed : PeakValidationStatus::failed;
                return result;
            }

            template<typename Curve>
            double integrateWindow(const Curve &curve, double left, double right, double apex, double total)
            {
                if (!(total > 0 && std::isfinite(total) && std::isfinite(apex)
                    && std::isfinite(left) && std::isfinite(right) && left < right))
                {
                    return unavailable;
                }
                // Split at the mode and the quadratic join so the integrator
                // sees the narrow peak even when the window is very wide.
                std::array<double, 4> points{{left, std::clamp(apex, left, right),
                    std::clamp(0.0, left, right), right}};
                std::sort(points.begin(), points.end());
                auto simpson = [](double width, double a, double m, double b)
                {
                    return width * (a + 4 * m + b) / 6;
                };
                size_t remaining = 65536;
                auto refine = [&](auto &self, double a, double b, double fa, double fm, double fb,
                    double previous, double tolerance, int depth) -> double
                {
                    if (remaining == 0)
                    {
                        return unavailable;
                    }
                    --remaining;
                    const double mid = a + (b - a) / 2;
                    const double fl = curve(a + (mid - a) / 2), fr = curve(mid + (b - mid) / 2);
                    const double l = simpson(mid - a, fa, fl, fm), r = simpson(b - mid, fm, fr, fb);
                    const double difference = l + r - previous;
                    if (!std::isfinite(difference))
                    {
                        return unavailable;
                    }
                    if (std::abs(difference) <= 15 * tolerance)
                    {
                        return l + r + difference / 15;
                    }
                    if (depth == 0 || mid == a || mid == b)
                    {
                        return unavailable;
                    }
                    const double first = self(self, a, mid, fa, fl, fm, l, tolerance / 2, depth - 1);
                    return std::isfinite(first) ? first + self(self, mid, b, fm, fr, fb, r, tolerance / 2, depth - 1)
                        : unavailable;
                };
                double window = 0;
                for (size_t i = 1; i < points.size(); ++i)
                {
                    const double a = points[i - 1], b = points[i];
                    if (a < b)
                    {
                        const double fa = curve(a), fm = curve(a + (b - a) / 2), fb = curve(b);
                        window += refine(refine, a, b, fa, fm, fb, simpson(b - a, fa, fm, fb), 1e-11 * total, 24);
                    }
                }
                const double fraction = window / total;
                return fraction >= 0 && fraction <= 1 + 1e-9 ? std::min(1.0, fraction) : unavailable;
            }

            double quadraticHalfArea(double slope, double curvature, double peakLog)
            {
                const double root = std::sqrt(std::abs(curvature));
                if (curvature == 0)
                {
                    return std::exp(-peakLog - std::log(-slope));
                }
                const double z = -slope / root / 2;
                if (curvature > 0)
                {
                    return std::exp(-peakLog + std::log(libcerf::dawson(z)) - std::log(root));
                }
                const double scaledErfc = z >= 0 ? libcerf::erfcx(z) : std::erfc(z);
                const double exponent = z >= 0 ? -peakLog : z * z - peakLog;
                return std::exp(exponent + std::log(scaledErfc) - std::log(root)
                    + std::log(std::sqrt(std::acos(-1.0)) / 2));
            }

            double quadraticCoverage(const RegCoeffs &c)
            {
                if (!hasPeakShape(c))
                {
                    return unavailable;
                }
                const double apex = Model<PeakModelKind::PiecewiseQuadratic>::apex(&c);
                const double peakLog = c.b1 * apex / 2;
                const double scale = c.scale;
                const double left = c.b2 > 0 ? std::max(-scale, -c.b1 / (2 * c.b2)) : -scale;
                const double right = c.b3 > 0 ? std::min(scale, -c.b1 / (2 * c.b3)) : scale;
                const double total = quadraticHalfArea(-c.b1, c.b2, peakLog)
                    + quadraticHalfArea(c.b1, c.b3, peakLog);
                return integrateWindow([&](double x)
                {
                    return std::exp((c.b1 + x * (x < 0 ? c.b2 : c.b3)) * x - peakLog);
                }, left, right, apex, total);
            }

            double hyperbolicCoverage(const HyperbolicModel::Coefficients &c, double slope, double edge)
            {
                if (!hasPeakShape(c, slope))
                {
                    return unavailable;
                }
                auto shape = c;
                shape.intercept = 0;
                const double apex = detail::hyperbolic::apexU(shape, slope);
                const double peakLog = detail::hyperbolic::predict(shape, slope, apex);
                double total;
                try
                {
                    total = std::exp(detail::hyperbolic::logAreaU(shape, slope) - peakLog);
                }
                catch (const std::domain_error &)
                {
                    return unavailable;
                }
                catch (const std::runtime_error &)
                {
                    return unavailable;
                }
                return integrateWindow([&](double u)
                {
                    return std::exp(detail::hyperbolic::predict(shape, slope, u) - peakLog);
                }, -edge, edge, apex, total);
            }

            template<typename Fraction, typename Covariance>
            PeakWindowCoverage testCoverage(const std::array<double, 3> &coefficients, const PeakFormNoise &noise,
                double minimum, Fraction fraction, Covariance covariance)
            {
                PeakWindowCoverage result;
                result.fraction = fraction(coefficients);
                if (!(result.fraction > 0 && result.fraction < 1) || !validNoise(noise.logMse, noise.residualDf, noise.alpha)
                    || !(noise.alpha < .5 && minimum > 0 && minimum < 1))
                {
                    return result;
                }
                std::array<double, 3> gradient{};
                for (size_t i = 0; i < gradient.size(); ++i)
                {
                    const double step = 1e-4 * std::max(std::abs(coefficients[i]), 1e-3);
                    auto differences = [&](double delta)
                    {
                        auto left = coefficients, right = coefficients;
                        left[i] -= delta;
                        right[i] += delta;
                        return std::array<double, 2>{{
                            (result.fraction - fraction(left)) / (coefficients[i] - left[i]),
                            (fraction(right) - result.fraction) / (right[i] - coefficients[i])}};
                    };
                    const auto coarse = differences(step), fine = differences(step / 2);
                    const double left = 2 * fine[0] - coarse[0], right = 2 * fine[1] - coarse[1];
                    const double derivative = (left + right) / 2;
                    const double tolerance = 1e-3 * std::max(std::abs(derivative), 1e-5);
                    // A clipping boundary at a window edge can have two
                    // different derivatives despite stable central differences.
                    if (!std::isfinite(derivative) || !(std::abs(left - right) <= tolerance))
                    {
                        return result;
                    }
                    gradient[i] = derivative;
                }
                const double variance = noise.logMse * covariance(gradient);
                if (!(variance > 0 && std::isfinite(variance)))
                {
                    return result;
                }
                const double critical = studentUpperQuantile(noise.residualDf, noise.alpha);
                if (!std::isfinite(critical))
                {
                    return result;
                }
                result.standardError = std::sqrt(variance);
                result.lowerConfidenceBound = std::max(0.0, result.fraction - critical * result.standardError);
                result.status = result.lowerConfidenceBound >= minimum ? PeakValidationStatus::passed : PeakValidationStatus::failed;
                return result;
            }

            template<typename CheckConfidence, typename CheckCoverage, typename CheckArea>
            PeakCandidateValidation checkCandidate(PeakGeometryValidation geometry, CheckConfidence checkConfidence,
                std::optional<double> minimumCoverage, CheckCoverage checkCoverage,
                std::optional<double> minimumArea, CheckArea checkArea)
            {
                PeakCandidateValidation result{geometry};
                if (geometry.apexSupport == PeakValidationStatus::passed)
                {
                    result.formConfidence = checkConfidence();
                }
                if (result.formConfidence != PeakValidationStatus::passed)
                {
                    return result;
                }
                if (minimumArea)
                {
                    result.areaEvidence = checkArea(*minimumArea);
                    if (result.areaEvidence.status != PeakValidationStatus::passed)
                    {
                        return result;
                    }
                }
                if (minimumCoverage)
                {
                    result.windowCoverage = checkCoverage(*minimumCoverage);
                }
                return result;
            }
        }

        bool hasPeakShape(const RegCoeffs &c)
        {
            return std::isfinite(c.b1) && std::isfinite(c.b2) && std::isfinite(c.b3)
                && Model<PeakModelKind::PiecewiseQuadratic>::hasLegacyPeakShape(c) && (c.b1 != 0 || c.b2 < 0);
        }

        bool hasPeakShape(const HyperbolicModel::Coefficients &c, double slope)
        {
            return std::isfinite(slope) && detail::hyperbolic::validIntegrablePeak(c);
        }

        PeakGeometryValidation validatePeakGeometry(const RegCoeffs &c, size_t minimum)
        {
            const bool shape = hasPeakShape(c);
            return checkGeometry(shape, shape ? Model<PeakModelKind::PiecewiseQuadratic>::apex(&c) : 0,
                c.scale, 1, minimum);
        }

        PeakGeometryValidation validatePeakGeometry(const HyperbolicModel::Coefficients &c, double slope,
            size_t scale, double step, size_t minimum)
        {
            const bool shape = hasPeakShape(c, slope);
            return checkGeometry(shape, shape ? detail::hyperbolic::apexU(c, slope) : 0, scale, step, minimum);
        }

        PeakValidationStatus validatePeakFormConfidence(const RegCoeffs &c, double mse, size_t df, double alpha)
        {
            if (!hasPeakShape(c))
            {
                return PeakValidationStatus::failed;
            }
            if (!validNoise(mse, df, alpha) || c.scale < 2 || c.scale > maxscale_global
                || df != 2 * static_cast<size_t>(c.scale) - 3)
            {
                return PeakValidationStatus::unavailable;
            }
            using Quadratic = Model<PeakModelKind::PiecewiseQuadratic>;
            const double slope[4] = {0, 1, 0, 0}, left[4] = {0, 0, 1, 0}, right[4] = {0, 0, 0, 1};
            const double slopeVariance = mse * Quadratic::covarianceProjection(slope, c.scale);
            const double leftP = positiveContrastP(-c.b2, mse * Quadratic::covarianceProjection(left, c.scale), df);
            const double rightP = positiveContrastP(-c.b3, mse * Quadratic::covarianceProjection(right, c.scale), df);
            const double negativeP = positiveContrastP(-c.b1, slopeVariance, df);
            const double positiveP = positiveContrastP(c.b1, slopeVariance, df);
            if (!std::isfinite(leftP) || !std::isfinite(rightP) || !std::isfinite(negativeP) || !std::isfinite(positiveP))
            {
                return PeakValidationStatus::unavailable;
            }
            // Simultaneous one-sided statements certify one branch of the
            // shape region without treating the estimated slope sign as fixed.
            const double limit = alpha / 4;
            const bool passes = (leftP <= limit && (rightP <= limit || negativeP <= limit))
                || (rightP <= limit && positiveP <= limit);
            return passes ? PeakValidationStatus::passed : PeakValidationStatus::failed;
        }

        PeakValidationStatus validatePeakFormConfidence(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, double mse, size_t df, double alpha)
        {
            if (trace >= fit.locals.size() || fit.scale < 2 || fit.scale > geometry.maxScale()
                || !fit.identifiable || !validNoise(mse, df, alpha))
            {
                return PeakValidationStatus::unavailable;
            }
            if (!matchesHyperbolicDf(fit, df))
            {
                return PeakValidationStatus::unavailable;
            }
            const auto &c = fit.locals[trace];
            if (!hasPeakShape(c, fit.sharedSlope))
            {
                return PeakValidationStatus::failed;
            }
            const auto &window = geometry.at(fit.scale);
            if (!window.identifiable)
            {
                return PeakValidationStatus::unavailable;
            }
            // C and D come from orthogonal even and odd regressors.
            const double variance = mse * (window.varCurvatureFactor() + window.varLocalSkewFactor(fit.locals.size()));
            const double leftP = positiveContrastP(c.curvature - c.skew, variance, df);
            const double rightP = positiveContrastP(c.curvature + c.skew, variance, df);
            if (!std::isfinite(leftP) || !std::isfinite(rightP))
            {
                return PeakValidationStatus::unavailable;
            }
            // The null is that at least one flank is invalid.
            return leftP <= alpha && rightP <= alpha ? PeakValidationStatus::passed : PeakValidationStatus::failed;
        }

        PeakCandidateValidation validatePeakCandidate(const RegCoeffs &c, const PeakFormNoise &noise, size_t minimum,
            std::optional<double> minimumCoverage, std::optional<double> minimumArea)
        {
            return checkCandidate(validatePeakGeometry(c, minimum), [&]()
            {
                return validatePeakFormConfidence(c, noise.logMse, noise.residualDf, noise.alpha);
            }, minimumCoverage, [&](double threshold) { return validatePeakWindowCoverage(c, noise, threshold); },
                minimumArea, [&](double threshold) { return validatePeakAreaEvidence(c, noise, threshold); });
        }

        PeakCandidateValidation validatePeakCandidate(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, size_t minimum,
            std::optional<double> minimumCoverage, std::optional<double> minimumArea)
        {
            if (trace >= fit.locals.size() || fit.scale < 2 || fit.scale > geometry.maxScale())
            {
                PeakCandidateValidation result;
                result.geometry.form = PeakValidationStatus::unavailable;
                return result;
            }
            return checkCandidate(validatePeakGeometry(fit.locals[trace], fit.sharedSlope, fit.scale, geometry.step(), minimum), [&]()
            {
                return validatePeakFormConfidence(fit, trace, geometry, noise.logMse, noise.residualDf, noise.alpha);
            }, minimumCoverage, [&](double threshold) { return validatePeakWindowCoverage(fit, trace, geometry, noise, threshold); },
                minimumArea, [&](double threshold) { return validatePeakAreaEvidence(fit, trace, geometry, noise, threshold); });
        }

        PeakWindowCoverage validatePeakWindowCoverage(const RegCoeffs &c, const PeakFormNoise &noise, double minimum)
        {
            if (c.scale < 2 || c.scale > maxscale_global)
            {
                return {};
            }
            auto effectiveNoise = noise;
            if (noise.residualDf != 2 * static_cast<size_t>(c.scale) - 3)
            {
                effectiveNoise.residualDf = 0;
            }
            return testCoverage({{c.b1, c.b2, c.b3}}, effectiveNoise, minimum, [&](const auto &p)
            {
                auto coefficients = c;
                coefficients.b1 = p[0];
                coefficients.b2 = p[1];
                coefficients.b3 = p[2];
                return quadraticCoverage(coefficients);
            }, [&](const auto &g)
            {
                const double jacobian[4] = {0, g[0], g[1], g[2]};
                return Model<PeakModelKind::PiecewiseQuadratic>::covarianceProjection(jacobian, c.scale);
            });
        }

        PeakWindowCoverage validatePeakWindowCoverage(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, double minimum)
        {
            if (trace >= fit.locals.size() || fit.scale < 2 || fit.scale > geometry.maxScale() || !fit.identifiable)
            {
                return {};
            }
            const auto &window = geometry.at(fit.scale);
            if (!window.identifiable)
            {
                return {};
            }
            const auto &c = fit.locals[trace];
            auto effectiveNoise = noise;
            if (!matchesHyperbolicDf(fit, noise.residualDf))
            {
                effectiveNoise.residualDf = 0;
            }
            return testCoverage({{fit.sharedSlope, c.curvature, c.skew}}, effectiveNoise, minimum, [&](const auto &p)
            {
                return hyperbolicCoverage({0, p[1], p[2]}, p[0], window.edgeU);
            }, [&](const auto &g)
            {
                const size_t traces = fit.locals.size();
                const double slope = g[0], curvature = g[1], skew = g[2];
                return slope * slope * window.varSharedSlopeFactor(traces)
                    + curvature * curvature * window.varCurvatureFactor()
                    + skew * skew * window.varLocalSkewFactor(traces)
                    + 2 * slope * skew * window.covSharedSlopeLocalSkewFactor(traces);
            });
        }

        PeakAreaEvidence validatePeakAreaEvidence(const RegCoeffs &c, const PeakFormNoise &noise, double minimumArea)
        {
            if (!(minimumArea > 0 && std::isfinite(minimumArea)) || !hasPeakShape(c)
                || c.scale < 2 || c.scale > maxscale_global || noise.residualDf != 2 * static_cast<size_t>(c.scale) - 3)
            {
                return {};
            }
            const auto box = confidenceBox({{c.b0, c.b1, c.b2, c.b3}}, noise, [&](size_t i)
            {
                double jacobian[4] = {};
                jacobian[i] = 1;
                return Model<PeakModelKind::PiecewiseQuadratic>::covarianceProjection(jacobian, c.scale);
            });
            if (!box)
            {
                return {};
            }
            const auto &[intercept, slope, leftCurve, rightCurve] = *box;
            double left = -static_cast<double>(c.scale), right = static_cast<double>(c.scale);
            if (leftCurve.upper >= 0)
            {
                if (!(slope.lower > 0 && rightCurve.upper < 0))
                {
                    return uncertifiedArea();
                }
                if (leftCurve.upper > 0)
                {
                    left = std::max(left, -slope.lower / leftCurve.upper / 2);
                }
            }
            if (rightCurve.upper >= 0)
            {
                if (!(slope.upper < 0 && leftCurve.upper < 0))
                {
                    return uncertifiedArea();
                }
                if (rightCurve.upper > 0)
                {
                    right = std::min(right, -slope.upper / rightCurve.upper / 2);
                }
            }
            return testArea(left, right, minimumArea, [&](double mid, double halfWidth)
            {
                const double linear = mid < 0 ? slope.upper : slope.lower;
                const double quadratic = mid < 0 ? leftCurve.lower : rightCurve.lower;
                const double meanSquared = mid * mid + halfWidth * halfWidth / 3;
                return intercept.lower + linear * mid + quadratic * meanSquared;
            });
        }

        PeakAreaEvidence validatePeakAreaEvidence(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, double minimumArea)
        {
            if (!(minimumArea > 0 && std::isfinite(minimumArea)) || !fit.identifiable || trace >= fit.locals.size()
                || fit.scale < 2 || fit.scale > geometry.maxScale() || !matchesHyperbolicDf(fit, noise.residualDf))
            {
                return {};
            }
            const auto &window = geometry.at(fit.scale);
            const auto &c = fit.locals[trace];
            if (!window.identifiable || !hasPeakShape(c, fit.sharedSlope))
            {
                return {};
            }
            const size_t traces = fit.locals.size();
            const std::array<double, 4> variances{{window.varInterceptFactor(), window.varSharedSlopeFactor(traces),
                window.varCurvatureFactor(), window.varLocalSkewFactor(traces)}};
            const auto box = confidenceBox({{c.intercept, fit.sharedSlope, c.curvature, c.skew}}, noise,
                [&](size_t i) { return variances[i]; });
            if (!box)
            {
                return {};
            }
            const auto &[intercept, slope, curvature, skew] = *box;
            if (!(curvature.lower > std::max(std::abs(skew.lower), std::abs(skew.upper))))
            {
                return uncertifiedArea();
            }
            return testArea(-window.edgeU, window.edgeU, minimumArea, [&](double mid, double halfWidth)
            {
                const double linear = mid < 0 ? slope.upper : slope.lower;
                const double skewTerm = mid < 0 ? skew.lower : skew.upper;
                const double factor = std::sinh(halfWidth) / halfWidth;
                return intercept.lower + linear * mid
                    - factor * (curvature.upper * std::cosh(mid) + skewTerm * std::sinh(mid));
            });
        }
    }
}
