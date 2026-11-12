#include "qalgorithms_peak_validation.h"
#include "cephes.h"

#include <cmath>

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

            template<typename CheckConfidence>
            PeakCandidateValidation checkCandidate(PeakGeometryValidation geometry, CheckConfidence checkConfidence)
            {
                PeakCandidateValidation result{geometry};
                if (geometry.apexSupport == PeakValidationStatus::passed)
                {
                    result.formConfidence = checkConfidence();
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
            const size_t max = std::numeric_limits<size_t>::max();
            if (fit.scale - 1 > max / 2)
            {
                return PeakValidationStatus::unavailable;
            }
            const size_t localDf = 2 * (fit.scale - 1);
            if (fit.locals.size() > max / localDf || df != fit.locals.size() * localDf - 1)
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

        PeakCandidateValidation validatePeakCandidate(const RegCoeffs &c, const PeakFormNoise &noise, size_t minimum)
        {
            return checkCandidate(validatePeakGeometry(c, minimum), [&]()
            {
                return validatePeakFormConfidence(c, noise.logMse, noise.residualDf, noise.alpha);
            });
        }

        PeakCandidateValidation validatePeakCandidate(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, size_t minimum)
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
            });
        }
    }
}
