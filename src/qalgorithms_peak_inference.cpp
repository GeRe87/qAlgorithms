#include "qalgorithms_peak_inference.h"
#include "detail/qalgorithms_peak_inference.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace qAlgorithms
{
    namespace peak_model
    {
        namespace
        {
            namespace hyperbolic = detail::hyperbolic;

            // For numerator DF=2m the F survival is I_x(a,m), a=nu/2.
            // Integer m gives a finite positive series. Log evaluation avoids
            // cancellation and preserves decisions beyond p-value underflow.
            double logComparisonSurvival(double statistic, double residualDf, size_t traces)
            {
                if (!(statistic > 0))
                {
                    return 0;
                }
                const double logRatio = std::log(2.0 * traces) + std::log(statistic) - std::log(residualDf);
                const double logX = logRatio > 0
                    ? -logRatio - std::log1p(std::exp(-logRatio))
                    : -std::log1p(std::exp(logRatio));
                const double logComplement = std::log(-std::expm1(logX));
                const double shape = residualDf / 2;
                double logTerm = 0, logSum = 0;
                for (size_t k = 1; k < traces; ++k)
                {
                    logTerm += std::log(shape + k - 1) - std::log(static_cast<double>(k)) + logComplement;
                    const double larger = std::max(logSum, logTerm);
                    logSum = larger + std::log1p(std::exp(-std::abs(logSum - logTerm)));
                }
                return std::min(0.0, shape * logX + logSum);
            }

            struct PropertyInfluence
            {
                std::vector<long double> local;
                long double shared = 0;
                long double localProjection = 0;
                double derivativeScale = 1;
            };

            double checkedLogIntegral(const hyperbolic::LocalFit &peak, double slope)
            {
                const double value = hyperbolic::logAreaU(peak, slope);
                if (!std::isfinite(value))
                {
                    throw std::runtime_error("Hyperbolic integral outside supported numeric range");
                }
                return value;
            }

            // Derivative of log integral with respect to (b,C,D). C,D use
            // Bessel recurrences; the order derivative uses checked Richardson
            // extrapolation. Failure to converge is reported, never set to zero.
            std::array<double, 3> integralGradient(const hyperbolic::LocalFit &peak, double slope)
            {
                const double z = std::sqrt((peak.curvature - peak.skew) * (peak.curvature + peak.skew));
                const double base = std::cyl_bessel_k(std::abs(slope), z);
                const double lower = std::cyl_bessel_k(std::abs(slope - 1), z);
                const double upper = std::cyl_bessel_k(std::abs(slope + 1), z);
                const double argumentDerivative = -0.5 * (lower + upper) / base;
                hyperbolic::LocalFit unit = peak;
                unit.intercept = 0;
                const double step = std::cbrt(std::numeric_limits<double>::epsilon()) * (1 + std::abs(slope));
                const double coarse = (checkedLogIntegral(unit, slope + step) - checkedLogIntegral(unit, slope - step)) / (2 * step);
                const double fine = (checkedLogIntegral(unit, slope + step / 2) - checkedLogIntegral(unit, slope - step / 2)) / step;
                if (!std::isfinite(fine) || std::abs(fine - coarse) > 1e-5 * (1 + std::abs(fine)))
                {
                    throw std::runtime_error("Hyperbolic order derivative did not converge");
                }
                const double denominator = z * z;
                std::array<double, 3> gradient = {{
                    (4 * fine - coarse) / 3,
                    slope * peak.skew / denominator + argumentDerivative * peak.curvature / z,
                    -slope * peak.curvature / denominator - argumentDerivative * peak.skew / z}};
                for (double value : gradient)
                {
                    if (!std::isfinite(value))
                    {
                        throw std::runtime_error("Nonfinite hyperbolic integral gradient");
                    }
                }
                return gradient;
            }

            // Full derivative through BOTH raw observations and log-OLS shape.
            // The direct term w_i belongs only to this trace. Shared b couples
            // all traces through r_i = u_i-kappa*sinh(u_i).
            PropertyInfluence makeInfluence(
                const std::vector<double> &u,
                const std::vector<long double> &weights,
                std::array<double, 3> gradient,
                bool includeAmplitude,
                size_t traces,
                const hyperbolic::WindowGeometry &geometry)
            {
                const size_t count = u.size();
                PropertyInfluence influence;
                influence.local.resize(count);
                if (includeAmplitude)
                {
                    long double weightU = 0, weightCosh = 0, weightSinh = 0;
                    long double squaredU = 0, squaredCosh = 0, squaredSinh = 0;
                    // weights hold normalized raw*f; shape-squared weights are
                    // constructed separately by the caller below.
                    for (size_t i = 0; i < count; ++i)
                    {
                        const long double direct = weights[i];
                        const long double squared = weights[count + i];
                        weightU += direct * u[i];
                        weightCosh += direct * std::cosh(u[i]);
                        weightSinh += direct * std::sinh(u[i]);
                        squaredU += squared * u[i];
                        squaredCosh += squared * std::cosh(u[i]);
                        squaredSinh += squared * std::sinh(u[i]);
                    }
                    gradient[0] += static_cast<double>(weightU - 2 * squaredU);
                    gradient[1] += static_cast<double>(-weightCosh + 2 * squaredCosh);
                    gradient[2] += static_cast<double>(-weightSinh + 2 * squaredSinh);
                }
                influence.shared = (gradient[0] + geometry.kappa * gradient[2]) * geometry.invOddResidualNorm2 / traces;
                for (size_t i = 0; i < count; ++i)
                {
                    const double cosine = std::cosh(u[i]);
                    const double sine = std::sinh(u[i]);
                    const double residual = u[i] - geometry.kappa * sine;
                    influence.local[i] = (includeAmplitude ? weights[i] : 0)
                        - gradient[1] * (cosine - geometry.meanCosh) * geometry.invCenteredCoshNorm2
                        - gradient[2] * sine * geometry.invSinhNorm2;
                    influence.localProjection += influence.local[i] * residual;
                }
                return influence;
            }
        }

        static HyperbolicInference inferWindow(
            const std::vector<const float *> &rawTraces,
            size_t length,
            size_t center,
            size_t scale,
            double fixedStep,
            const HyperbolicInferenceOptions &options, const HyperbolicModel::Geometry *prepared)
        {
            if (rawTraces.empty() || center >= length || scale < 2 || scale > center || scale > length - 1 - center
                || options.familySize == 0 || !(options.familyAlpha > 0 && options.familyAlpha < 1)
                || !std::isfinite(options.familyAlpha)
                || (options.selection != WindowSelection::exploratory && options.selection != WindowSelection::prespecified))
            {
                throw std::invalid_argument("Invalid selected-window inference request");
            }
            const size_t traces = rawTraces.size();
            if (scale > (std::numeric_limits<size_t>::max() - 1) / 2)
            {
                throw std::invalid_argument("Inference window size overflow");
            }
            const size_t count = 2 * scale + 1;
            if (traces > (std::numeric_limits<size_t>::max() - 1) / count
                || traces > std::numeric_limits<size_t>::max() / 3)
            {
                throw std::invalid_argument("Inference trace count overflow");
            }
            std::optional<hyperbolic::GeometryCache> localCache;
            if (!prepared)
            {
                localCache.emplace(scale, fixedStep);
            }
            const auto &cache = prepared ? *prepared : *localCache;
            const auto &geometry = cache.at(scale);
            if (!geometry.identifiable)
            {
                HyperbolicInference output;
                output.residualDegreesOfFreedom = traces * (count - 3) - 1;
                output.comparisonDegreesOfFreedom = 2 * traces;
                output.status = InferenceStatus::invalidFit;
                return output;
            }
            std::vector<std::vector<float>> logs(traces, std::vector<float>(count));
            std::vector<const float *> logPointers;
            std::vector<const float *> selectedRaw;
            for (size_t j = 0; j < traces; ++j)
            {
                if (rawTraces[j] == nullptr)
                {
                    throw std::invalid_argument("Null raw trace");
                }
                const float *raw = rawTraces[j] + center - scale;
                selectedRaw.push_back(raw);
                for (size_t i = 0; i < count; ++i)
                {
                    if (!(raw[i] > 0 && std::isfinite(raw[i])))
                    {
                        throw std::invalid_argument("Inference requires positive finite uncensored raw observations");
                    }
                    logs[j][i] = static_cast<float>(std::log(static_cast<double>(raw[i])));
                }
                logPointers.push_back(logs[j].data());
            }
            hyperbolic::FitView fit;
            hyperbolic::forEachRunning(logPointers, count, cache, [&](const hyperbolic::FitView &candidate)
            {
                if (candidate.center == scale && candidate.scale == scale)
                {
                    fit = candidate;
                }
            });
            return inferFittedHyperbolicWindow(selectedRaw, logPointers, count, fit, cache, options, false);
        }

        HyperbolicInference inferFittedHyperbolicWindow(
            const std::vector<const float *> &rawTraces, const std::vector<const float *> &logPointers,
            size_t length, const HyperbolicModel::Fit &fit, const HyperbolicModel::Geometry &cache,
            const HyperbolicInferenceOptions &options, bool supportedOnly)
        {
            const size_t traces = rawTraces.size(), scale = fit.scale, count = 2 * scale + 1;
            const size_t offset = fit.center - scale;
            const auto &geometry = cache.at(scale);
            const double fixedStep = cache.step();
            HyperbolicInference output;
            output.residualDegreesOfFreedom = traces * (count - 3) - 1;
            output.comparisonDegreesOfFreedom = 2 * traces;
            if (!geometry.identifiable)
            {
                output.status = InferenceStatus::invalidFit;
                return output;
            }
            std::vector<double> u(count);
            for (size_t i = 0; i < count; ++i)
            {
                u[i] = fixedStep * (static_cast<double>(i) - static_cast<double>(scale));
            }
            if (supportedOnly)
            {
                for (const auto &peak : fit.locals)
                {
                    const double apex = hyperbolic::apexU(peak, fit.sharedSlope);
                    const auto left = std::lower_bound(u.begin(), u.end(), apex) - u.begin();
                    const auto right = u.end() - std::upper_bound(u.begin(), u.end(), apex);
                    if (!std::isfinite(apex) || std::abs(apex) > geometry.edgeU
                        || static_cast<size_t>(left) < options.minimumSamplesPerSide
                        || static_cast<size_t>(right) < options.minimumSamplesPerSide)
                    {
                        output.status = InferenceStatus::screenedOut;
                        return output;
                    }
                }
            }
            std::vector<hyperbolic::Moments> moments(traces);
            for (size_t j = 0; j < traces; ++j)
            {
                for (size_t i = 0; i < count; ++i)
                {
                    moments[j].sumY += logPointers[j][offset + i];
                    moments[j].sumUY += u[i] * logPointers[j][offset + i];
                }
            }
            long double sumU2 = 0, nullSlopeNumerator = 0;
            for (double coordinate : u)
            {
                sumU2 += coordinate * coordinate;
            }
            for (size_t j = 0; j < traces; ++j)
            {
                nullSlopeNumerator += moments[j].sumUY;
            }
            const long double nullSlope = nullSlopeNumerator / (traces * sumU2);
            long double residual = 0, nullResidual = 0;
            for (size_t j = 0; j < traces; ++j)
            {
                for (size_t i = 0; i < count; ++i)
                {
                    const long double error = logPointers[j][offset + i] - hyperbolic::predict(fit.locals[j], fit.sharedSlope, u[i]);
                    const long double nullError = logPointers[j][offset + i] - moments[j].sumY / count - nullSlope * u[i];
                    residual += error * error;
                    nullResidual += nullError * nullError;
                }
            }
            output.residualSumSquares = static_cast<double>(residual);
            output.nullResidualSumSquares = static_cast<double>(nullResidual);
            if (!(residual > 0 && std::isfinite(residual) && std::isfinite(nullResidual)))
            {
                output.status = InferenceStatus::degenerateNoise;
                return output;
            }
            output.logNoiseVariance = static_cast<double>(residual / output.residualDegreesOfFreedom);
            if (!(output.logNoiseVariance > 0 && std::isfinite(output.logNoiseVariance)))
            {
                output.status = InferenceStatus::degenerateNoise;
                return output;
            }
            const long double improvement = nullResidual - residual;
            const long double roundingTolerance = 128 * std::numeric_limits<double>::epsilon() * (1 + nullResidual);
            if (improvement < -roundingTolerance)
            {
                output.status = InferenceStatus::numericalFailure;
                return output;
            }
            output.comparisonF = static_cast<double>(std::max(0.0L, improvement) / output.comparisonDegreesOfFreedom / output.logNoiseVariance);
            output.logNominalP = logComparisonSurvival(output.comparisonF,
                static_cast<double>(output.residualDegreesOfFreedom), traces);
            output.logAdjustedP = std::min(0.0, output.logNominalP + std::log(static_cast<double>(options.familySize)));
            output.nominalP = std::exp(output.logNominalP);
            output.adjustedP = std::exp(output.logAdjustedP);
            if (std::isnan(output.logNominalP) || !(output.nominalP >= 0 && output.nominalP <= 1))
            {
                output.status = InferenceStatus::numericalFailure;
                return output;
            }
            bool supported = true;
            std::vector<PropertyInfluence> influences;
            influences.reserve(3 * traces);
            try
            {
                for (size_t j = 0; j < traces; ++j)
                {
                    auto property = evaluateHyperbolicWindow(fit, j, cache, rawTraces[j], length);
                    output.properties.push_back(property);
                    if (!property.finiteCoefficients || !property.integrable || !property.amplitude.valid
                        || property.apexU.status != PropertyStatus::Available
                        || property.height.status != PropertyStatus::Available || property.areaU.status != PropertyStatus::Available)
                    {
                        output.status = InferenceStatus::invalidFit;
                        return output;
                    }
                    size_t left = 0, right = 0;
                    for (double coordinate : u)
                    {
                        left += coordinate < property.apexU.value;
                        right += coordinate > property.apexU.value;
                    }
                    supported = supported && property.apexInObservedSupport
                        && left >= options.minimumSamplesPerSide && right >= options.minimumSamplesPerSide;
                    const auto &peak = fit.locals[j];
                    std::vector<long double> weights(2 * count);
                    double maxShape = -std::numeric_limits<double>::infinity();
                    for (double coordinate : u)
                    {
                        maxShape = std::max(maxShape, hyperbolic::predict(peak, fit.sharedSlope, coordinate) - peak.intercept);
                    }
                    long double numerator = 0, denominator = 0;
                    for (size_t i = 0; i < count; ++i)
                    {
                        const long double shape = std::exp(hyperbolic::predict(peak, fit.sharedSlope, u[i]) - peak.intercept - maxShape);
                        weights[i] = rawTraces[j][offset + i] * shape;
                        weights[count + i] = shape * shape;
                        numerator += weights[i];
                        denominator += weights[count + i];
                    }
                    for (size_t i = 0; i < count; ++i)
                    {
                        weights[i] /= numerator;
                        weights[count + i] /= denominator;
                    }
                    const double apex = property.apexU.value;
                    const double apexCurvature = peak.curvature * std::cosh(apex) + peak.skew * std::sinh(apex);
                    const std::array<double, 3> apexGradient = {{
                        1 / apexCurvature,
                        -std::sinh(apex) / apexCurvature,
                        -std::cosh(apex) / apexCurvature}};
                    const std::array<double, 3> heightGradient = {{
                        apex, -std::cosh(apex), -std::sinh(apex)}};
                    influences.push_back(makeInfluence(u, weights, apexGradient, false, traces, geometry));
                    influences.push_back(makeInfluence(u, weights, heightGradient, true, traces, geometry));
                    influences.back().derivativeScale = property.height.value;
                    influences.push_back(makeInfluence(u, weights, integralGradient(peak, fit.sharedSlope), true, traces, geometry));
                    influences.back().derivativeScale = property.areaU.value;
                }
                const size_t dimension = influences.size();
                if (dimension > std::numeric_limits<size_t>::max() / dimension)
                {
                    throw std::runtime_error("Property covariance size overflow");
                }
                output.covariance.resize(dimension * dimension);
                const long double oddNorm = 1.0L / geometry.invOddResidualNorm2;
                for (size_t row = 0; row < dimension; ++row)
                {
                    for (size_t column = row; column < dimension; ++column)
                    {
                        const auto &a = influences[row];
                        const auto &b = influences[column];
                        long double factor = traces * a.shared * b.shared * oddNorm
                            + a.shared * b.localProjection + b.shared * a.localProjection;
                        if (row / 3 == column / 3)
                        {
                            for (size_t i = 0; i < count; ++i)
                            {
                                factor += a.local[i] * b.local[i];
                            }
                        }
                        const double value = static_cast<double>(factor * output.logNoiseVariance * a.derivativeScale * b.derivativeScale);
                        if (!std::isfinite(value) || (row == column && value < 0))
                        {
                            throw std::runtime_error("Invalid property covariance");
                        }
                        output.covariance[row * dimension + column] = value;
                        output.covariance[column * dimension + row] = value;
                    }
                }
                for (size_t j = 0; j < traces; ++j)
                {
                    const size_t areaIndex = 3 * j + 2;
                    output.properties[j].areaUncertainty = {
                        std::sqrt(output.covariance[areaIndex * dimension + areaIndex]), PropertyStatus::Available};
                }
                output.status = InferenceStatus::available;
                output.observedSupportSufficient = supported;
                if (options.selection == WindowSelection::prespecified)
                {
                    output.groupDecision = supported && output.logAdjustedP <= std::log(options.familyAlpha)
                        ? ConditionalDecision::passed : ConditionalDecision::rejected;
                }
            }
            catch (const std::domain_error &)
            {
                output.status = InferenceStatus::numericalFailure;
                output.covariance.clear();
            }
            catch (const std::runtime_error &)
            {
                output.status = InferenceStatus::numericalFailure;
                output.covariance.clear();
            }
            return output;
        }

        HyperbolicInference inferHyperbolicWindow(
            const std::vector<const float *> &rawTraces, size_t length, size_t center, size_t scale,
            double fixedStep, const HyperbolicInferenceOptions &options)
        {
            return inferWindow(rawTraces, length, center, scale, fixedStep, options, nullptr);
        }

        HyperbolicInference inferHyperbolicWindow(
            const std::vector<const float *> &rawTraces, size_t length, size_t center, size_t scale,
            const HyperbolicModel::Geometry &cache, const HyperbolicInferenceOptions &options)
        {
            return inferWindow(rawTraces, length, center, scale, cache.step(), options, &cache);
        }

        MappedHyperbolicProperties mapHyperbolicAffine(
            const HyperbolicInference &inference,
            const AffinePeakCoordinates &coordinates)
        {
            if (inference.status != InferenceStatus::available || !std::isfinite(coordinates.origin)
                || !std::isfinite(coordinates.unitsPerU) || coordinates.unitsPerU == 0)
            {
                throw std::invalid_argument("Affine mapping requires available inference and finite nonzero scale");
            }
            MappedHyperbolicProperties mapped;
            mapped.axisUnit = coordinates.unit;
            const size_t dimension = 3 * inference.properties.size();
            if (dimension == 0 || inference.covariance.size() / dimension != dimension
                || inference.covariance.size() % dimension != 0)
            {
                throw std::invalid_argument("Invalid inference covariance shape");
            }
            for (const auto &property : inference.properties)
            {
                const std::array<double, 3> values = {{
                    coordinates.origin + coordinates.unitsPerU * property.apexU.value,
                    property.height.value,
                    std::abs(coordinates.unitsPerU) * property.areaU.value}};
                mapped.values.push_back(values);
            }
            mapped.covariance.resize(inference.covariance.size());
            const double factors[3] = {coordinates.unitsPerU, 1, std::abs(coordinates.unitsPerU)};
            for (size_t row = 0; row < dimension; ++row)
            {
                for (size_t column = 0; column < dimension; ++column)
                {
                    mapped.covariance[row * dimension + column] = inference.covariance[row * dimension + column]
                        * factors[row % 3] * factors[column % 3];
                    if (!std::isfinite(mapped.covariance[row * dimension + column]))
                    {
                        throw std::overflow_error("Affine covariance outside numeric range");
                    }
                }
                if (!std::isfinite(mapped.values[row / 3][row % 3]))
                {
                    throw std::overflow_error("Affine property outside numeric range");
                }
            }
            return mapped;
        }
    }
}
