#include "detail/qalgorithms_peak_background.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

// Raw-domain least squares and rotation calibration. The caller supplies
// the observed template and the same template procedure for null draws.
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
                double xNorm, const std::vector<double> &observedShape, bool supported, const PeakBackgroundReplay &replay,
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
                    if (replay(simulated, shape))
                    {
                        if (shape.size() != simulated.size())
                        {
                            throw std::invalid_argument("Replayed template must match the raw window length");
                        }
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

        PeakBackgroundValidation runPeakBackground(const std::vector<double> &values, const std::vector<double> &x,
            const std::vector<double> &shape, bool supported, const PeakBackgroundReplay &replay, const PeakBackgroundOptions &options)
        {
            if (values.empty() || x.size() != values.size() || (supported && (shape.size() != values.size() || !replay)))
            {
                throw std::invalid_argument("Invalid background window or template");
            }
            PeakBackgroundValidation result;
            result.supportedPeakShape = supported;
            std::mt19937_64 engine(options.seed);
            result.constant = compareBackground(values, x, 0, shape, supported, replay, options.rotations, engine);
            result.linear = compareBackground(values, x, dot(x, x), shape, supported, replay, options.rotations, engine);
            if (std::isfinite(result.constant.pValue) && std::isfinite(result.linear.pValue))
            {
                result.status = std::max(result.constant.pValue, result.linear.pValue) <= options.alpha
                    ? PeakValidationStatus::passed : PeakValidationStatus::failed;
            }
            return result;
        }
    }
}
