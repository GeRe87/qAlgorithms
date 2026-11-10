#include "qalgorithms_peak_inference.h"

#define main researchReferenceTestsC
#include "../experiments/hyperbolic_running/test.cpp"
#undef main

#include <random>
#include "cephes.h"

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

namespace
{
    struct ReferenceProperties
    {
        std::vector<double> values;
        double rss = 0;
    };

    // Independent dense solve, direct scalar raw LS, derivative bisection and
    // Simpson whole-line integration. No production property/gradient calls.
    ReferenceProperties referenceProperties(const std::vector<std::vector<float>> &raw, double step)
    {
        const size_t m = raw.size(), n = raw[0].size(), center = n / 2;
        std::vector<std::vector<float>> logs(m, std::vector<float>(n));
        std::vector<const float *> pointers;
        for (size_t j = 0; j < m; ++j)
        {
            for (size_t i = 0; i < n; ++i)
            {
                logs[j][i] = static_cast<float>(std::log(static_cast<double>(raw[j][i])));
            }
            pointers.push_back(logs[j].data());
        }
        const auto fit = solveReference(pointers, center, center, step);
        const double slope = static_cast<double>(fit.beta[3 * m]);
        ReferenceProperties output;
        for (size_t j = 0; j < m; ++j)
        {
            const double intercept = static_cast<double>(fit.beta[3 * j]);
            const double curvature = static_cast<double>(fit.beta[3 * j + 1]);
            const double skew = static_cast<double>(fit.beta[3 * j + 2]);
            long double numerator = 0, denominator = 0;
            for (size_t i = 0; i < n; ++i)
            {
                const double u = step * (static_cast<double>(i) - center);
                const double shape = std::exp(slope * u - curvature * std::cosh(u) - skew * std::sinh(u));
                numerator += raw[j][i] * shape;
                denominator += shape * shape;
                const double error = logs[j][i] - intercept - std::log(shape);
                output.rss += error * error;
            }
            const double amplitude = static_cast<double>(numerator / denominator);
            double lo = -12, hi = 12;
            for (int k = 0; k < 100; ++k)
            {
                const double mid = (lo + hi) / 2;
                if (slope - curvature * std::sinh(mid) - skew * std::cosh(mid) > 0)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
            }
            const double apex = (lo + hi) / 2;
            const double height = amplitude * std::exp(slope * apex - curvature * std::cosh(apex) - skew * std::sinh(apex));
            constexpr int intervals = 12000;
            constexpr double du = 16.0 / intervals;
            double integral = 0;
            for (int k = 0; k <= intervals; ++k)
            {
                const double u = -8 + k * du;
                integral += (k == 0 || k == intervals ? 1 : (k % 2 ? 4 : 2))
                    * amplitude * std::exp(slope * u - curvature * std::cosh(u) - skew * std::sinh(u));
            }
            output.values.insert(output.values.end(), {apex, height, integral * du / 3});
        }
        return output;
    }

    void testFullInfluenceCovariance()
    {
        for (size_t m : {1u, 3u})
        {
            constexpr size_t n = 17, scale = 8;
            constexpr double step = .18;
            std::vector<std::vector<float>> raw(m, std::vector<float>(n));
            std::vector<const float *> pointers;
            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i = 0; i < n; ++i)
                {
                    const double u = step * (static_cast<double>(i) - scale);
                    raw[j][i] = static_cast<float>(std::exp(7 + .1 * j + .6 * u - (2 + .1 * j) * std::cosh(u)
                        - .2 * std::sinh(u) + .015 * std::sin(i * 1.4 + j)));
                }
                pointers.push_back(raw[j].data());
            }
            HyperbolicInferenceOptions options;
            options.selection = WindowSelection::prespecified;
            options.familySize = 9;
            const auto actual = inferHyperbolicWindow(pointers, n, scale, scale, step, options);
            require(actual.status == InferenceStatus::available, "Available selected-window inference");
            require(actual.groupDecision == ConditionalDecision::passed, "Prespecified curved group passes conditional gate");
            require(actual.residualDegreesOfFreedom == m * (n - 3) - 1, "Hyperbolic residual DF");
            require(actual.comparisonDegreesOfFreedom == 2 * m, "Hyperbolic comparison DF");
            near(actual.adjustedP, std::min(1.0, 9 * actual.nominalP), 1e-14, "Explicit Bonferroni family adjustment");
            require(actual.adjustedP > 0, "Nonzero p-value exercises significance threshold");
            options.familyAlpha = actual.adjustedP / 2;
            const auto insignificant = inferHyperbolicWindow(pointers, n, scale, scale, step, options);
            require(insignificant.groupDecision == ConditionalDecision::rejected, "Declared significance threshold gates comparison");
            options.familyAlpha = .05;
            const auto reference = referenceProperties(raw, step);
            near(actual.residualSumSquares, reference.rss, 1e-9, "Full-window RSS vs dense solve");
            const size_t dimension = 3 * m, observations = n * m;
            std::vector<std::vector<double>> jacobian(dimension, std::vector<double>(observations));
            constexpr double epsilon = .002;
            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i = 0; i < n; ++i)
                {
                    auto plus = raw, minus = raw;
                    plus[j][i] = static_cast<float>(raw[j][i] * std::exp(epsilon));
                    minus[j][i] = static_cast<float>(raw[j][i] * std::exp(-epsilon));
                    const auto high = referenceProperties(plus, step);
                    const auto low = referenceProperties(minus, step);
                    const double difference = std::log(static_cast<double>(plus[j][i]) / minus[j][i]);
                    for (size_t row = 0; row < dimension; ++row)
                    {
                        jacobian[row][j * n + i] = (high.values[row] - low.values[row]) / difference;
                    }
                }
            }
            std::vector<double> expected(dimension * dimension, 0);
            for (size_t row = 0; row < dimension; ++row)
            {
                for (size_t column = 0; column < dimension; ++column)
                {
                    for (size_t k = 0; k < observations; ++k)
                    {
                        expected[row * dimension + column] += jacobian[row][k] * jacobian[column][k] * actual.logNoiseVariance;
                    }
                }
            }
            for (size_t row = 0; row < dimension; ++row)
            {
                const double value = row % 3 == 0 ? actual.properties[row / 3].apexU.value
                    : (row % 3 == 1 ? actual.properties[row / 3].height.value : actual.properties[row / 3].areaU.value);
                near(value, reference.values[row], 1e-7, "Corrected properties vs independent reference");
                for (size_t column = 0; column < dimension; ++column)
                {
                    const double normalization = std::sqrt(expected[row * dimension + row] * expected[column * dimension + column]);
                    if (std::abs(actual.covariance[row * dimension + column] - expected[row * dimension + column]) > .003 * normalization + 1e-10)
                    {
                        throw std::runtime_error("Full postfit covariance differs from independent observational finite differences");
                    }
                    require(actual.covariance[row * dimension + column] == actual.covariance[column * dimension + row], "Symmetric property covariance");
                }
            }
            if (m > 1)
            {
                require(std::abs(expected[2 * dimension + 5]) > 1e-10, "Nonzero cross-trace corrected-area covariance exercised");
            }
            const auto mapped = mapHyperbolicAffine(actual, {100, -2.5, AxisUnit::MassToCharge});
            for (size_t j = 0; j < m; ++j)
            {
                near(mapped.values[j][0], 100 - 2.5 * actual.properties[j].apexU.value, 0, "Signed affine position");
                near(mapped.values[j][2], 2.5 * actual.properties[j].areaU.value, 0, "Affine area Jacobian");
                near(mapped.covariance[(3 * j) * dimension + 3 * j + 2],
                    -6.25 * actual.covariance[(3 * j) * dimension + 3 * j + 2], 1e-12, "Signed affine cross covariance");
                near(actual.properties[j].areaUncertainty.value,
                    std::sqrt(actual.covariance[(3 * j + 2) * dimension + 3 * j + 2]), 0, "Area standard error");
            }
            const auto exploratory = inferHyperbolicWindow(pointers, n, scale, scale, step);
            require(exploratory.groupDecision == ConditionalDecision::unavailable, "Exploratory windows do not acquire statistical acceptance");
            options.minimumSamplesPerSide = n;
            const auto unsupported = inferHyperbolicWindow(pointers, n, scale, scale, step, options);
            require(unsupported.groupDecision == ConditionalDecision::rejected, "Observed support is a separate gate");
        }
    }

    void testInvalidRequests()
    {
        std::vector<float> raw(17, 100);
        auto reject = [&](const HyperbolicInferenceOptions &options, double step)
        {
            bool caught = false;
            try
            {
                (void)inferHyperbolicWindow({raw.data()}, raw.size(), 8, 8, step, options);
            }
            catch (const std::invalid_argument &)
            {
                caught = true;
            }
            require(caught, "Invalid inference request rejected");
        };
        HyperbolicInferenceOptions options;
        options.familySize = 0;
        reject(options, .18);
        options.familySize = 1;
        options.familyAlpha = 1;
        reject(options, .18);
        options.familyAlpha = .05;
        reject(options, 0);
        raw[8] = 0;
        reject(options, .18);
        raw[8] = -1;
        reject(options, .18);
        raw[8] = std::numeric_limits<float>::infinity();
        reject(options, .18);
        HyperbolicInference unavailable;
        bool caught = false;
        try
        {
            (void)mapHyperbolicAffine(unavailable, {});
        }
        catch (const std::invalid_argument &)
        {
            caught = true;
        }
        require(caught, "Unavailable inference cannot acquire physical properties");
    }

    void testNullCalibrationAndTailProbability()
    {
        std::mt19937 random(18754);
        std::normal_distribution<double> noise(0, .04);
        size_t rejected = 0;
        constexpr size_t trials = 1500, n = 17, scale = 8;
        for (size_t trial = 0; trial < trials; ++trial)
        {
            const size_t m = trial % 3 + 1;
            std::vector<std::vector<float>> raw(m, std::vector<float>(n));
            std::vector<const float *> pointers;
            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i = 0; i < n; ++i)
                {
                    raw[j][i] = static_cast<float>(std::exp(6 + .1 * j + .3 * .18 * (static_cast<double>(i) - scale) + noise(random)));
                }
                pointers.push_back(raw[j].data());
            }
            const auto result = inferHyperbolicWindow(pointers, n, scale, scale, .18);
            require(std::isfinite(result.nominalP), "Nominal pooled test remains diagnostic for invalid peak shapes");
            // Independent finite polynomial I_x(a,m), since numerator DF=2m.
            const double a = result.residualDegreesOfFreedom / 2.0;
            const double x = result.residualDegreesOfFreedom / (result.residualDegreesOfFreedom + 2 * m * result.comparisonF);
            double series = 1, term = 1;
            for (size_t k = 1; k < m; ++k)
            {
                term *= (a + k - 1) * (1 - x) / k;
                series += term;
            }
            near(result.nominalP, std::pow(x, a) * series, 2e-12, "F survival vs ordinary-scale beta polynomial");
            const double continuedFraction = x == 0 ? 0 : (x == 1 ? 1 : cephes::incbeta(a, static_cast<double>(m), x));
            near(result.nominalP, continuedFraction, 2e-7, "F survival vs independent incomplete-beta continued fraction");
            rejected += result.nominalP < .05;
        }
        const double rate = static_cast<double>(rejected) / trials;
        require(rate > .025 && rate < .075, "Fixed-window null Monte Carlo calibration");
        std::cout << "Fixed-window nominal null rejection rate: " << rate << '\n';
    }
}

int main()
{
    try
    {
        testFullInfluenceCovariance();
        testInvalidRequests();
        testNullCalibrationAndTailProbability();
        std::cout << "Stage C independent inference and covariance checks passed\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
