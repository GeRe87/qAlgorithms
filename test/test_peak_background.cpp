#include "qalgorithms_peak_validation.h"

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

static void require(bool condition, const char *message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

static void testBackground(PeakModelKind model, double step)
{
    constexpr size_t length = 17;
    std::vector<float> raw(length);
    PeakBackgroundOptions options;
    options.rotations = 199;
    options.seed = 741;
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 8;
        raw[i] = static_cast<float>(model == PeakModelKind::Hyperbolic
            ? 800 * std::exp(1.8 * (1 - std::cosh(step * x)))
            : 800 * std::exp(-.055 * x * x));
    }
    const auto peak = validatePeakBackground(model, raw.data(), length, step, options);
    const auto repeat = validatePeakBackground(model, raw.data(), length, step, options);
    require(peak.status == PeakValidationStatus::passed && peak.supportedPeakShape
        && peak.constant.explainedFraction > .999999 && peak.linear.explainedFraction > .999999,
        "Both raw-domain tests detect a matched peak and recover its raw RSS improvement");
    require(repeat.constant.pValue == peak.constant.pValue && repeat.linear.pValue == peak.linear.pValue,
        "Explicit seed reproduces both rank tests");
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 8;
        raw[i] += static_cast<float>(10 + .5 * x);
    }
    const auto onTrend = validatePeakBackground(model, raw.data(), length, step, options);
    require(onTrend.status == PeakValidationStatus::passed && onTrend.linear.explainedFraction > .98,
        "A peak on a sloped raw background passes the affine comparison");

    for (size_t i = 0; i < length; ++i)
    {
        raw[i] = 50 + 2 * static_cast<float>(i);
    }
    const auto affine = validatePeakBackground(model, raw.data(), length, step, options);
    require(affine.status == PeakValidationStatus::failed && affine.linear.pValue == 1
        && affine.linear.backgroundRSS == 0 && affine.constant.backgroundRSS == 1632,
        "Independent raw affine and constant RSS on an exact trend");
    raw.front() = 0;
    const auto zero = validatePeakBackground(model, raw.data(), length, step, options);
    require(!zero.supportedPeakShape && zero.constant.pValue == 1 && zero.linear.pValue == 1,
        "A nonpositive observation gives no log-OLS template");
    if (model == PeakModelKind::Hyperbolic)
    {
        require(validatePeakBackground(model, raw.data(), length, 1e-12, options).status == PeakValidationStatus::unavailable,
            "Unidentifiable hyperbolic geometry has no background decision");
    }

    std::mt19937_64 engine(1527);
    std::normal_distribution<double> normal(0, 5);
    options.rotations = 99;
    constexpr size_t trials = 256;
    for (double slope : {0., 2.})
    {
        size_t rejected = 0, componentRejected = 0;
        for (size_t trial = 0; trial < trials; ++trial)
        {
            for (size_t i = 0; i < length; ++i)
            {
                raw[i] = static_cast<float>(20 + slope * (static_cast<double>(i) - 8) + normal(engine));
            }
            options.seed = engine();
            const auto null = validatePeakBackground(model, raw.data(), length, step, options);
            require(std::isfinite(null.constant.pValue) && std::isfinite(null.linear.pValue), "Identifiable null comparisons produce p-values");
            rejected += null.status == PeakValidationStatus::passed;
            componentRejected += (slope == 0 ? null.constant.pValue : null.linear.pValue) <= options.alpha;
        }
        require(rejected <= 26 && componentRejected <= 26,
            "Neither the joint decision nor the relevant null component may strongly over-reject");
        std::cout << (model == PeakModelKind::Hyperbolic ? "Hyperbolic" : "Quadratic")
                  << ", null slope=" << slope << ", rejected=" << rejected << '/' << trials
                  << ", null component rejected=" << componentRejected << '/' << trials << '\n';
    }
}

int main()
{
    try
    {
        testBackground(PeakModelKind::PiecewiseQuadratic, 1);
        testBackground(PeakModelKind::Hyperbolic, .25);
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
