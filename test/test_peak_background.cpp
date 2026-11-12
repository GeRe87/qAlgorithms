#include "qalgorithms_peak_validation.h"
#include "detail/qalgorithms_peak_background.h"

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
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

static void testBackgroundReplay()
{
    const std::vector<double> grid{-3, -2, -1, 0, 1, 2, 3};
    const std::vector<double> shape{0, .1, .4, 1, .4, .1, 0};
    std::vector<double> observed(grid.size());
    for (size_t i = 0; i < grid.size(); ++i)
    {
        observed[i] = 20 + grid[i] + 30 * shape[i];
    }
    PeakBackgroundOptions options;
    options.rotations = 99;
    size_t replayed = 0;
    const auto result = runPeakBackground(observed, grid, shape, true,
        [&](const auto &draw, auto &replayedShape)
        {
            require(draw.size() == grid.size(), "Replay receives one complete raw window");
            ++replayed;
            replayedShape = shape;
            return true;
        }, options);
    require(replayed == 2 * options.rotations && result.linear.explainedFraction > .999999,
        "The raw comparison accepts a model-independent template procedure");
    bool wrongLengthRejected = false;
    try
    {
        (void)runPeakBackground(observed, grid, shape, true, [](const auto &, auto &replayedShape)
        {
            replayedShape.clear();
            return true;
        }, options);
    }
    catch (const std::invalid_argument &)
    {
        wrongLengthRejected = true;
    }
    require(wrongLengthRejected, "A successful replay must return one value per raw sample");
}

static void testCompleteCascade(PeakModelKind model, double step)
{
    constexpr size_t length = 17;
    std::vector<float> raw(length);
    FullPeakCandidateOptions options{100, .75, .05, {}};
    options.background.rotations = 199;
    options.background.seed = 741;
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 8;
        const double u = step * x;
        const double shape = model == PeakModelKind::Hyperbolic ? .9 * (1 - std::cosh(u)) + .07 * u - .12 * std::sinh(u)
            : -.055 * x * x;
        raw[i] = static_cast<float>(800 * std::exp(shape + .003 * std::sin(2 * x)));
    }
    const auto evaluate = [&]() { return validatePeakCandidate(model, raw.data(), raw.size(), step, options); };
    bool missingThresholdRejected = false;
    try
    {
        FullPeakCandidateOptions missing;
        missing.minimumArea = 100;
        (void)validatePeakCandidate(model, raw.data(), raw.size(), step, missing);
    }
    catch (const std::invalid_argument &)
    {
        missingThresholdRejected = true;
    }
    require(missingThresholdRejected, "Full validation requires both explicit thresholds");
    const auto complete = evaluate();
    const auto background = validatePeakBackground(model, raw.data(), raw.size(), step, options.background);
    require(complete.status == PeakValidationStatus::passed && complete.stages.geometry.form == PeakValidationStatus::passed
        && complete.stages.geometry.apexSupport == PeakValidationStatus::passed
        && complete.stages.formConfidence == PeakValidationStatus::passed
        && complete.stages.areaEvidence.status == PeakValidationStatus::passed
        && complete.stages.windowCoverage.status == PeakValidationStatus::passed
        && complete.background.status == PeakValidationStatus::passed && complete.logNoise.residualDf == 13,
        "Both models pass the complete prespecified cascade on a supported peak");
    require(complete.background.constant.pValue == background.constant.pValue
        && complete.background.linear.pValue == background.linear.pValue
        && complete.background.linear.explainedFraction == background.linear.explainedFraction,
        "The reused observed fit yields the same raw comparison as the standalone test");
    // Orthogonalize the four design columns independently of either
    // running kernel and project the observed log trace onto their span.
    std::array<std::array<double, length>, 4> basis{};
    std::array<double, length> residual{};
    const auto dot = [](const auto &a, const auto &b)
    {
        return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
    };
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 8, u = step * x;
        basis[0][i] = 1;
        basis[1][i] = u;
        basis[2][i] = model == PeakModelKind::Hyperbolic ? std::cosh(u) : (x < 0 ? x * x : 0);
        basis[3][i] = model == PeakModelKind::Hyperbolic ? std::sinh(u) : (x > 0 ? x * x : 0);
        residual[i] = static_cast<float>(std::log(static_cast<double>(raw[i])));
    }
    for (size_t col = 0; col < basis.size(); ++col)
    {
        auto &column = basis[col];
        for (size_t prev = 0; prev < col; ++prev)
        {
            const double projection = dot(column, basis[prev]);
            for (size_t i = 0; i < length; ++i)
            {
                column[i] -= projection * basis[prev][i];
            }
        }
        const double norm = std::sqrt(dot(column, column));
        for (double &value : column)
        {
            value /= norm;
        }
        const double projection = dot(residual, column);
        for (size_t i = 0; i < length; ++i)
        {
            residual[i] -= projection * column[i];
        }
    }
    require(std::abs(complete.logNoise.logMse - dot(residual, residual) / 13) < 1e-10,
        "Log residual variance agrees with an independent full-window projection");

    options.logAlpha = 1e-80;
    const auto uncertain = evaluate();
    require(uncertain.status == PeakValidationStatus::failed && uncertain.stages.geometry.apexSupport == PeakValidationStatus::passed
        && uncertain.stages.formConfidence == PeakValidationStatus::failed
        && uncertain.stages.areaEvidence.status == PeakValidationStatus::notEvaluated
        && uncertain.background.status == PeakValidationStatus::notEvaluated, "Uncertain shape skips the remaining stages");
    options.logAlpha = .05;
    options.minimumArea = 1e10;
    const auto small = evaluate();
    require(small.status == PeakValidationStatus::failed && small.stages.areaEvidence.status == PeakValidationStatus::failed
        && small.stages.windowCoverage.status == PeakValidationStatus::notEvaluated
        && small.background.status == PeakValidationStatus::notEvaluated, "Area failure stops before coverage and background");
    options.minimumArea = 100;
    options.minimumCoverage = .999999;
    const auto uncovered = evaluate();
    require(uncovered.status == PeakValidationStatus::failed && uncovered.stages.windowCoverage.status == PeakValidationStatus::failed
        && uncovered.background.status == PeakValidationStatus::notEvaluated, "Coverage failure stops before background");
    options.minimumCoverage = .75;
    options.background.alpha = .001; // below the smallest attainable p for 199 rotations
    const auto rejected = evaluate();
    require(rejected.status == PeakValidationStatus::failed && rejected.stages.windowCoverage.status == PeakValidationStatus::passed
        && rejected.background.status == PeakValidationStatus::failed, "Raw background rejection determines the final verdict");

    raw.front() = 0;
    const auto zero = evaluate();
    require(zero.status == PeakValidationStatus::unavailable && zero.stages.geometry.form == PeakValidationStatus::unavailable
        && zero.background.status == PeakValidationStatus::notEvaluated, "No log fit on a window containing zero");
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 8;
        raw[i] = static_cast<float>(std::exp(4 + .02 * x * x));
    }
    const auto convex = evaluate();
    require(convex.status == PeakValidationStatus::failed && convex.stages.geometry.form == PeakValidationStatus::failed
        && convex.stages.geometry.apexSupport == PeakValidationStatus::notEvaluated
        && convex.stages.formConfidence == PeakValidationStatus::notEvaluated
        && convex.background.status == PeakValidationStatus::notEvaluated, "Shape failure skips all later stages");
    for (size_t i = 0; i < length; ++i)
    {
        const double x = static_cast<double>(i) - 18;
        raw[i] = static_cast<float>(std::exp(7 + (model == PeakModelKind::Hyperbolic
            ? 1.1 * (1 - std::cosh(step * x)) : -.03 * x * x)));
    }
    const auto offWindow = evaluate();
    require(offWindow.status == PeakValidationStatus::failed && offWindow.stages.geometry.form == PeakValidationStatus::passed
        && offWindow.stages.geometry.apexSupport == PeakValidationStatus::failed
        && offWindow.stages.formConfidence == PeakValidationStatus::notEvaluated
        && offWindow.background.status == PeakValidationStatus::notEvaluated, "Apex failure skips statistical tests");
}

int main()
{
    try
    {
        testBackgroundReplay();
        testBackground(PeakModelKind::PiecewiseQuadratic, 1);
        testBackground(PeakModelKind::Hyperbolic, .25);
        testCompleteCascade(PeakModelKind::PiecewiseQuadratic, 1);
        testCompleteCascade(PeakModelKind::Hyperbolic, .25);
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
