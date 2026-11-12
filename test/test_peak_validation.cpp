#include "qalgorithms_peak_validation.h"

#include <cmath>
#include <iostream>
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

template<typename Curve>
static double midpointIntegral(const Curve &curve, double left, double right)
{
    constexpr size_t points = 60000;
    const double step = (right - left) / points;
    double sum = 0;
    for (size_t i = 0; i < points; ++i)
    {
        sum += curve(left + (static_cast<double>(i) + .5) * step);
    }
    return sum * step;
}

static void testWindowCoverage()
{
    using Quadratic = Model<PeakModelKind::PiecewiseQuadratic>;
    const PeakFormNoise noise{.01, 5};
    RegCoeffs gaussian;
    gaussian.b2 = gaussian.b3 = -.05;
    gaussian.scale = 4;
    const auto result = validatePeakWindowCoverage(gaussian, noise, .7);
    const double fraction = std::erf(4 * std::sqrt(.05));
    const double derivative = -2 * std::exp(-16 * .05) / std::sqrt(std::acos(-1.0) * .05);
    const double jacobian[4] = {0, 0, derivative, derivative};
    const double deviation = std::sqrt(.01 * Quadratic::covarianceProjection(jacobian, gaussian.scale));
    require(std::abs(result.fraction - fraction) < 1e-9 && std::abs(result.standardError - deviation) < 1e-7
        && std::abs(result.lowerConfidenceBound - (fraction - 2.015048373 * deviation)) < 1e-6,
        "Gaussian coverage, analytic gradient and one-sided t bound");
    require(result.status == PeakValidationStatus::passed
        && validatePeakWindowCoverage(gaussian, noise, .9).status == PeakValidationStatus::failed,
        "Minimum coverage is a caller policy");
    gaussian.b0 = 1e3;
    require(validatePeakWindowCoverage(gaussian, noise, .7).fraction == result.fraction,
        "Amplitude cancels even when its exponential is unrepresentable");
    require(validatePeakWindowCoverage(gaussian, {.01, 6}, .7).status == PeakValidationStatus::unavailable,
        "Coverage uncertainty requires matching residual degrees of freedom");
    require(validatePeakCandidate(gaussian, {1e-4, 5}, 2, .7).windowCoverage.status == PeakValidationStatus::passed,
        "The common cascade applies an explicitly requested coverage threshold");
    require(validatePeakCandidate(gaussian, {1, 5}, 2, .7).windowCoverage.status == PeakValidationStatus::notEvaluated,
        "Coverage is skipped when shape confidence fails");

    RegCoeffs mixed;
    mixed.b1 = .4;
    mixed.b2 = .2;
    mixed.b3 = -.15;
    mixed.scale = 4;
    const auto clipped = validatePeakWindowCoverage(mixed, noise, .5);
    const auto mixedCurve = [&](double x) { return std::exp(Quadratic::logValue(&mixed, x)); };
    const double expectedClipped = midpointIntegral(mixedCurve, -1, 4) / midpointIntegral(mixedCurve, -1, 50);
    require(std::abs(clipped.fraction - expectedClipped) < 1e-7,
        "Quadratic coverage excludes the rising branch beyond the valley");
    mixed.b2 = .05;
    const auto atValley = validatePeakWindowCoverage(mixed, noise, .5);
    require(std::isfinite(atValley.fraction) && atValley.status == PeakValidationStatus::unavailable,
        "A valley at the window edge has no unique first-order derivative");
    mixed.b2 = 0;
    const auto linearFlank = validatePeakWindowCoverage(mixed, noise, .5);
    const double rightArea = midpointIntegral(mixedCurve, 0, 50);
    const double expectedLinear = ((1 - std::exp(-4 * mixed.b1)) / mixed.b1
        + midpointIntegral(mixedCurve, 0, 4)) / (1 / mixed.b1 + rightArea);
    require(std::abs(linearFlank.fraction - expectedLinear) < 1e-8,
        "A decaying zero-curvature flank has the exponential area limit");

    const HyperbolicModel::Geometry geometry(4, .35);
    HyperbolicModel::Fit fit;
    fit.scale = 4;
    fit.identifiable = true;
    fit.sharedSlope = .2;
    fit.locals = {{0, 1.2, .3}, {0, 1.5, -.1}};
    const auto hyperbolic = validatePeakWindowCoverage(fit, 0, geometry, {.01, 11}, .7);
    const auto referenceFraction = [](double b, double c, double d)
    {
        const auto curve = [&](double u) { return std::exp(b * u - c * std::cosh(u) - d * std::sinh(u)); };
        return midpointIntegral(curve, -1.4, 1.4) / midpointIntegral(curve, -12, 12);
    };
    constexpr double delta = 1e-4;
    const double slope = (referenceFraction(.2 + delta, 1.2, .3) - referenceFraction(.2 - delta, 1.2, .3)) / (2 * delta);
    const double curvature = (referenceFraction(.2, 1.2 + delta, .3) - referenceFraction(.2, 1.2 - delta, .3)) / (2 * delta);
    const double skew = (referenceFraction(.2, 1.2, .3 + delta) - referenceFraction(.2, 1.2, .3 - delta)) / (2 * delta);
    const auto &g = geometry.at(4);
    const double covariance = slope * slope * hyperbolicCovarianceFactor(g, 2, 0, 0)
        + curvature * curvature * hyperbolicCovarianceFactor(g, 2, 2, 2)
        + skew * skew * hyperbolicCovarianceFactor(g, 2, 3, 3)
        + 2 * slope * skew * hyperbolicCovarianceFactor(g, 2, 0, 3);
    require(std::abs(hyperbolic.fraction - referenceFraction(.2, 1.2, .3)) < 1e-8
        && std::abs(hyperbolic.standardError - std::sqrt(.01 * covariance)) < 1e-7,
        "Hyperbolic coverage and pooled coefficient covariance match independent integration");
    fit.locals.front().skew = 1.2;
    require(validatePeakWindowCoverage(fit, 0, geometry, {.01, 11}, .7).status == PeakValidationStatus::unavailable,
        "A nonintegrable hyperbolic fit has no coverage decision");
}

static void testAreaEvidence()
{
    using Quadratic = Model<PeakModelKind::PiecewiseQuadratic>;
    const PeakFormNoise noise{1e-6, 1};
    const double critical = std::tan((.5 - noise.alpha / 8) * std::acos(-1.0));
    RegCoeffs quadratic;
    quadratic.b0 = 5;
    quadratic.b2 = quadratic.b3 = -.3;
    quadratic.scale = 2;
    const auto q = validatePeakAreaEvidence(quadratic, noise, 100);
    const auto radius = [&](size_t i)
    {
        double gradient[4] = {};
        gradient[i] = 1;
        return critical * std::sqrt(noise.logMse * Quadratic::covarianceProjection(gradient, quadratic.scale));
    };
    const double intercept = quadratic.b0 - radius(0), slope = radius(1), curve = quadratic.b2 - radius(2);
    const double qEnvelope = midpointIntegral([&](double x)
    {
        return std::exp(intercept - slope * std::abs(x) + curve * x * x);
    }, -2, 2);
    const double qBound = std::exp(q.logLowerBound);
    require(q.status == PeakValidationStatus::passed && q.shapeCertified && q.left == -2 && q.right == 2
        && qBound < qEnvelope && qBound / qEnvelope > .999, "Quadratic lower bound matches an independently integrated confidence envelope");
    require(validatePeakAreaEvidence(quadratic, noise, 1000).status == PeakValidationStatus::failed,
        "Area below a larger practical minimum is not established");
    const auto gated = validatePeakCandidate(quadratic, noise, 2, .5, 1000);
    require(gated.areaEvidence.status == PeakValidationStatus::failed
        && gated.windowCoverage.status == PeakValidationStatus::notEvaluated, "Area evidence precedes coverage in the common cascade");
    quadratic.b0 += 1000;
    require(std::abs(validatePeakAreaEvidence(quadratic, noise, 100).logLowerBound - q.logLowerBound - 1000) < 1e-12,
        "Log area stays usable when the raw amplitude overflows");

    RegCoeffs mixed;
    mixed.b0 = 5;
    mixed.b1 = .4;
    mixed.b2 = .2;
    mixed.b3 = -.15;
    mixed.scale = 4;
    const auto clipped = validatePeakAreaEvidence(mixed, {1e-7, 5}, 50);
    require(clipped.status == PeakValidationStatus::passed && clipped.left > -1 && clipped.left < 0 && clipped.right == 4,
        "The certified mixed-tail interval ends before every plausible quadratic valley");

    const HyperbolicModel::Geometry geometry(2, .7);
    HyperbolicModel::Fit fit;
    fit.scale = 2;
    fit.identifiable = true;
    fit.sharedSlope = .2;
    fit.locals = {{6, 1.2, .1}};
    const auto h = validatePeakAreaEvidence(fit, 0, geometry, noise, 50);
    const auto &g = geometry.at(2);
    const auto hyperRadius = [&](size_t i)
    {
        return critical * std::sqrt(noise.logMse * hyperbolicCovarianceFactor(g, 1, i, i));
    };
    const double hEnvelope = midpointIntegral([&](double u)
    {
        return std::exp(6 - hyperRadius(1) + .2 * u - hyperRadius(0) * std::abs(u)
            - (1.2 + hyperRadius(2)) * std::cosh(u) - .1 * std::sinh(u) - hyperRadius(3) * std::abs(std::sinh(u)));
    }, -1.4, 1.4);
    const double hBound = std::exp(h.logLowerBound);
    require(h.status == PeakValidationStatus::passed && h.shapeCertified && hBound < hEnvelope && hBound / hEnvelope > .999,
        "Hyperbolic lower bound matches an independently integrated confidence envelope");
    fit.locals.push_back({8, 1.5, -.25});
    const auto pooled = validatePeakAreaEvidence(fit, 1, geometry, {1e-6, 3}, 50);
    const auto pooledRadius = [&](size_t i)
    {
        constexpr double quantile = 5.391949071934051; // independent Student t survival inverse, DF=3
        return quantile * std::sqrt(noise.logMse * hyperbolicCovarianceFactor(g, 2, i, i));
    };
    const double pooledEnvelope = midpointIntegral([&](double u)
    {
        return std::exp(8 - pooledRadius(4) + .2 * u - pooledRadius(0) * std::abs(u)
            - (1.5 + pooledRadius(5)) * std::cosh(u) + .25 * std::sinh(u) - pooledRadius(6) * std::abs(std::sinh(u)));
    }, -1.4, 1.4);
    const double pooledBound = std::exp(pooled.logLowerBound);
    require(pooled.status == PeakValidationStatus::passed && pooledBound < pooledEnvelope && pooledBound / pooledEnvelope > .999,
        "Pooled area evidence uses the selected trace and shared slope covariance");
    fit.locals.pop_back();
    fit.locals.front().skew = 1.199;
    const auto uncertain = validatePeakAreaEvidence(fit, 0, geometry, noise, 1);
    require(uncertain.status == PeakValidationStatus::failed && !uncertain.shapeCertified
        && uncertain.logLowerBound < 0 && std::isinf(uncertain.logLowerBound),
        "An uncertain tail cannot certify area even with a positive fitted shape");
    require(validatePeakAreaEvidence(fit, 0, geometry, {1e-6, 2}, 1).status == PeakValidationStatus::unavailable,
        "Area inference requires degrees of freedom from the same fit");
}

int main()
{
    try
    {
        testWindowCoverage();
        testAreaEvidence();
        RegCoeffs quadratic;
        quadratic.b0 = std::numeric_limits<double>::quiet_NaN();
        quadratic.b2 = -2;
        quadratic.b3 = -3;
        quadratic.scale = 4;
        const auto centered = validatePeakCandidate(quadratic, {.01, 5});
        require(centered.geometry.form == PeakValidationStatus::passed
            && centered.geometry.apexSupport == PeakValidationStatus::passed
            && centered.geometry.apex == 0 && centered.formConfidence == PeakValidationStatus::passed,
            "Quadratic shape and apex ignore the intercept");
        quadratic.b2 = 1;
        require(!hasPeakShape(quadratic), "Mixed flank at a stationary join has no strict apex");
        quadratic.b1 = 1;
        require(hasPeakShape(quadratic) && validatePeakFormConfidence(quadratic, 1e-4, 5) == PeakValidationStatus::passed,
            "A supported slope sign can establish the clipped quadratic branch");
        quadratic.b1 = 18;
        const auto boundary = validatePeakCandidate(quadratic, {.01, 5});
        require(boundary.geometry.form == PeakValidationStatus::passed
            && boundary.geometry.apexSupport == PeakValidationStatus::failed
            && boundary.formConfidence == PeakValidationStatus::notEvaluated,
            "Apex at the second sample from the edge leaves one bracketing sample");

        const HyperbolicModel::Geometry cache(2, .35);
        HyperbolicModel::Fit fit;
        fit.scale = 2;
        fit.identifiable = true;
        fit.locals = {{std::numeric_limits<double>::quiet_NaN(), 1, 1}};
        const auto invalid = validatePeakCandidate(fit, 0, cache, {.01, 1});
        require(invalid.geometry.form == PeakValidationStatus::failed
            && invalid.geometry.apexSupport == PeakValidationStatus::notEvaluated,
            "A zero hyperbolic flank fails before apex evaluation");
        fit.locals.front().skew = 0;
        const auto &g = cache.at(2);
        const double deviation = std::sqrt(g.varCurvatureFactor() + g.varLocalSkewFactor(1));
        // Student t with one degree of freedom is a standard Cauchy.
        const double critical = std::tan(.45 * std::acos(-1.0));
        fit.locals.front().curvature = .99 * critical * deviation;
        require(validatePeakCandidate(fit, 0, cache, {1, 1}).formConfidence == PeakValidationStatus::failed,
            "Hyperbolic contrast below the one-sided Cauchy critical value");
        fit.locals.front().curvature = 1.01 * critical * deviation;
        require(validatePeakCandidate(fit, 0, cache, {1, 1}).formConfidence == PeakValidationStatus::passed,
            "Both contrasts pass the one-sided Cauchy critical value");
        fit.locals.front().skew = .5 * fit.locals.front().curvature;
        require(validatePeakFormConfidence(fit, 0, cache, 1, 1) == PeakValidationStatus::failed,
            "One uncertain flank fails the joint shape claim");
        require(validatePeakFormConfidence(fit, 0, cache, 0, 1) == PeakValidationStatus::unavailable,
            "Zero residual variance supplies no noise estimate");
        require(validatePeakFormConfidence(fit, 0, cache, 1, 5) == PeakValidationStatus::unavailable,
            "Degrees of freedom must match the fitted window");
        fit.locals.front().skew = 0;
        fit.locals.front().curvature = 1e155 * deviation;
        require(validatePeakFormConfidence(fit, 0, cache, 1, 1, 1e-200) == PeakValidationStatus::failed,
            "A representable Cauchy tail remains above extremely small alpha");
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
