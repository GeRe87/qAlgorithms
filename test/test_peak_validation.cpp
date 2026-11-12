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

int main()
{
    try
    {
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
