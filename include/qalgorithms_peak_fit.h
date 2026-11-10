#ifndef QALGORITHMS_PEAK_FIT_H
#define QALGORITHMS_PEAK_FIT_H
#include "qalgorithms_qpeaks.h"
#include "detail/qalgorithms_hyperbolic.h"
#include <variant>
#include <limits>

namespace qAlgorithms {
namespace peak_model {
// Mathematical fitting only: never statistical acceptance or physical calibration.
template<> struct Model<PeakModelKind::Hyperbolic> {
    using Coefficients = detail::hyperbolic::LocalFit;
    using Geometry = detail::hyperbolic::GeometryCache;
    using Fit = detail::hyperbolic::FitView;
    static void runningFit(const std::vector<const float*> &logs, size_t n,
                           const Geometry &geometry, const detail::hyperbolic::FitVisitor &visit);
};
using HyperbolicModel = Model<PeakModelKind::Hyperbolic>;

// Unit-variance log-OLS factors, NOT standard errors or postfit covariance.
// Parameter order: [b, a_0,C_0,D_0, a_1,C_1,D_1, ...].
// Homoskedastic independent log observations on the fixed aligned grid only.
double hyperbolicCovarianceFactor(const detail::hyperbolic::WindowGeometry &,
                                 size_t traces, size_t row, size_t column);
enum class PropertyStatus { Unavailable, Available, NonIntegrable, NumericalFailure };
struct FitProperty {
    double value = std::numeric_limits<double>::quiet_NaN();
    PropertyStatus status = PropertyStatus::Unavailable;
};
struct HyperbolicProperties {
    bool finiteCoefficients = false, identifiable = false, integrable = false;
    bool apexInObservedSupport = false;
    FitProperty apexU, height, areaU;
    FitProperty widthU, areaUncertainty, statisticalAcceptance; // unavailable in B
    detail::hyperbolic::AmplitudeRefit amplitude;
    // Corrected intercept is separate; the original Fit is never mutated.
    double correctedIntercept = std::numeric_limits<double>::quiet_NaN();
};
// Evaluate one explicitly selected trace/window, outside the running hot loop.
// nullptr raw means no amplitude refit; otherwise raw has exactly length samples.
// Area is intensity*u; apex/width are dimensionless u. No m/z or RT conversion.
HyperbolicProperties evaluateHyperbolicWindow(const HyperbolicModel::Fit &fit,
    size_t trace, const HyperbolicModel::Geometry &geometry,
    const float *raw, size_t length);
}

struct PeakFitRequest {
    PeakModelKind model; // explicit selection, independent of legacy detection API
    std::vector<const float*> logTraces;
    size_t length = 0;
    size_t maxHalfWidth = 0;
    double fixedStep = 0; // hyperbolic u step; PiecewiseQuadratic requires 1
};
// Borrowed pointers valid only during callback; reused every window. Copy the
// pointed-to fit to retain it. Both paths visit center-major/increasing scale.
// PiecewiseQuadratic uses sample-index coordinates and existing shared-shape fit.
using PeakFitView = std::variant<const SharedShapeRegCoeffs*,
                                const peak_model::HyperbolicModel::Fit*>;
using PeakFitVisitor = std::function<void(const PeakFitView &)>;
// Validate all finite log samples before any callback. No acceptance, postfit,
// Bessel evaluation, or per-window allocation in this streaming diagnostics API.
// qpeaks_find(Hyperbolic) remains unsupported until statistical gates exist.
void forEachPeakFit(const PeakFitRequest &, const PeakFitVisitor &);
}
#endif
