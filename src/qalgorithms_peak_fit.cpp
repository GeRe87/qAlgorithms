#include "qalgorithms_peak_fit.h"
#include <cmath>
#include <stdexcept>

// Built-in fit dispatch and the bridge from hyperbolic coefficients to
// fixed-window fit properties.
namespace qAlgorithms {
namespace {
void validateLogs(const std::vector<const float*> &logs, size_t n) {
    if (logs.empty() || n<5) throw std::invalid_argument("Expected aligned log traces with length >=5");
    for (const auto *p:logs) {
        if (!p) throw std::invalid_argument("Null log trace");
        for (size_t i=0;i<n;++i) if (!std::isfinite(p[i]))
            throw std::invalid_argument("Nonfinite log sample");
    }
}
peak_model::FitProperty finiteProperty(double value) {
    if (std::isfinite(value)) return {value,peak_model::PropertyStatus::Available};
    return {std::numeric_limits<double>::quiet_NaN(),peak_model::PropertyStatus::NumericalFailure};
}
}
namespace peak_model {
void HyperbolicModel::runningFit(const std::vector<const float*> &logs, size_t n,
    const Geometry &geometry, const detail::hyperbolic::FitVisitor &visit) {
    validateLogs(logs,n);
    if (!visit) throw std::invalid_argument("Empty fit visitor");
    detail::hyperbolic::forEachRunning(logs,n,geometry,visit);
}
double hyperbolicCovarianceFactor(const detail::hyperbolic::WindowGeometry &g,
    size_t m, size_t row, size_t col) {
    if (!m || m>(std::numeric_limits<size_t>::max()-1)/3 || row>3*m || col>3*m)
        throw std::invalid_argument("Covariance index/trace count");
    if (row>col) std::swap(row,col);
    if (!row) {
        if (!col) return g.varSharedSlopeFactor(m);
        return (col-1)%3==2 ? g.covSharedSlopeLocalSkewFactor(m) : 0;
    }
    const size_t j=(row-1)/3,k=(col-1)/3, a=(row-1)%3,b=(col-1)%3;
    if(a==2 && b==2) return j==k ? g.varLocalSkewFactor(m) : g.kappa*g.kappa*g.varSharedSlopeFactor(m);
    if(j!=k) return 0;
    if(a==0 && b==0) return g.varInterceptFactor();
    if(a==0 && b==1) return g.covInterceptCurvatureFactor();
    if(a==1 && b==1) return g.varCurvatureFactor();
    return 0;
}
HyperbolicProperties evaluateHyperbolicWindow(const HyperbolicModel::Fit &fit,
    size_t trace, const HyperbolicModel::Geometry &geometry, const float *raw, size_t length) {
    if (trace>=fit.locals.size() || fit.center>=length || fit.scale>fit.center ||
        fit.scale>length-1-fit.center) throw std::invalid_argument("Invalid selected fit window");
    const auto &g=geometry.at(fit.scale);
    HyperbolicProperties out;
    const auto &original=fit.locals[trace];
    out.identifiable=g.identifiable && fit.identifiable;
    out.finiteCoefficients=std::isfinite(original.intercept) && std::isfinite(original.curvature) &&
        std::isfinite(original.skew) && std::isfinite(fit.sharedSlope);
    out.integrable=out.finiteCoefficients && detail::hyperbolic::validIntegrablePeak(original);
    auto corrected=original;
    if (raw) {
        out.amplitude=detail::hyperbolic::refitRawAmplitude(raw,length,fit.center,g,
            geometry.step(),original,fit.sharedSlope,fit.locals.size());
        if(out.amplitude.valid) {
            corrected.intercept=out.amplitude.logAmplitude;
            out.correctedIntercept=corrected.intercept;
        }
    }
    if (!out.integrable) {
        out.apexU.status=out.height.status=out.areaU.status=PropertyStatus::NonIntegrable;
        return out;
    }
    out.apexU=finiteProperty(detail::hyperbolic::apexU(original,fit.sharedSlope));
    out.apexInObservedSupport=out.apexU.status==PropertyStatus::Available &&
        std::abs(out.apexU.value)<=g.edgeU;
    // An invalid requested postfit must not silently expose uncorrected quantities.
    if(raw && !out.amplitude.valid) return out;
    if(out.apexU.status==PropertyStatus::Available)
        out.height=finiteProperty(std::exp(detail::hyperbolic::predict(corrected,fit.sharedSlope,out.apexU.value)));
    try { out.areaU=finiteProperty(std::exp(detail::hyperbolic::logAreaU(corrected,fit.sharedSlope))); }
    catch(const std::exception &) { out.areaU.status=PropertyStatus::NumericalFailure; }
    if(out.areaU.status==PropertyStatus::Available && !(out.areaU.value>0))
        out.areaU={std::numeric_limits<double>::quiet_NaN(),PropertyStatus::NumericalFailure};
    if(out.height.status==PropertyStatus::Available && !(out.height.value>0))
        out.height={std::numeric_limits<double>::quiet_NaN(),PropertyStatus::NumericalFailure};
    return out;
}
}
void forEachPeakFit(const PeakFitRequest &request, const PeakFitVisitor &visit) {
    if(request.model!=PeakModelKind::PiecewiseQuadratic && request.model!=PeakModelKind::Hyperbolic)
        throw std::invalid_argument("Unknown fit model");
    if(!visit || request.maxHalfWidth<2 || request.maxHalfWidth>(request.length ? (request.length-1)/2 : 0))
        throw std::invalid_argument("Invalid visitor or window range");
    validateLogs(request.logTraces,request.length);
    switch(request.model) {
    case PeakModelKind::Hyperbolic: {
        peak_model::HyperbolicModel::Geometry geometry(request.maxHalfWidth,request.fixedStep);
        // Input was checked once above; the original allocation-free kernel follows.
        detail::hyperbolic::forEachRunning(request.logTraces,request.length,geometry,
            [&](const auto &fit) { visit(PeakFitView{&fit}); });
        return;
    }
    case PeakModelKind::PiecewiseQuadratic:
        if(request.fixedStep!=1 || request.maxHalfWidth>maxscale_global)
            throw std::invalid_argument("Legacy diagnostics require sample-index step 1 and cached scales");
        forEachCoefficientsMultiTrace(request.logTraces,request.length,[&](const auto &fit) {
            if(fit.scale<=request.maxHalfWidth) visit(PeakFitView{&fit});
        });
        return;
    }
}
}
