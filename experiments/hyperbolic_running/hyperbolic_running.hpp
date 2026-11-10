#pragma once

// Research-only, aligned-grid, log-hyperbolic peak regression.
// No production qPeaks include, state change, or numerical fit dependency.
// Model: y_j(u) = a_j + b*u - C_j*cosh(u) - D_j*sinh(u).
// Units: u = step*(sample_index - window_centre); step is fixed externally.
// A common b is estimated by pooled OLS, local (a,C,D) for every trace.
// Unconstrained linear fits are retained for diagnostics; a finite full-area
// interpretation additionally requires (C+D)/2 > 0 and (C-D)/2 > 0.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hyperbolic_research {

struct WindowGeometry {
    std::size_t scale = 0;
    double edgeU = 0, edgeCosh = 0, edgeSinh = 0;
    double invN = 0, meanCosh = 0;
    double invCenteredCoshNorm2 = 0, invSinhNorm2 = 0;
    double kappa = 0, invOddResidualNorm2 = 0;
    double relativeEvenInformation = 0, relativeOddInformation = 0;
    bool identifiable = false;

    // Unit-variance OLS covariance factors. Multiply by pooled residual MSE.
    double varInterceptFactor() const {
        return invN + meanCosh*meanCosh*invCenteredCoshNorm2;
    }
    double varCurvatureFactor() const { return invCenteredCoshNorm2; }
    double covInterceptCurvatureFactor() const {
        return meanCosh*invCenteredCoshNorm2;
    }
    double varSharedSlopeFactor(std::size_t m) const {
        return invOddResidualNorm2/static_cast<double>(m);
    }
    double varLocalSkewFactor(std::size_t m) const {
        return invSinhNorm2 + kappa*kappa*varSharedSlopeFactor(m);
    }
    double covSharedSlopeLocalSkewFactor(std::size_t m) const {
        return kappa*varSharedSlopeFactor(m);
    }
};

class GeometryCache {
public:
    GeometryCache(std::size_t maxScale, double step) : step_(step), geometries_(maxScale+1) {
        if (maxScale < 2 || !(std::isfinite(step) && step > 0.0)
            || !std::isfinite(step*static_cast<double>(maxScale))
            || step*static_cast<double>(maxScale) > 30.0)
            throw std::invalid_argument("Invalid step or geometry cache range");
        for (std::size_t h = 2; h <= maxScale; ++h) {
            auto &g = geometries_[h];
            g.scale = h;
            g.edgeU = step*static_cast<double>(h);
            g.edgeCosh = std::cosh(g.edgeU);
            g.edgeSinh = std::sinh(g.edgeU);
            const std::size_t n = 2*h+1;
            g.invN = 1.0/static_cast<double>(n);
            long double cs=0, ss=0, us=0, uu=0, cc=0;
            for (std::size_t i=0; i<n; ++i) {
                const long double u = static_cast<long double>(step) *
                      (static_cast<long double>(i) - static_cast<long double>(h));
                const long double c=std::cosh(u), s=std::sinh(u);
                cs+=c; ss+=s*s; us+=u*s; uu+=u*u; cc+=c*c;
            }
            const long double cm=cs/static_cast<long double>(n);
            const long double k=us/ss;
            // Stable direct residual norms, not differences of nearly equal sums.
            long double varC=0, residual=0;
            for (std::size_t i=0; i<n; ++i) {
                const long double u = static_cast<long double>(step) *
                      (static_cast<long double>(i) - static_cast<long double>(h));
                const long double dC=std::cosh(u)-cm;
                const long double r=u-k*std::sinh(u);
                varC+=dC*dC;
                residual+=r*r;
            }
            if (!(varC > 0 && residual > 0 && ss > 0))
                throw std::invalid_argument("Unidentifiable hyperbolic regression geometry");
            g.meanCosh=static_cast<double>(cm);
            g.invCenteredCoshNorm2=1.0/static_cast<double>(varC);
            g.invSinhNorm2=1.0/static_cast<double>(ss);
            g.kappa=static_cast<double>(k);
            g.invOddResidualNorm2=1.0/static_cast<double>(residual);
            g.relativeEvenInformation=static_cast<double>(varC/cc);
            g.relativeOddInformation=static_cast<double>(residual/uu);
            g.identifiable=(g.relativeEvenInformation > 1e-12 &&
                            g.relativeOddInformation > 1e-12);
        }
    }
    const WindowGeometry &at(std::size_t scale) const {
        if (scale < 2 || scale >= geometries_.size())
            throw std::out_of_range("Window scale outside prepared cache");
        return geometries_[scale];
    }
    double step() const { return step_; }
    std::size_t maxScale() const { return geometries_.size()-1; }

private:
    double step_;
    std::vector<WindowGeometry> geometries_;
};

struct Moments {
    double sumY=0, sumCoshY=0, sumSinhY=0, sumUY=0;
};

struct LocalFit {
    double intercept=0, curvature=0, skew=0;
};

struct FitView {
    std::size_t center=0, scale=0;
    double sharedSlope=0;
    bool identifiable=false;
    std::vector<LocalFit> locals; // re-used between visitor calls
};

using FitVisitor = std::function<void(const FitView&)>;

inline void forEachRunning(const std::vector<const float*> &traces,
                           std::size_t length, const GeometryCache &cache,
                           const FitVisitor &visitor) {
    if (traces.empty() || length < 5 || !visitor ||
        std::any_of(traces.begin(),traces.end(),[](const float *p){return p==nullptr;}))
        throw std::invalid_argument("Expected nonempty aligned log-intensity traces and a visitor");
    const std::size_t m=traces.size();
    const double invM=1.0/static_cast<double>(m);
    const std::size_t maxH=std::min(cache.maxScale(), (length-1)/2);
    const double step=cache.step();
    const double c1=std::cosh(step), s1=std::sinh(step);
    std::vector<Moments> moments(m);
    FitView fit;
    fit.locals.resize(m);
    for (std::size_t center=2; center+2<length; ++center) {
        double totalUY=0, totalSinhY=0;
        for (std::size_t j=0; j<m; ++j) {
            const float *y=traces[j]+center;
            const double l=y[-1], mid=y[0], r=y[1];
            const double d=r-l;
            moments[j]={l+mid+r,c1*(l+r)+mid,s1*d,step*d};
            totalUY+=step*d;
            totalSinhY+=s1*d;
        }
        const std::size_t localMaxH=std::min(maxH,std::min(center,length-center-1));
        for (std::size_t h=2; h<=localMaxH; ++h) {
            const auto &g=cache.at(h);
            for (std::size_t j=0; j<m; ++j) {
                const double l=traces[j][center-h], r=traces[j][center+h];
                const double s=r+l, d=r-l;
                auto &v=moments[j];
                v.sumY+=s;
                v.sumCoshY+=g.edgeCosh*s;
                v.sumSinhY+=g.edgeSinh*d;
                v.sumUY+=g.edgeU*d;
                totalUY+=g.edgeU*d;
                totalSinhY+=g.edgeSinh*d;
            }
            fit.sharedSlope=(totalUY-g.kappa*totalSinhY)*g.invOddResidualNorm2*invM;
            fit.identifiable=g.identifiable;
            fit.center=center;
            fit.scale=h;
            for (std::size_t j=0; j<m; ++j) {
                const auto &v=moments[j];
                auto &p=fit.locals[j];
                p.curvature=-(v.sumCoshY-g.meanCosh*v.sumY)*g.invCenteredCoshNorm2;
                p.intercept=v.sumY*g.invN+p.curvature*g.meanCosh;
                p.skew=g.kappa*fit.sharedSlope-v.sumSinhY*g.invSinhNorm2;
            }
            visitor(fit);
        }
    }
}

inline bool validIntegrablePeak(const LocalFit &p) {
    return std::isfinite(p.curvature) && std::isfinite(p.skew) &&
           p.curvature > std::abs(p.skew);
}

inline double apexU(const LocalFit &peak, double sharedSlope) {
    if (!validIntegrablePeak(peak))
        return std::numeric_limits<double>::quiet_NaN();
    const double p=(peak.curvature+peak.skew)*0.5;
    const double q=(peak.curvature-peak.skew)*0.5;
    return std::asinh(sharedSlope/(2*std::sqrt(p*q))) + 0.5*std::log(q/p);
}

inline double logAreaU(const LocalFit &peak, double sharedSlope) {
    if (!validIntegrablePeak(peak))
        return std::numeric_limits<double>::quiet_NaN();
    const double p=(peak.curvature+peak.skew)*0.5;
    const double q=(peak.curvature-peak.skew)*0.5;
    const double k=std::cyl_bessel_k(std::abs(sharedSlope),2*std::sqrt(p*q));
    if (!(k > 0.0 && std::isfinite(k)))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(2.0)+peak.intercept+0.5*sharedSlope*std::log(q/p)+std::log(k);
}

inline double predict(const LocalFit &p, double b, double u) {
    return p.intercept+b*u-p.curvature*std::cosh(u)-p.skew*std::sinh(u);
}

} // namespace hyperbolic_research
