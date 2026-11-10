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

namespace qAlgorithms::detail::hyperbolic {

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

    
// POST-FIT ONLY, not part of the O(1) running-regression hot loop.
// For fixed (b,C,D) and raw observations I_i >= 0, minimize
//     sum_i (I_i - A*exp(b*u_i-C*cosh(u_i)-D*sinh(u_i)))^2.
// The exact scalar OLS solution is A=sum(I_i*f_i)/sum(f_i*f_i).
// It changes only a_j=log(A), never the linear log-domain shape estimates.
// The log-amplitude variance factor below includes the *dependence* of the
// refit on the first-stage pooled log-OLS shape coefficients. It assumes the
// raw intensities are exp(the same independent, homoskedastic log-observations).
// Multiply by the first-stage log-residual MSE for a first-order variance.
// This is NOT a complete covariance/area-uncertainty calculation.
struct AmplitudeRefit {
    bool valid=false;
    double logAmplitude=std::numeric_limits<double>::quiet_NaN();
    double logCorrection=std::numeric_limits<double>::quiet_NaN();
    double rawResidualSumSquares=std::numeric_limits<double>::quiet_NaN();
    double logAmplitudeVarianceFactor=std::numeric_limits<double>::quiet_NaN();
};

inline AmplitudeRefit refitRawAmplitude(
    const float *rawIntensity, std::size_t length, std::size_t center,
    const WindowGeometry &geometry, double step,
    const LocalFit &peak, double sharedSlope, std::size_t pooledTraceCount)
{
    if (rawIntensity==nullptr || pooledTraceCount==0 || center>=length ||
        geometry.scale<2 || geometry.scale>center ||
        geometry.scale>length-1-center || !(step>0 && std::isfinite(step)))
        throw std::invalid_argument("Invalid raw-amplitude refit window or trace count");

    AmplitudeRefit result;
    if (!geometry.identifiable || !std::isfinite(peak.intercept) ||
        !std::isfinite(peak.curvature) || !std::isfinite(peak.skew) ||
        !std::isfinite(sharedSlope))
        return result;
    const std::size_t half=geometry.scale;
    const std::size_t start=center-half;
    const std::size_t end=center+half+1;

    // Normalize the shape before exponentiation to avoid under/overflow
    // when the unnormalized exponential contains very large coefficients.
    double maxLogShape=-std::numeric_limits<double>::infinity();
    for(std::size_t i=start;i<end;++i) {
        if (!(std::isfinite(rawIntensity[i]) && rawIntensity[i]>=0))
            return result;
        const double u=step*(static_cast<double>(i)-static_cast<double>(center));
        const double v=sharedSlope*u-peak.curvature*std::cosh(u)
                                         -peak.skew*std::sinh(u);
        if(!std::isfinite(v))return result;
        maxLogShape=std::max(maxLogShape,v);
    }

    long double numerator=0., denominator=0.;
    long double nwU=0., nwC=0., nwS=0.;
    long double dwU=0., dwC=0., dwS=0.;
    for(std::size_t i=start;i<end;++i) {
        const double u=step*(static_cast<double>(i)-static_cast<double>(center));
        const double c=std::cosh(u),s=std::sinh(u);
        const double z=sharedSlope*u-peak.curvature*c-peak.skew*s;
        const long double g=std::exp(z-maxLogShape);
        const long double raw=static_cast<long double>(rawIntensity[i]);
        const long double weighted=raw*g, squared=g*g;
        numerator+=weighted;
        denominator+=squared;
        nwU+=weighted*u; nwC+=weighted*c; nwS+=weighted*s;
        dwU+=squared*u;  dwC+=squared*c;  dwS+=squared*s;
    }
    if (!(numerator>0 && denominator>0) ||
        !std::isfinite(numerator) || !std::isfinite(denominator))
        return result;
    const long double fittedA=numerator/denominator;
    const long double fittedLogA=std::log(fittedA)-maxLogShape;
    if (!std::isfinite(fittedLogA) ||
        !std::isfinite(static_cast<double>(fittedLogA)))
        return result;

    // Derivatives of the fitted log-amplitude with respect to shape b,C,D,
    // keeping observed raw intensities fixed.
    const long double gb=nwU/numerator-2*dwU/denominator;
    const long double gc=-nwC/numerator+2*dwC/denominator;
    const long double gd=-nwS/numerator+2*dwS/denominator;
    const long double invVc=geometry.invCenteredCoshNorm2;
    const long double invSs=geometry.invSinhNorm2;
    const long double R=1.0L/geometry.invOddResidualNorm2;
    const long double t=(gb+geometry.kappa*gd)/
                           (static_cast<long double>(pooledTraceCount)*R);

    // First-order influence on log-amplitude from ALL pooled log observations:
    // local trace: w_i - gc*(cosh_i-cbar)/Vc - gd*sinh_i/Ss + t*r_i
    // other traces: t*r_i, where r_i=u_i-kappa*sinh_i.
    // The cross-trace term is essential: the pooled b is uncertain.
    long double variance=static_cast<long double>(pooledTraceCount-1)*t*t*R;
    long double rawRSS=0.;
    for(std::size_t i=start;i<end;++i) {
        const double u=step*(static_cast<double>(i)-static_cast<double>(center));
        const double c=std::cosh(u),s=std::sinh(u);
        const double z=sharedSlope*u-peak.curvature*c-peak.skew*s;
        const long double g=std::exp(z-maxLogShape);
        const long double observed=rawIntensity[i];
        const long double predicted=fittedA*g;
        const long double residual=observed-predicted;
        rawRSS+=residual*residual;
        const long double wi=observed*g/numerator;
        const long double ri=u-geometry.kappa*s;
        const long double influence=wi-gc*(c-geometry.meanCosh)*invVc
                                      -gd*s*invSs+t*ri;
        variance+=influence*influence;
    }
    if (!(variance>=0) || !std::isfinite(variance) ||
        !std::isfinite(static_cast<double>(variance)) ||
        !std::isfinite(static_cast<double>(rawRSS)))
        return result;
    result.valid=true;
    result.logAmplitude=static_cast<double>(fittedLogA);
    result.logCorrection=result.logAmplitude-peak.intercept;
    result.rawResidualSumSquares=static_cast<double>(rawRSS);
    result.logAmplitudeVarianceFactor=static_cast<double>(variance);
    return result;
}

// Convenience call for a retained multi-trace candidate. One refit per
// trace, no modification of the shape coefficients or original FitView.
inline std::vector<AmplitudeRefit> refitSelectedAmplitudes(
    const std::vector<const float*> &rawTraces, std::size_t length,
    const FitView &fit, const GeometryCache &cache)
{
    if (rawTraces.empty() || rawTraces.size()!=fit.locals.size())
        throw std::invalid_argument("Raw traces do not match pooled log fit");
    const auto &g=cache.at(fit.scale);
    std::vector<AmplitudeRefit> corrected;
    corrected.reserve(rawTraces.size());
    for(std::size_t j=0;j<rawTraces.size();++j)
        corrected.push_back(refitRawAmplitude(rawTraces[j],length,fit.center,
                         g,cache.step(),fit.locals[j],fit.sharedSlope,rawTraces.size()));
    return corrected;
}

} // namespace qAlgorithms::detail::hyperbolic
