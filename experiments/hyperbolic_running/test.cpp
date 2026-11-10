#include "hyperbolic_running.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace hyperbolic_research;

namespace {
using Matrix=std::vector<std::vector<long double>>;

void require(bool ok,const char *msg) {
    if (!ok) throw std::runtime_error(msg);
}
void near(double x,long double ref, double tol, const char *msg) {
    if (!std::isfinite(x) || std::abs(static_cast<long double>(x)-ref) >
        tol*(1.+std::abs(ref))) {
        std::cerr << msg << ": got " << x << ", reference " << static_cast<double>(ref) << '\n';
        throw std::runtime_error(msg);
    }
}

struct DenseResult {std::vector<long double> beta; Matrix inverse;};
DenseResult solveReference(const std::vector<const float*> &data, std::size_t center,
                           std::size_t h,double step) {
    const std::size_t m=data.size(), n=3*m+1;
    Matrix gram(n,std::vector<long double>(n,0));
    std::vector<long double> rhs(n,0);
    std::vector<long double> x(n,0);
    for (std::size_t j=0;j<m;++j)
        for (int index=-static_cast<int>(h);index<=static_cast<int>(h);++index) {
            std::fill(x.begin(),x.end(),0);
            const long double u=static_cast<long double>(step)*index;
            x[3*j]=1;
            x[3*j+1]=-std::cosh(u);
            x[3*j+2]=-std::sinh(u);
            x[n-1]=u;
            const long double y=data[j][center+index];
            for (std::size_t a=0;a<n;++a) {
                rhs[a]+=x[a]*y;
                for (std::size_t b=0;b<n;++b) gram[a][b]+=x[a]*x[b];
            }
        }
    Matrix inv(n,std::vector<long double>(n,0));
    for(std::size_t i=0;i<n;++i) inv[i][i]=1;
    for(std::size_t i=0;i<n;++i) {
        std::size_t pivot=i;
        for(std::size_t j=i+1;j<n;++j)
            if(std::abs(gram[j][i])>std::abs(gram[pivot][i]))pivot=j;
        require(std::abs(gram[pivot][i])>1e-23L,"Singular independent matrix");
        std::swap(gram[i],gram[pivot]);std::swap(rhs[i],rhs[pivot]);std::swap(inv[i],inv[pivot]);
        const long double div=gram[i][i];
        for(std::size_t k=0;k<n;++k){gram[i][k]/=div;inv[i][k]/=div;}
        rhs[i]/=div;
        for(std::size_t j=0;j<n;++j){
            if(j==i)continue;
            const long double factor=gram[j][i];
            for(std::size_t k=0;k<n;++k){gram[j][k]-=factor*gram[i][k];inv[j][k]-=factor*inv[i][k];}
            rhs[j]-=factor*rhs[i];
        }
    }
    return {rhs,inv};
}

void testDenseReference() {
    constexpr std::size_t length=117, center=58;
    const double step=0.12;
    GeometryCache cache(45,step);
    for (std::size_t m : {std::size_t(1),std::size_t(3),std::size_t(10)}) {
        std::vector<std::vector<float>> storage(m,std::vector<float>(length));
        std::vector<const float*> data;
        for (std::size_t j=0;j<m;++j){
            for(std::size_t i=0;i<length;++i){
                const double u=(static_cast<double>(i)-center)*step;
                const double y=9.0+.08*j+.66*u-(3.4+.04*j)*std::cosh(u)
                     -(-.2+.025*j)*std::sinh(u)+.006*std::sin(.31*(i+3)*(j+1));
                storage[j][i]=static_cast<float>(y);
            }
            data.push_back(storage[j].data());
        }
        std::vector<FitView> chosen;
        forEachRunning(data,length,cache,[&](const FitView &fit){
            if(fit.center==center && (fit.scale==8 || fit.scale==16 || fit.scale==30))
                chosen.push_back(fit);
        });
        require(chosen.size()==3,"Expected three selected windows");
        for(const auto &fit:chosen){
            const auto ref=solveReference(data,center,fit.scale,step);
            const auto &g=cache.at(fit.scale);
            near(fit.sharedSlope,ref.beta[3*m],1e-7,"Shared b vs full OLS");
            for(std::size_t j=0;j<m;++j){
                near(fit.locals[j].intercept,ref.beta[3*j],1e-7,"a vs full OLS");
                near(fit.locals[j].curvature,ref.beta[3*j+1],1e-7,"C vs full OLS");
                near(fit.locals[j].skew,ref.beta[3*j+2],1e-7,"D vs full OLS");
                near(g.varInterceptFactor(),ref.inverse[3*j][3*j],1e-6,"var a");
                near(g.varCurvatureFactor(),ref.inverse[3*j+1][3*j+1],1e-6,"var C");
                near(g.covInterceptCurvatureFactor(),ref.inverse[3*j][3*j+1],1e-6,"cov a,C");
                near(g.varLocalSkewFactor(m),ref.inverse[3*j+2][3*j+2],1e-6,"var D");
                near(g.covSharedSlopeLocalSkewFactor(m),ref.inverse[3*m][3*j+2],1e-6,"cov b,D");
                if(j>0)near(g.kappa*g.kappa*g.varSharedSlopeFactor(m),
                            ref.inverse[2][3*j+2],1e-6,"cov cross-D");
            }
            near(g.varSharedSlopeFactor(m),ref.inverse[3*m][3*m],1e-6,"var b");
            require(fit.identifiable,"Selected geometry should be identifiable");
        }
    }
}

void testExactSyntheticApexAndArea() {
    const double step=.1;
    GeometryCache cache(32,step);
    for(double b : {-1.2,0.0,1.2}) {
        const std::size_t mid=38;
        const double p=2.9,q=3.7;
        std::vector<float> samples(77);
        for(std::size_t i=0;i<samples.size();++i) {
            double u=step*(static_cast<double>(i)-mid);
            samples[i]=static_cast<float>(9.2+b*u-p*std::exp(u)-q*std::exp(-u));
        }
        bool found=false;
        forEachRunning({samples.data()},samples.size(),cache,[&](const FitView &fit){
            if(fit.center!=mid || fit.scale!=24)return;
            found=true;
            const auto &pk=fit.locals.front();
            const double theoreticalApex=std::asinh(b/(2*std::sqrt(p*q)))+0.5*std::log(q/p);
            near(apexU(pk,fit.sharedSlope),theoreticalApex,5e-5,"Apex from sampled profile");
            // Independent Simpson quadrature on the real u axis.
            const int n=60000;
            const double lo=-8.,hi=8.,h=(hi-lo)/n;
            double sum=0.;
            for(int k=0;k<=n;++k){
                double u=lo+k*h;
                double f=std::exp(9.2+b*u-p*std::exp(u)-q*std::exp(-u));
                sum+=(k==0||k==n?1:(k&1?4:2))*f;
            }
            const double referenceArea=sum*h/3.;
            near(std::exp(logAreaU(pk,fit.sharedSlope)),referenceArea,5e-5,"Bessel area vs quadrature");
        });
        require(found,"Missing apex and area test fit");
    }
}


void testClosedRawAmplitudeRefit() {
    constexpr std::size_t n=77, center=38, h=24;
    const double step=.1, b=.8, C=4.2, D=.35, correctA=9.3;
    GeometryCache cache(30,step);
    std::vector<float> raw(n);
    for(std::size_t i=0;i<n;++i){
        const double u=(static_cast<double>(i)-center)*step;
        raw[i]=static_cast<float>(std::exp(correctA+b*u-C*std::cosh(u)-D*std::sinh(u)));
    }
    const LocalFit original{correctA+.4,C,D};
    const auto ref=refitRawAmplitude(raw.data(),n,center,cache.at(h),step,original,b,1);
    require(ref.valid,"Exact amplitude refit should be valid");
    near(ref.logAmplitude,correctA,5e-7,"Exact amplitude intercept");
    near(ref.logCorrection,-.4,5e-7,"Exact amplitude correction");
    near(ref.rawResidualSumSquares,0,1e-7,"Perfect shape raw SSE");
    require(ref.logAmplitudeVarianceFactor>0,"Positive refit influence variance");
    const LocalFit corrected{ref.logAmplitude,C,D};
    near(std::exp(logAreaU(corrected,b)-logAreaU(original,b)),std::exp(-.4),
         1e-6,"Corrected full Bessel area scales with intercept only");

    std::vector<float> highCurvature(n);
    const double highC=80., highA=85.;
    for(std::size_t i=0;i<n;++i){
        const double u=(static_cast<double>(i)-center)*step;
        highCurvature[i]=static_cast<float>(std::exp(highA-highC*std::cosh(u)));
    }
    const LocalFit high{0.,highC,0.};
    const auto stable=refitRawAmplitude(highCurvature.data(),n,center,cache.at(h),
                                        step,high,0.,1);
    require(stable.valid,"Scaled exponential refit must avoid overflow and underflow");
    near(stable.logAmplitude,highA,5e-7,"Stable high-dynamic-range amplitude refit");

    std::vector<float> invalid(n,0.f);
    require(!refitRawAmplitude(invalid.data(),n,center,cache.at(h),step,original,b,1).valid,
            "Zero total intensity must not create a valid amplitude");
    invalid[center]=-1.f;
    require(!refitRawAmplitude(invalid.data(),n,center,cache.at(h),step,original,b,1).valid,
            "Negative raw intensities must be rejected");
    invalid[center]=std::numeric_limits<float>::quiet_NaN();
    require(!refitRawAmplitude(invalid.data(),n,center,cache.at(h),step,original,b,1).valid,
            "Nonfinite raw intensities must be rejected");
    bool thrown=false;
    try{(void)refitRawAmplitude(raw.data(),n,n-2,cache.at(h),step,original,b,1);}
    catch(const std::invalid_argument&){thrown=true;}
    require(thrown,"Out-of-range window must be rejected");
}

void testAmplitudeInfluenceAgainstFiniteDifferences() {
    const std::size_t m=3,n=25,center=12,h=9,target=1;
    const double step=.16, eps=.002;
    GeometryCache cache(11,step);
    std::vector<std::vector<float>> logs(m,std::vector<float>(n));
    for(std::size_t j=0;j<m;++j)
        for(std::size_t i=0;i<n;++i) {
            const double u=step*(static_cast<double>(i)-center);
            logs[j][i]=static_cast<float>(
                9.0+.13*j+.7*u-(4.4+.04*j)*std::cosh(u)
                -(.12+.04*j)*std::sinh(u)
                +.012*std::sin(.47*i+.27*j));
        }
    auto solveCorrected=[&](std::size_t changedTrace,std::size_t changedIndex,double delta){
        std::vector<std::vector<float>> values=logs;
        if(changedTrace<m) values[changedTrace][changedIndex]+=
                                     static_cast<float>(delta);
        std::vector<const float*> ptrs;
        for(auto &v:values)ptrs.push_back(v.data());
        const auto dense=solveReference(ptrs,center,h,step);
        std::vector<float> raw(n);
        for(std::size_t i=0;i<n;++i)raw[i]=std::exp(values[target][i]);
        const LocalFit fitted{static_cast<double>(dense.beta[3*target]),
                              static_cast<double>(dense.beta[3*target+1]),
                              static_cast<double>(dense.beta[3*target+2])};
        return refitRawAmplitude(raw.data(),n,center,cache.at(h),step,fitted,
                                  static_cast<double>(dense.beta[3*m]),m);
    };
    const auto analytical=solveCorrected(m,0,0.); // no change
    require(analytical.valid,"Reference influence calculation not valid");
    long double sumsq=0;
    for(std::size_t j=0;j<m;++j)
        for(std::size_t i=center-h;i<=center+h;++i){
            const auto positive=solveCorrected(j,i,eps);
            const auto negative=solveCorrected(j,i,-eps);
            require(positive.valid&&negative.valid,"Perturbed fit invalid");
            const long double derivative=(positive.logAmplitude-negative.logAmplitude)/
                                          (2*eps);
            sumsq+=derivative*derivative;
        }
    near(analytical.logAmplitudeVarianceFactor,sumsq,3e-3,
         "Full multitrace first-order intercept variance factor vs finite differences");
}

void testUnstableAndInvalidInputs(){
    bool caught=false;
    try {GeometryCache c(1,.1); (void)c;} catch(const std::invalid_argument&){caught=true;}
    require(caught,"Invalid geometry should fail");
    GeometryCache cache(12,.001);
    require(!cache.at(2).identifiable,"Small narrow window must flag ill-conditioning");
    LocalFit invalid;invalid.curvature=.1;invalid.skew=.3;
    require(!validIntegrablePeak(invalid),"Divergent tails must fail area check");
    require(std::isnan(apexU(invalid,.8)),"Invalid apex must return NaN");
}
}

int main(){
    try{
        testDenseReference();
        testExactSyntheticApexAndArea();
        testUnstableAndInvalidInputs();
        testClosedRawAmplitudeRefit();
        testAmplitudeInfluenceAgainstFiniteDifferences();
        std::cout<<"PASS: dense OLS coefficients and covariance; pooled traces 1/3/10;\n"
                 <<"  fixed-geometry running windows; closed Bessel area vs quadrature;\n"
                 <<"  apex and identifiability guards.\n";
    }catch(const std::exception &error){std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
}
