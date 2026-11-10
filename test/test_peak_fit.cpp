#include "qalgorithms_peak_fit.h"
// Reuse independent dense OLS, quadrature and finite-difference references,
// rather than copying their implementation or using the kernel as its oracle.
#define main researchReferenceTests
#include "../experiments/hyperbolic_running/test.cpp"
#undef main
using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

int main() {
 try {
    if(researchReferenceTests()!=0) return 1;
    for(double step:{.06,.12,.25}) for(size_t m:{1u,3u,10u}) {
        constexpr size_t n=77, center=38;
        std::vector<std::vector<float>> storage(m,std::vector<float>(n)),raw=storage;
        std::vector<const float*> logs;
        for(size_t j=0;j<m;++j) {
            for(size_t i=0;i<n;++i) {
                const double u=step*(double(i)-center);
                storage[j][i]=float(7+.02*j+.6*u-(2+.03*j)*std::cosh(u)-.2*std::sinh(u)+.001*std::sin(i*.7));
                raw[j][i]=std::exp(storage[j][i]);
            }
            logs.push_back(storage[j].data());
        }
        HyperbolicModel::Geometry cache(24,step);
        size_t previousCenter=0,previousScale=0,count=0;
        const void *borrowed=nullptr;
        PeakFitRequest request{PeakModelKind::Hyperbolic,logs,n,24,step};
        forEachPeakFit(request,[&](const PeakFitView &view) {
            const auto &fit=*std::get<const HyperbolicModel::Fit*>(view);
            require(fit.center>previousCenter || (fit.center==previousCenter && fit.scale>previousScale),"window order");
            previousCenter=fit.center; previousScale=fit.scale; ++count;
            if(borrowed) require(borrowed==fit.locals.data(),"borrowed storage reused");
            borrowed=fit.locals.data();
            if(fit.center!=center || (fit.scale!=8 && fit.scale!=16 && fit.scale!=24)) return;
            const auto ref=solveReference(logs,center,fit.scale,step);
            near(fit.sharedSlope,ref.beta[3*m],1e-7,"binding shared slope");
            for(size_t j=0;j<m;++j) {
                near(fit.locals[j].intercept,ref.beta[3*j],1e-7,"binding a");
                near(fit.locals[j].curvature,ref.beta[3*j+1],1e-7,"binding C");
                near(fit.locals[j].skew,ref.beta[3*j+2],1e-7,"binding D");
            }
            for(size_t i=0;i<3*m+1;++i) for(size_t j=0;j<3*m+1;++j)
                near(hyperbolicCovarianceFactor(cache.at(fit.scale),m,i,j),
                     ref.inverse[i?i-1:3*m][j?j-1:3*m],1e-6,"full binding covariance");
            const auto original=fit.locals[0];
            const auto p=evaluateHyperbolicWindow(fit,0,cache,raw[0].data(),n);
            require(p.finiteCoefficients && p.identifiable && p.integrable && p.amplitude.valid,"property validity");
            require(p.widthU.status==PropertyStatus::Unavailable && p.areaUncertainty.status==PropertyStatus::Unavailable && p.statisticalAcceptance.status==PropertyStatus::Unavailable,"no invented capabilities");
            require(fit.locals[0].intercept==original.intercept,"original log fit immutable");
            const auto before=evaluateHyperbolicWindow(fit,0,cache,nullptr,n);
            near(p.areaU.value/before.areaU.value,std::exp(p.amplitude.logCorrection),1e-8,"postfit area scaling");
            near(p.apexU.value,before.apexU.value,0,"postfit apex invariant");
            // Independent derivative root, not the analytic apex expression.
            double lo=-20,hi=20;
            for(int k=0;k<100;++k) {double mid=(lo+hi)/2; double d=fit.sharedSlope-original.curvature*std::sinh(mid)-original.skew*std::cosh(mid); if(d>0)lo=mid;else hi=mid;}
            near(p.apexU.value,(lo+hi)/2,1e-10,"apex vs derivative root");
            // Simpson integral of corrected shape in u.
            double integral=0; constexpr int intervals=20000; const double du=16./intervals;
            for(int k=0;k<=intervals;++k) {double u=-8+k*du; double value=std::exp(p.correctedIntercept+fit.sharedSlope*u-original.curvature*std::cosh(u)-original.skew*std::sinh(u)); integral+=(k==0||k==intervals?1:(k%2?4:2))*value;}
            near(p.areaU.value,integral*du/3,1e-7,"bound area vs quadrature");
        });
        require(count>0,"binding emits fits");
    }
    std::vector<float> logs(17,1),raw(17,0);
    PeakFitRequest request{PeakModelKind::Hyperbolic,{logs.data()},17,6,.1};
    int calls=0;
    auto reject=[&] {bool caught=false;try {forEachPeakFit(request,[&](const auto&){++calls;});}catch(const std::invalid_argument&){caught=true;} require(caught && calls==0,"reject before callbacks");};
    logs.back()=std::numeric_limits<float>::quiet_NaN(); reject(); logs.back()=1;
    request.model=static_cast<PeakModelKind>(99);reject(); request.model=PeakModelKind::Hyperbolic;
    request.fixedStep=0;reject();request.fixedStep=.1;
    HyperbolicModel::Geometry cache(6,.1);
    HyperbolicModel::Fit f;f.center=8;f.scale=6;f.sharedSlope=.5;f.identifiable=true;f.locals={{7,2,.2}};
    auto p=evaluateHyperbolicWindow(f,0,cache,raw.data(),17);
    require(!p.amplitude.valid && p.areaU.status==PropertyStatus::Unavailable,"invalid selected raw does not fall back");
    raw[8]=-1;p=evaluateHyperbolicWindow(f,0,cache,raw.data(),17);require(!p.amplitude.valid,"negative raw");
    raw[8]=std::numeric_limits<float>::infinity();p=evaluateHyperbolicWindow(f,0,cache,raw.data(),17);require(!p.amplitude.valid,"nonfinite raw");
    f.locals[0]={7,.2,.3};p=evaluateHyperbolicWindow(f,0,cache,nullptr,17);require(!p.integrable,"invalid tails");
    f.locals[0]={7,2,.2};f.sharedSlope=40;p=evaluateHyperbolicWindow(f,0,cache,nullptr,17);require(p.integrable && !p.apexInObservedSupport && p.apexU.status==PropertyStatus::Available,"support separate from integrability");
    f.sharedSlope=0;f.locals[0]={7,2,2-1e-8};p=evaluateHyperbolicWindow(f,0,cache,nullptr,17);require(p.integrable,"near tail boundary");
    f.locals[0]={7,1e6,0};p=evaluateHyperbolicWindow(f,0,cache,nullptr,17);require(p.areaU.status==PropertyStatus::NumericalFailure,"Bessel failure explicit");
    HyperbolicModel::Geometry tiny(6,.001);f.scale=2;f.locals[0]={7,2,.2};p=evaluateHyperbolicWindow(f,0,tiny,nullptr,17);require(!p.identifiable,"unidentifiable geometry not accepted");
    request={PeakModelKind::PiecewiseQuadratic,{logs.data()},17,6,1};
    std::vector<SharedShapeRegCoeffs> expected;forEachCoefficientsMultiTrace(request.logTraces,17,[&](const auto &fit){if(fit.scale<=6)expected.push_back(fit);});
    size_t i=0;forEachPeakFit(request,[&](const auto &view){const auto &fit=*std::get<const SharedShapeRegCoeffs*>(view);const auto &e=expected.at(i++);require(fit.b0==e.b0 && fit.b1==e.b1 && fit.b2==e.b2 && fit.b3==e.b3 && fit.x0==e.x0 && fit.scale==e.scale,"legacy diagnostic equivalence");});require(i==expected.size(),"legacy diagnostic count");
    std::cout<<"Stage B independent binding checks passed\n";
 } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
