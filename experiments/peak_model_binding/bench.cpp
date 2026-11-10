#include "qalgorithms_qpeaks.h"
#include "../hyperbolic_running/hyperbolic_running.hpp"
#ifdef STAGE_B
#include "qalgorithms_peak_fit.h"
#endif
#include <chrono>
#include <iostream>
#include <cmath>
using namespace qAlgorithms;
using Clock=std::chrono::steady_clock;
volatile double sink=0;
template<class F> void measure(const char *kind,size_t m,F run) {
    size_t count=0;double sum=0;auto start=Clock::now();
    for(int pass=0;pass<20;++pass)run(count,sum);
    double ns=std::chrono::duration<double,std::nano>(Clock::now()-start).count();sink=sum;
    std::cout<<kind<<','<<m<<','<<ns/count<<','<<count<<'\n';
}
int main() {
    constexpr size_t n=129;
    for(size_t m:{1u,4u,16u}) {
        std::vector<std::vector<float>> storage(m,std::vector<float>(n)),raw=storage;std::vector<const float*> logs;
        for(size_t j=0;j<m;++j){for(size_t i=0;i<n;++i){double u=.1*(double(i)-64);storage[j][i]=float(8+.03*j+.5*u-2*std::cosh(u)-.2*std::sinh(u));raw[j][i]=std::exp(storage[j][i]);}logs.push_back(storage[j].data());}
        if(m==1) measure("legacy_collect",m,[&](size_t &count,double &sum){std::vector<RegCoeffs> out;findCoefficients(logs[0],n,&out);count+=out.size();for(const auto &c:out)sum+=c.b0;});
        measure("legacy_stream",m,[&](size_t &count,double &sum){forEachCoefficientsMultiTrace(logs,n,[&](const auto &fit){++count;sum+=fit.b1;});});
        measure("geometry",m,[&](size_t &count,double &sum){hyperbolic_research::GeometryCache cache(64,.1);++count;sum+=cache.at(64).meanCosh;});
        hyperbolic_research::GeometryCache cache(64,.1);
        measure("research_hot",m,[&](size_t &count,double &sum){hyperbolic_research::forEachRunning(logs,n,cache,[&](const auto &fit){++count;sum+=fit.sharedSlope;});});
        hyperbolic_research::FitView selected;
        hyperbolic_research::forEachRunning(logs,n,cache,[&](const auto &fit){if(fit.center==64 && fit.scale==24)selected=fit;});
        measure("research_postfit",m,[&](size_t &count,double &sum){for(size_t j=0;j<m;++j){auto p=hyperbolic_research::refitRawAmplitude(raw[j].data(),n,64,cache.at(24),.1,selected.locals[j],selected.sharedSlope,m);++count;sum+=p.logAmplitude;}});
#ifdef STAGE_B
        measure("binding_stream",m,[&](size_t &count,double &sum){forEachPeakFit({PeakModelKind::Hyperbolic,logs,n,64,.1},[&](const PeakFitView &v){const auto &fit=*std::get<const peak_model::HyperbolicModel::Fit*>(v);++count;sum+=fit.sharedSlope;});});
        measure("binding_prepared",m,[&](size_t &count,double &sum){peak_model::HyperbolicModel::runningFit(logs,n,cache,[&](const auto &fit){++count;sum+=fit.sharedSlope;});});
        measure("binding_collect",m,[&](size_t &count,double &sum){std::vector<peak_model::HyperbolicModel::Fit> out;forEachPeakFit({PeakModelKind::Hyperbolic,logs,n,64,.1},[&](const PeakFitView &v){out.push_back(*std::get<const peak_model::HyperbolicModel::Fit*>(v));});count+=out.size();for(auto &fit:out)sum+=fit.sharedSlope;});
        measure("binding_properties",m,[&](size_t &count,double &sum){for(size_t j=0;j<m;++j){auto p=peak_model::evaluateHyperbolicWindow(selected,j,cache,raw[j].data(),n);++count;sum+=p.correctedIntercept;}});
#endif
    }
}
