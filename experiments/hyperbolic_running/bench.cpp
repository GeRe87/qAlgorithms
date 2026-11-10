// Research-only paired benchmark: production qPeaks running multi-trace OLS
// versus cached hyperbolic running OLS. Same log spectra and same windows.
#include "hyperbolic_running.hpp"
#include "qalgorithms_qpeaks.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using hyperbolic_research::GeometryCache;
using hyperbolic_research::FitView;
using hyperbolic_research::LocalFit;

namespace {
volatile double measuredSink=0.;
constexpr double PI=3.1415926535897932384626433832795;
constexpr double STEP=0.08;
constexpr std::size_t LENGTH=129, CENTER=LENGTH/2, EVAL_SCALE=24;

struct Scenario {
    std::string shape;
    std::vector<std::vector<float>> noisy;
    std::vector<std::vector<float>> rawIntensity;
    std::vector<std::vector<double>> noiseless;
    std::vector<const float*> ptrs;
    std::vector<const float*> rawPtrs;
    std::vector<double> trueAreas;
    double trueApex=0;
};
Scenario makeScenario(std::string shape,std::size_t m) {
    Scenario s;
    s.shape=std::move(shape);
    s.noisy.resize(m,std::vector<float>(LENGTH));
    s.rawIntensity.resize(m,std::vector<float>(LENGTH));
    s.noiseless.resize(m,std::vector<double>(LENGTH));
    std::mt19937_64 rng(20261110);
    std::normal_distribution<double> noise(0.,.015);
    const double b=1.2, C=6., D=.4;
    const double p=(C+D)/2.,q=(C-D)/2.;
    const double mu=std::asinh(b/(2.*std::sqrt(p*q)))+.5*std::log(q/p);
    s.trueApex=(s.shape=="gaussian")?0.:mu;
    for(std::size_t j=0;j<m;++j) {
        const double height=900.+75.*static_cast<double>(j);
        const double logH=std::log(height);
        const double a=logH-b*mu+C*std::cosh(mu)+D*std::sinh(mu);
        if(s.shape=="gaussian")s.trueAreas.push_back(height*std::sqrt(2.*PI/6.));
        else s.trueAreas.push_back(2.*std::exp(a)*std::pow(q/p,b/2.)*
                                   std::cyl_bessel_k(b,2.*std::sqrt(p*q)));
        for(std::size_t i=0;i<LENGTH;++i) {
            const double u=STEP*(static_cast<double>(i)-static_cast<double>(CENTER));
            const double clean=(s.shape=="gaussian") ? logH-3.*u*u
                        : a+b*u-C*std::cosh(u)-D*std::sinh(u);
            s.noiseless[j][i]=clean;
            s.noisy[j][i]=static_cast<float>(clean+noise(rng));
            s.rawIntensity[j][i]=static_cast<float>(std::exp(s.noisy[j][i]));
        }
        s.ptrs.push_back(s.noisy[j].data());
        s.rawPtrs.push_back(s.rawIntensity[j].data());
    }
    return s;
}

struct Outcome {double rmse=0, rawRmse=0, apexMae=0, areaMape=0;std::size_t valid=0;};


double legacyRefitLogAmplitude(const float* raw,
                               const qAlgorithms::SharedShapeRegCoeffs &fit,
                               std::size_t center,std::size_t h) {
    // Reference implementation of the same scalar intensity-domain OLS
    // estimator as production correctB0(), limited to this selected window.
    double maxShape=-std::numeric_limits<double>::infinity();
    for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k) {
        const double z=fit.b1*k+(k<0?fit.b2:fit.b3)*k*k;
        maxShape=std::max(maxShape,z);
    }
    long double numerator=0,denominator=0;
    for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k) {
        const double z=fit.b1*k+(k<0?fit.b2:fit.b3)*k*k;
        const long double g=std::exp(z-maxShape);
        const double observed=raw[center+k];
        if(!(std::isfinite(observed)&&observed>=0))return
            std::numeric_limits<double>::quiet_NaN();
        numerator+=observed*g;
        denominator+=g*g;
    }
    if(!(numerator>0&&denominator>0))return
        std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(std::log(numerator/denominator)-maxShape);
}

double legacyArea(const qAlgorithms::SharedShapeRegCoeffs &fit,double intercept) {
    if(!(fit.b2<0 && fit.b3<0)) return std::numeric_limits<double>::quiet_NaN();
    const double l=-fit.b2, r=-fit.b3, b=fit.b1;
    const double L=std::sqrt(PI)/(2.*std::sqrt(l))*
        std::exp(b*b/(4.*l))*std::erfc(b/(2.*std::sqrt(l)));
    const double R=std::sqrt(PI)/(2.*std::sqrt(r))*
        std::exp(b*b/(4.*r))*std::erfc(-b/(2.*std::sqrt(r)));
    // Legacy coordinates are integer sample indices. STEP converts to u units.
    return std::exp(intercept)*(L+R)*STEP;
}

Outcome legacyOutcome(const Scenario &s,
                      const qAlgorithms::SharedShapeRegCoeffs &fit,bool rawRefit) {
    Outcome out;
    double squared=0.,squaredRaw=0.;
    for(std::size_t j=0;j<s.ptrs.size();++j) {
        const double intercept=rawRefit ?
            legacyRefitLogAmplitude(s.rawPtrs[j],fit,CENTER,EVAL_SCALE) : fit.b0[j];
        for(int k=-static_cast<int>(EVAL_SCALE);k<=static_cast<int>(EVAL_SCALE);++k){
            const double guess=intercept+fit.b1*k+(k<0?fit.b2:fit.b3)*k*k;
            const double truth=s.noiseless[j][CENTER+k];
            squared+=(guess-truth)*(guess-truth);
            const double d=std::exp(guess)-s.rawIntensity[j][CENTER+k];
            squaredRaw+=d*d;
        }
        const double curv=fit.b1>=0?fit.b3:fit.b2;
        const double apex=(curv<0)? -fit.b1/(2.*curv)*STEP :
                                      std::numeric_limits<double>::quiet_NaN();
        const double area=legacyArea(fit,intercept);
        if(std::isfinite(apex)&&std::isfinite(area)){
            out.apexMae+=std::abs(apex-s.trueApex);
            out.areaMape+=100.*std::abs(area/s.trueAreas[j]-1.);
            ++out.valid;
        }
    }
    out.rmse=std::sqrt(squared/(s.ptrs.size()*(2*EVAL_SCALE+1)));
    out.rawRmse=std::sqrt(squaredRaw/(s.ptrs.size()*(2*EVAL_SCALE+1)));
    if(out.valid){out.apexMae/=out.valid;out.areaMape/=out.valid;}
    return out;
}

Outcome hyperOutcome(const Scenario &s,const FitView &fit,
                     const GeometryCache &cache,bool rawRefit) {
    Outcome out;
    double squared=0.,squaredRaw=0.;
    const auto corrected=rawRefit ?
        hyperbolic_research::refitSelectedAmplitudes(s.rawPtrs,LENGTH,fit,cache) :
        std::vector<hyperbolic_research::AmplitudeRefit>{};
    for(std::size_t j=0;j<s.ptrs.size();++j) {
        LocalFit local=fit.locals[j];
        if(rawRefit) {
            if(!corrected[j].valid)continue;
            local.intercept=corrected[j].logAmplitude;
        }
        for(int k=-static_cast<int>(EVAL_SCALE);k<=static_cast<int>(EVAL_SCALE);++k){
            const double u=k*STEP;
            const double guess=hyperbolic_research::predict(local,fit.sharedSlope,u);
            const double truth=s.noiseless[j][CENTER+k];
            squared+=(guess-truth)*(guess-truth);
            const double d=std::exp(guess)-s.rawIntensity[j][CENTER+k];
            squaredRaw+=d*d;
        }
        const double apex=hyperbolic_research::apexU(local,fit.sharedSlope);
        const double logArea=hyperbolic_research::logAreaU(local,fit.sharedSlope);
        if(std::isfinite(apex)&&std::isfinite(logArea)){
            out.apexMae+=std::abs(apex-s.trueApex);
            out.areaMape+=100.*std::abs(std::exp(logArea)/s.trueAreas[j]-1.);
            ++out.valid;
        }
    }
    out.rmse=std::sqrt(squared/(s.ptrs.size()*(2*EVAL_SCALE+1)));
    out.rawRmse=std::sqrt(squaredRaw/(s.ptrs.size()*(2*EVAL_SCALE+1)));
    if(out.valid){out.apexMae/=out.valid;out.areaMape/=out.valid;}
    return out;
}

std::size_t windowCount(const GeometryCache &cache) {
    std::size_t n=0;
    for(std::size_t c=2;c+2<LENGTH;++c)
        for(std::size_t h=2;h<=std::min(cache.maxScale(),std::min(c,LENGTH-c-1));++h)
            ++n;
    return n;
}

template<class F> double seconds(F&& run) {
    const auto start=std::chrono::steady_clock::now();
    run();
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}

double median(std::vector<double> x){
    std::sort(x.begin(),x.end());
    return x[x.size()/2];
}

void runBenchmark(const std::string &shape,std::size_t m, int trials,int blocks) {
    auto data=makeScenario(shape,m);
    GeometryCache cache((LENGTH-1)/2,STEP);
    qAlgorithms::SharedShapeRegCoeffs legacySelected;
    FitView hyperSelected;
    std::size_t legacyWindows=0, hyperWindows=0;
    qAlgorithms::forEachCoefficientsMultiTrace(data.ptrs,LENGTH,
        [&](const auto &fit){
            ++legacyWindows;
            if(fit.x0==CENTER && fit.scale==EVAL_SCALE)legacySelected=fit;
        });
    hyperbolic_research::forEachRunning(data.ptrs,LENGTH,cache,
        [&](const auto &fit){
            ++hyperWindows;
            if(fit.center==CENTER && fit.scale==EVAL_SCALE)hyperSelected=fit;
        });
    const auto count=windowCount(cache);
    if(legacyWindows!=count || hyperWindows!=count)
        throw std::runtime_error("Window mismatch: benchmark is not paired");
    const auto oldOutcome=legacyOutcome(data,legacySelected,false);
    const auto newOutcome=hyperOutcome(data,hyperSelected,cache,false);
    const auto oldRefit=legacyOutcome(data,legacySelected,true);
    const auto newRefit=hyperOutcome(data,hyperSelected,cache,true);
    if(oldRefit.rawRmse>oldOutcome.rawRmse*(1.+1e-9) ||
       newRefit.rawRmse>newOutcome.rawRmse*(1.+1e-9))
        throw std::runtime_error("Raw height refit increased intensity RSS");
    auto legacyRun=[&](){
        double checksum=0.;
        for(int block=0;block<blocks;++block)
            qAlgorithms::forEachCoefficientsMultiTrace(data.ptrs,LENGTH,
                [&](const auto &fit){
                    checksum+=fit.b1+fit.b2+fit.b3;
                    for(double a:fit.b0)checksum+=a*.001;
                });
        measuredSink=checksum;
    };
    auto hyperRun=[&](){
        double checksum=0.;
        for(int block=0;block<blocks;++block)
            hyperbolic_research::forEachRunning(data.ptrs,LENGTH,cache,
                [&](const auto &fit){
                    checksum+=fit.sharedSlope;
                    for(const auto &p:fit.locals)
                        checksum+=p.intercept*.001+p.curvature*.0001+p.skew*.00001;
                });
        measuredSink=checksum;
    };
    legacyRun();hyperRun(); // warm both kernels, not timed
    std::vector<double> oldTimes,newTimes;
    for(int k=0;k<trials;++k){
        if((k&1)==0){oldTimes.push_back(seconds(legacyRun));newTimes.push_back(seconds(hyperRun));}
        else{newTimes.push_back(seconds(hyperRun));oldTimes.push_back(seconds(legacyRun));}
    }
    const double oldNs=median(oldTimes)*1e9/(blocks*count);
    const double newNs=median(newTimes)*1e9/(blocks*count);
    const std::size_t good=hyperSelected.identifiable ? static_cast<std::size_t>(
        std::count_if(hyperSelected.locals.begin(),hyperSelected.locals.end(),
            [](const LocalFit &p){return hyperbolic_research::validIntegrablePeak(p);})) : 0;
    std::cout<<shape<<','<<m<<','<<count<<','<<std::fixed<<std::setprecision(5)
             <<oldNs<<','<<newNs<<','<<newNs/oldNs<<','
             <<oldOutcome.rmse<<','<<newOutcome.rmse<<','
             <<oldOutcome.apexMae<<','<<newOutcome.apexMae<<','
             <<oldOutcome.areaMape<<','<<newOutcome.areaMape<<','
             <<oldRefit.areaMape<<','<<newRefit.areaMape<<','
             <<oldOutcome.rawRmse<<','<<newOutcome.rawRmse<<','
             <<oldRefit.rawRmse<<','<<newRefit.rawRmse<<','
             <<oldOutcome.valid<<','<<newOutcome.valid<<','<<good<<'\n';
}
} // namespace

int main(int argc,char**argv) {
    try {
        int trials=9,blocks=8;
        if(argc>1)trials=std::stoi(argv[1]);
        if(argc>2)blocks=std::stoi(argv[2]);
        if(trials<3||blocks<1||trials>1000||blocks>10000)
            throw std::invalid_argument("Usage: bench [trials>=3] [blocks>=1]");
        std::cout<<"shape,traces,windows_per_scan,legacy_ns_per_window,hyper_ns_per_window,"
                 <<"hyper_over_legacy_time,legacy_log_rmse_to_truth,hyper_log_rmse_to_truth,"
                 <<"legacy_apex_mae_u,hyper_apex_mae_u,legacy_area_mape_pct,"
                 <<"hyper_area_mape_pct,legacy_rawrefit_area_mape_pct,"
                 <<"hyper_rawrefit_area_mape_pct,legacy_raw_rmse,hyper_raw_rmse,"
                 <<"legacy_rawrefit_rmse,hyper_rawrefit_rmse,"
                 <<"legacy_valid,hyper_valid,hyper_integrable\n";
        for(const std::string shape: {"gaussian","hyperbolic"})
            for(std::size_t m : {std::size_t(1),std::size_t(4),std::size_t(16)})
                runBenchmark(shape,m,trials,blocks);
        if(!std::isfinite(measuredSink))throw std::runtime_error("Benchmark sink invalid");
        return 0;
    }catch(const std::exception &e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
