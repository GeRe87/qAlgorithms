// Real-data research comparator: same un-interpolated profile windows, two
// different pooled linear models, amplitude correction both ways.
// This is descriptive evidence; no certified mass/area reference exists.
#include "hyperbolic_running.hpp"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double STEP=0.35;         // fixed BEFORE inspecting any fit
constexpr double MIN_APEX=100.0;   // matches existing research mzML CLI
constexpr double MAX_SPACING_RATIO=1.45, MIN_SPACING_RATIO=.55;
constexpr double MAX_SHOULDER_FRACTION=.15;
constexpr int REPEATS=120;
volatile double sink=0.;

struct PeakWindow {
    size_t spectrum=0,anchorIndex=0;
    double mz=0,observedApex=0,spacingCV=0;
    std::vector<float> raw,log;
};
struct Batch {
    size_t spectrum=0,half=0,m=0;
    std::vector<PeakWindow> peaks;
    std::vector<const float*> rawPtrs, logPtrs;
    void bind() {
        rawPtrs.clear();logPtrs.clear();
        for(const auto &pk:peaks){
            rawPtrs.push_back(pk.raw.data());
            logPtrs.push_back(pk.log.data());
        }
    }
};
struct FitPair {
    qAlgorithms::SharedShapeRegCoeffs legacy;
    hyperbolic_research::FitView hyper;
    bool haveOld=false,haveNew=false;
};

double checkedNative(const float mz) {
    if(!(std::isfinite(mz)&&mz>0))throw std::runtime_error("Invalid m/z");
    return std::sqrt(static_cast<double>(mz));
}

std::vector<PeakWindow> candidates(const std::vector<float> &mz,
                                    const std::vector<float> &intensity,
                                    size_t spectrum,size_t half) {
    std::vector<PeakWindow> out;
    if(mz.size()!=intensity.size() || mz.size()<2*half+17)return out;
    const size_t n=mz.size();
    std::vector<double> native(n),diff;
    diff.reserve(n-1);
    for(size_t i=0;i<n;++i) {
        native[i]=checkedNative(mz[i]);
        if(!std::isfinite(intensity[i]))throw std::runtime_error("Invalid intensity");
        if(i>0) {
            const double step=native[i]-native[i-1];
            if(!(step>0))throw std::runtime_error("Non-increasing native coordinate");
            diff.push_back(step);
        }
    }
    // Conservative lower-third quantile, matching the existing research CLI
    // selection policy; used for gap filtering, NOT the fitting basis scale.
    const size_t qi=diff.size()/3;
    std::nth_element(diff.begin(),diff.begin()+static_cast<std::ptrdiff_t>(qi),diff.end());
    const double globalStep=diff[qi];
    auto isMaximum=[&](size_t i){
        return i>0&&i+1<n&&intensity[i]>intensity[i-1]&&
               intensity[i]>=intensity[i+1];
    };
    for(size_t k=8+half;k+8+half<n;++k) {
        if(!isMaximum(k)||intensity[k]<MIN_APEX)continue;
        const double top=intensity[k];
        bool reject=false;
        // Same model-blind isolation screen for every candidate.
        for(size_t j=k-8;j<=k+8;++j)
            if(j!=k&&isMaximum(j)&&intensity[j]>MAX_SHOULDER_FRACTION*top) {
                reject=true;break;
            }
        if(reject)continue;
        std::vector<double> gaps;
        gaps.reserve(2*half);
        for(int side : {-1,1}) {
            size_t prev=k;
            double previous=top;
            for(size_t d=1;d<=half;++d) {
                const size_t idx=side<0?k-d:k+d;
                const double raw=intensity[idx];
                const double delta=std::abs(native[idx]-native[prev]);
                if(!(raw>0)||raw>previous*1.12 ||
                   delta<MIN_SPACING_RATIO*globalStep ||
                   delta>MAX_SPACING_RATIO*globalStep) {reject=true;break;}
                gaps.push_back(delta);
                prev=idx;previous=raw;
            }
            if(reject)break;
        }
        if(reject)continue;
        const double mean=std::accumulate(gaps.begin(),gaps.end(),0.0)/gaps.size();
        double variance=0;
        for(double d:gaps)variance+=(d-mean)*(d-mean);
        const double cv=std::sqrt(variance/gaps.size())/mean;
        if(cv>.20)continue; // no interpolation to an artificially regular grid
        PeakWindow pk;
        pk.spectrum=spectrum;pk.anchorIndex=k;pk.mz=mz[k];
        pk.observedApex=top;pk.spacingCV=cv;
        for(size_t i=k-half;i<=k+half;++i) {
            pk.raw.push_back(intensity[i]);
            pk.log.push_back(static_cast<float>(std::log(intensity[i])));
        }
        out.push_back(std::move(pk));
    }
    return out;
}

std::vector<Batch> selectedBatches(const std::filesystem::path &path,
                                   std::size_t maxScans) {
    qAlgorithms::XML_File file(path.c_str(),qAlgorithms::mzML);
    if(file.defective)throw std::runtime_error("Source mzML could not be read");
    auto selected=qAlgorithms::filter_spectra(&file,true,
                  qAlgorithms::Polarities::positive,true);
    if(selected.empty())throw std::runtime_error("No positive profile MS1 spectra");
    if(selected.size()>maxScans) {
        // Predeclared, approximately uniform retention-time coverage, rather
        // than only the earliest spectra. Never choose scans by fitted outcome.
        std::vector<uint32_t> stratified;
        stratified.reserve(maxScans);
        for(size_t k=0;k<maxScans;++k)
            stratified.push_back(selected[
                ((2*k+1)*selected.size())/(2*maxScans)]);
        selected=std::move(stratified);
    }
    std::vector<Batch> batches;
    std::vector<float> mz,intensity;
    size_t extracted=0;
    for(uint32_t scan:selected){
        if(qAlgorithms::get_spectrum(&file,&mz,&intensity,scan)!=0)
            throw std::runtime_error("mzML binary decode failed on selected scan");
        for(size_t h:{size_t(2),size_t(3),size_t(4)}) {
            const auto found=candidates(mz,intensity,scan,h);
            extracted+=found.size();
            for(size_t m:{size_t(1),size_t(2),size_t(4)}) {
                // Non-overlapping, sorted-by-m/z groups, fixed before fitting.
                for(size_t pos=0;pos+m<=found.size();pos+=m) {
                    Batch b;
                    b.spectrum=scan;b.half=h;b.m=m;
                    b.peaks.insert(b.peaks.end(),found.begin()+pos,
                                              found.begin()+pos+m);
                    batches.push_back(std::move(b));
                }
            }
        }
    }
    for(auto &b:batches)b.bind();
    std::cerr<<"decoded_scans="<<selected.size()
             <<" screened_peak_windows="<<extracted
             <<" grouped_model_pairs="<<batches.size()<<'\n';
    if(batches.empty())throw std::runtime_error("No pre-specified eligible real peak windows");
    return batches;
}

FitPair fitModels(const Batch &batch,const hyperbolic_research::GeometryCache &cache){
    FitPair p;
    const size_t n=2*batch.half+1;
    qAlgorithms::forEachCoefficientsMultiTrace(batch.logPtrs,n,[&](const auto &f){
        if(f.x0==batch.half&&f.scale==batch.half){p.legacy=f;p.haveOld=true;}
    });
    hyperbolic_research::forEachRunning(batch.logPtrs,n,cache,[&](const auto &f){
        if(f.center==batch.half&&f.scale==batch.half){p.hyper=f;p.haveNew=true;}
    });
    if(!p.haveOld||!p.haveNew)throw std::runtime_error("Paired comparison window missing");
    return p;
}

double legacyPred(const qAlgorithms::SharedShapeRegCoeffs &f,double a,int k){
    return a+f.b1*k+(k<0?f.b2:f.b3)*k*k;
}
double hyperPred(const hyperbolic_research::FitView &f,double a,
                 size_t j,double u){
    const auto &p=f.locals[j];
    return a+f.sharedSlope*u-p.curvature*std::cosh(u)-p.skew*std::sinh(u);
}

double legacyAmplitude(const qAlgorithms::SharedShapeRegCoeffs &fit,
                       const PeakWindow &peak,std::size_t h){
    double high=-std::numeric_limits<double>::infinity();
    for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k)
        high=std::max(high,legacyPred(fit,0,k));
    long double num=0,den=0;
    for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k){
        const long double shape=std::exp(legacyPred(fit,0,k)-high);
        num+=peak.raw[h+k]*shape;den+=shape*shape;
    }
    if(!(num>0&&den>0))return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(std::log(num/den)-high);
}

double legacyArea(const qAlgorithms::SharedShapeRegCoeffs &fit,double a){
    if(!(fit.b2<0&&fit.b3<0))return std::numeric_limits<double>::quiet_NaN();
    const double l=-fit.b2,r=-fit.b3,b=fit.b1;
    const double L=std::sqrt(std::acos(-1.0))/(2*std::sqrt(l))*
           std::exp(b*b/(4*l))*std::erfc(b/(2*std::sqrt(l)));
    const double R=std::sqrt(std::acos(-1.0))/(2*std::sqrt(r))*
           std::exp(b*b/(4*r))*std::erfc(-b/(2*std::sqrt(r)));
    return std::exp(a)*(L+R)*STEP;
}

double aicc(double rss,size_t n,size_t k){
    if(!(rss>0)||n<=k+1)return std::numeric_limits<double>::quiet_NaN();
    return n*std::log(rss/n)+2.0*k+2.0*k*(k+1)/(n-k-1);
}
void writeRows(const std::vector<Batch> &batches,std::ostream &stream){
    stream<<"scan,half_width,trace_count,anchor_mz,anchor_intensity,spacing_cv,"
             "old_log_rss,new_log_rss,old_AICc,new_AICc,"
             "old_raw_rss_pre,new_raw_rss_pre,old_raw_rss_post,new_raw_rss_post,"
             "old_apex_offset_samples,new_apex_offset_samples,"
             "old_area_u_pre,new_area_u_pre,old_area_u_post,new_area_u_post,"
             "old_integrable,new_integrable,raw_refit_valid\n";
    stream<<std::setprecision(12);
    size_t accepted=0,postfit=0;
    hyperbolic_research::GeometryCache cache(4,STEP);
    for(const auto &batch:batches) {
        const FitPair fit=fitModels(batch,cache);
        const size_t h=batch.half,n=2*h+1,m=batch.m;
        double oldLog=0,newLog=0;
        for(size_t j=0;j<m;++j)
            for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k) {
                const double observed=batch.peaks[j].log[h+k];
                const double u=STEP*k;
                const double old=legacyPred(fit.legacy,fit.legacy.b0[j],k);
                const double neu=hyperPred(fit.hyper,fit.hyper.locals[j].intercept,j,u);
                oldLog+=(old-observed)*(old-observed);
                newLog+=(neu-observed)*(neu-observed);
            }
        const double oldAicc=aicc(oldLog,n*m,m+3);
        const double newAicc=aicc(newLog,n*m,3*m+1);
        const auto corrected=hyperbolic_research::refitSelectedAmplitudes(
                                                batch.rawPtrs,n,fit.hyper,cache);
        for(size_t j=0;j<m;++j){
            const auto &pk=batch.peaks[j];
            const bool oldValid=fit.legacy.b2<0&&fit.legacy.b3<0;
            const bool newValid=fit.hyper.identifiable&&
                hyperbolic_research::validIntegrablePeak(fit.hyper.locals[j]);
            const double oldA=fit.legacy.b0[j];
            const double newA=fit.hyper.locals[j].intercept;
            const double oldCorrected=legacyAmplitude(fit.legacy,pk,h);
            const bool correctionValid=corrected[j].valid&&std::isfinite(oldCorrected);
            const double newCorrected=corrected[j].logAmplitude;
            const double oldApex=oldValid ? -fit.legacy.b1/
                     (2*(fit.legacy.b1<0?fit.legacy.b2:fit.legacy.b3)) :
                     std::numeric_limits<double>::quiet_NaN();
            const double newApex=newValid ?
                     hyperbolic_research::apexU(fit.hyper.locals[j],fit.hyper.sharedSlope)/STEP :
                     std::numeric_limits<double>::quiet_NaN();
            const double oldAreaPre=oldValid ? legacyArea(fit.legacy,oldA):
                     std::numeric_limits<double>::quiet_NaN();
            const double oldAreaPost=oldValid&&correctionValid?
                     legacyArea(fit.legacy,oldCorrected):
                     std::numeric_limits<double>::quiet_NaN();
            const double newAreaPre=newValid?
                     std::exp(hyperbolic_research::logAreaU(fit.hyper.locals[j],fit.hyper.sharedSlope)):
                     std::numeric_limits<double>::quiet_NaN();
            auto correctedLocal=fit.hyper.locals[j];
            correctedLocal.intercept=newCorrected;
            const double newAreaPost=newValid&&correctionValid?
                     std::exp(hyperbolic_research::logAreaU(correctedLocal,fit.hyper.sharedSlope)):
                     std::numeric_limits<double>::quiet_NaN();
            double oldRawBefore=0,newRawBefore=0,oldRawAfter=0,newRawAfter=0;
            for(int k=-static_cast<int>(h);k<=static_cast<int>(h);++k) {
                const double observed=pk.raw[h+k];
                const double beforeOld=std::exp(legacyPred(fit.legacy,oldA,k));
                const double beforeNew=std::exp(hyperPred(fit.hyper,newA,j,k*STEP));
                const double afterOld=correctionValid?
                     std::exp(legacyPred(fit.legacy,oldCorrected,k)):beforeOld;
                const double afterNew=correctionValid?
                     std::exp(hyperPred(fit.hyper,newCorrected,j,k*STEP)):beforeNew;
                oldRawBefore+=(beforeOld-observed)*(beforeOld-observed);
                newRawBefore+=(beforeNew-observed)*(beforeNew-observed);
                oldRawAfter+=(afterOld-observed)*(afterOld-observed);
                newRawAfter+=(afterNew-observed)*(afterNew-observed);
            }
            if(correctionValid&&(oldRawAfter>oldRawBefore*(1+1e-6) ||
                               newRawAfter>newRawBefore*(1+1e-6)))
                throw std::runtime_error("Raw amplitude correction increased RSS");
            accepted+=static_cast<size_t>(oldValid&&newValid);
            postfit+=static_cast<size_t>(correctionValid);
            stream<<batch.spectrum<<','<<h<<','<<m<<','<<pk.mz<<','
                  <<pk.observedApex<<','<<pk.spacingCV<<','
                  <<oldLog<<','<<newLog<<','<<oldAicc<<','<<newAicc<<','
                  <<oldRawBefore<<','<<newRawBefore<<','
                  <<oldRawAfter<<','<<newRawAfter<<','
                  <<oldApex<<','<<newApex<<','
                  <<oldAreaPre<<','<<newAreaPre<<','
                  <<oldAreaPost<<','<<newAreaPost<<','
                  <<oldValid<<','<<newValid<<','<<correctionValid<<'\n';
        }
    }
    std::cerr<<"matched_integrable_traces="<<accepted
             <<" matched_raw_refits="<<postfit<<'\n';
}
void benchmark(const std::vector<Batch> &batches) {
    hyperbolic_research::GeometryCache cache(4,STEP);
    size_t oldFits=0,newFits=0;
    for(int which=0;which<2;++which) {
        double checksum=0;
        const auto start=std::chrono::steady_clock::now();
        for(int iter=0;iter<REPEATS;++iter)
            for(const auto &batch:batches){
                const size_t n=2*batch.half+1;
                if(which==0)
                    qAlgorithms::forEachCoefficientsMultiTrace(batch.logPtrs,n,[&](const auto& f){
                        checksum+=f.b1;
                        ++oldFits;
                    });
                else
                    hyperbolic_research::forEachRunning(batch.logPtrs,n,cache,[&](const auto& f){
                        checksum+=f.sharedSlope;
                        ++newFits;
                    });
            }
        const double elapsed=std::chrono::duration<double>(
             std::chrono::steady_clock::now()-start).count();
        sink=checksum;
        std::cerr<<(which==0?"old":"hyper")<<"_all_windows_ns="
                 <<std::setprecision(6)<<elapsed*1e9/
                       (which==0?oldFits:newFits)<<'\n';
    }
    if(oldFits!=newFits)throw std::runtime_error("Mismatched candidate windows in benchmark");
}

} // namespace

int main(int argc,char **argv){
    try{
        if(argc!=4)throw std::invalid_argument(
           "Usage: real_qtof INPUT_PROFILE.mzML OUTPUT.csv MAX_MS1_SCANS");
        const std::filesystem::path input(argv[1]),output(argv[2]);
        const auto maxScans=static_cast<size_t>(std::stoul(argv[3]));
        if(maxScans<1||maxScans>1000)throw std::invalid_argument("MAX_MS1_SCANS outside [1,1000]");
        const auto batches=selectedBatches(input,maxScans);
        std::ofstream result(output);
        if(!result)throw std::runtime_error("Unable to open output CSV");
        writeRows(batches,result);
        result.close();
        if(!result)throw std::runtime_error("Cannot finalize real-data CSV");
        benchmark(batches);
    }catch(const std::exception &error){
        std::cerr<<"REAL_QTOF_RESEARCH_FAILED: "<<error.what()<<'\n';return 1;
    }
    return 0;
}
