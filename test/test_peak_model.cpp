#include "qalgorithms_qpeaks.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <stdexcept>
using namespace qAlgorithms;
static void coeff(const RegCoeffs &c) { std::cout<<c.b0<<' '<<c.b1<<' '<<c.b2<<' '<<c.b3<<' '<<c.x0<<' '<<c.scale<<'\n'; }
int main(int argc, char **) {
 std::cout<<std::hexfloat;
 for(int seed=0;seed<40;++seed) {
  size_t n=17+seed*3; std::vector<float> x(n),y(n),logs(n); std::vector<uint16_t> df(n);
  for(size_t i=0;i<n;++i) { double u=(double(i)-n*.47)/(2+seed%7); x[i]=float(100+i*.03); y[i]=float(4+800*std::exp(-u*u*.5)*(1+.015*std::sin(i*3.7+seed))); if(seed%5==0) y[i]+=float(300*std::exp(-std::pow((double(i)-n*.7)/3,2))); logs[i]=std::log(y[i]); df[i]=uint16_t(i+1); }
  std::vector<RegCoeffs> c; findCoefficients(logs.data(),n,&c); for(auto &v:c) coeff(v);
  for(bool collective:{false,true}) {
   std::vector<RegressionGauss> r; std::vector<UnresolvedPeak> groups;
   auto count=qpeaks_find(y.data(),x.data(),seed%2?df.data():nullptr,n,&r,collective?&groups:nullptr);
#ifdef STAGE_A
   std::vector<RegressionGauss> explicitResult; std::vector<UnresolvedPeak> eg;
   if(count!=qpeaks_find(y.data(),x.data(),seed%2?df.data():nullptr,n,&explicitResult,PeakFitOptions{},collective?&eg:nullptr)||r.size()!=explicitResult.size()) return 2;
   // Exact field serialization below is also applied to explicit results.
   for(size_t i=0;i<r.size();++i) { auto a=peak_model::fromLegacy(r[i],peak_model::AxisUnit::MassToCharge); auto b=peak_model::fromLegacy(explicitResult[i],peak_model::AxisUnit::MassToCharge); if(r[i].coeffs.b0!=explicitResult[i].coeffs.b0||r[i].coeffs.b1!=explicitResult[i].coeffs.b1||r[i].coeffs.b2!=explicitResult[i].coeffs.b2||r[i].coeffs.b3!=explicitResult[i].coeffs.b3||r[i].coeffs.x0!=explicitResult[i].coeffs.x0||r[i].coeffs.scale!=explicitResult[i].coeffs.scale||r[i].span.startIdx!=explicitResult[i].span.startIdx||r[i].span.length!=explicitResult[i].span.length||r[i].jaccard!=explicitResult[i].jaccard||r[i].dqs!=explicitResult[i].dqs||r[i].df!=explicitResult[i].df||r[i].numCompetitors!=explicitResult[i].numCompetitors||r[i].apexGroup!=explicitResult[i].apexGroup||a.position!=b.position||a.area!=b.area||a.height!=b.height||a.width!=b.width||a.positionUncertainty!=b.positionUncertainty||a.areaUncertainty!=b.areaUncertainty||a.heightUncertainty!=b.heightUncertainty||a.valid!=b.valid) return 3; }
   if(groups.size()!=eg.size()) return 5; for(size_t i=0;i<groups.size();++i) if(groups[i].resultIndex!=eg[i].resultIndex||groups[i].reason!=eg[i].reason||groups[i].start!=eg[i].start||groups[i].end!=eg[i].end||groups[i].maximaIndices!=eg[i].maximaIndices) return 6;
   for(auto kind:{PeakModelKind::Hyperbolic,static_cast<PeakModelKind>(99)}) { auto before=r.size(); bool threw=false; try {qpeaks_find(nullptr,nullptr,nullptr,0,&r,PeakFitOptions{kind},&groups);} catch(const std::invalid_argument &) {threw=true;} if(!threw||r.size()!=before) return 4; }
#endif
   std::cout<<count<<' '<<groups.size()<<'\n';
   for(auto &v:r) { coeff(v.coeffs); std::cout<<v.span.startIdx<<' '<<v.span.length<<' '<<v.area<<' '<<v.position<<' '<<v.height<<' '<<v.area_unc<<' '<<v.position_unc<<' '<<v.height_unc<<' '<<v.fwhm<<' '<<v.jaccard<<' '<<v.dqs<<' '<<v.df<<' '<<int(v.numCompetitors)<<' '<<int(v.apexGroup)<<' '<<v.isValid<<'\n'; }
   for(auto &g:groups) { std::cout<<g.resultIndex<<' '<<int(g.reason)<<' '<<g.start<<' '<<g.end; for(auto i:g.maximaIndices) std::cout<<' '<<i; std::cout<<'\n'; }
  }
 }
 (void)argc;
}
