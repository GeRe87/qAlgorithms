#include "qalgorithms_qpeaks.h"
#include <chrono>
#include <cmath>
#include <iostream>
using namespace qAlgorithms;
void benchmark() { // Identical separately compiled benchmark for both source versions.
  std::vector<float> y(512); for(size_t i=0;i<y.size();++i)y[i]=float(std::sin(i*.07));
  for(int trial=0;trial<11;++trial) { auto start=std::chrono::steady_clock::now(); size_t windows=0; double sum=0;
   for(int pass=0;pass<30;++pass) {std::vector<RegCoeffs> c; 
findCoefficients(y.data(),y.size(),&c);
windows+=c.size(); for(auto &v:c) sum+=v.b0;}
   auto end=std::chrono::steady_clock::now(); std::cerr<<std::chrono::duration<double,std::nano>(end-start).count()/windows<<' '<<sum<<'\n';
  }
 }

int main() { benchmark(); }
