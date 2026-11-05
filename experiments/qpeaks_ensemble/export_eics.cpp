// Read-only diagnostic export of all EICs using the unchanged production path.
#include "qalgorithms_datatypes.h"
#include "qalgorithms_qbin.h"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_read_file.h"
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

template<class T> void writeValues(std::ofstream& out,const T* values,size_t n)
{
    out.write(reinterpret_cast<const char*>(values),std::streamsize(sizeof(T)*n));
    if(!out)throw std::runtime_error("EIC binary write failed.");
}

int main(int argc,char** argv)
{
    try {
        using namespace qAlgorithms;
        if(argc!=3)throw std::invalid_argument("Usage: export_eics INPUT.mzML OUTPUT_DIRECTORY");
        if(std::endian::native!=std::endian::little)throw std::runtime_error("This export requires little endian.");
        const std::filesystem::path path(argv[1]),folder(argv[2]);
        std::filesystem::create_directories(folder);
        XML_File file(path.c_str(),mzML);
        if(file.defective)throw std::runtime_error("Cannot read mzML.");
        std::ofstream binary(folder/"eics.bin",std::ios::binary);
        std::ofstream inventory(folder/"inventory.csv"),summary(folder/"polarities.csv");
        if(!binary||!inventory||!summary)throw std::runtime_error("Cannot open diagnostic output.");
        binary.write("QEIC0001",8);
        inventory<<std::setprecision(17)<<"polarity,binID,length,first_rt,last_rt,minimum_area,maximum_area,increasing_rt\n";
        summary<<"polarity,selected_ms1_profile_scans,centroids,eics,eic_observations\n";
        size_t records=0,selectedTotal=0;
        for(const auto polarity:{Polarities::positive,Polarities::negative}){
            const auto indices=filter_spectra(&file,true,polarity,true);
            const uint32_t code=polarity==Polarities::positive?1:2;
            std::cout<<"polarity="<<code<<" selected="<<indices.size()<<std::endl;
            selectedTotal+=indices.size();
            if(indices.empty()){
                summary<<code<<",0,0,0,0\n";
                continue;
            }
            std::vector<float> rt;
            get_spectra_RT(&file,&indices,&rt);
            if(rt.size()!=indices.size())throw std::runtime_error("RT selection mismatch.");
            std::vector<CentroidPeak> centroids;
            findCentroids(&file,&indices,&centroids);
            std::vector<EIC> bins;
            if(centroids.size()>=5)bins=performQbinning(&centroids,&rt);
            size_t observations=0;
            for(size_t id=0;id<bins.size();++id){
                const auto& e=bins[id];const size_t n=e.RT.size();
                if(n<1||n>UINT16_MAX||e.ints_area.size()!=n||e.mz.size()!=n||e.df.size()!=n||e.scanNumbers.size()!=n)
                    throw std::runtime_error("Inconsistent EIC vectors.");
                float lo=INFINITY,hi=0;bool increasing=true;
                for(size_t j=0;j<n;++j){
                    if(!std::isfinite(e.ints_area[j])||e.ints_area[j]<=0||!std::isfinite(e.RT[j])||!std::isfinite(e.mz[j])||e.df[j]!=j+1)
                        throw std::runtime_error("Invalid EIC observation; export cannot skip it.");
                    if(j&&e.RT[j]<=e.RT[j-1])increasing=false;
                    lo=std::min(lo,e.ints_area[j]);hi=std::max(hi,e.ints_area[j]);
                }
                const uint32_t header[]={code,uint32_t(id),uint32_t(n)};
                writeValues(binary,header,3);writeValues(binary,e.RT.data(),n);
                writeValues(binary,e.ints_area.data(),n);writeValues(binary,e.mz.data(),n);
                writeValues(binary,e.df.data(),n);writeValues(binary,e.scanNumbers.data(),n);
                inventory<<code<<','<<id<<','<<n<<','<<e.RT.front()<<','<<e.RT.back()<<','<<lo<<','<<hi<<','<<increasing<<'\n';
                observations+=n;++records;
            }
            summary<<code<<','<<indices.size()<<','<<centroids.size()<<','<<bins.size()<<','<<observations<<'\n';
            std::cout<<"polarity="<<code<<" centroids="<<centroids.size()<<" eics="<<bins.size()<<" observations="<<observations<<std::endl;
        }
        if(!selectedTotal)throw std::runtime_error("No supported MS1 profile scans selected.");
        for(auto* out:{&binary,&inventory,&summary}){out->close();if(!*out)throw std::runtime_error("Output finalization failed.");}
        std::cout<<"Complete EIC records="<<records<<std::endl;
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
