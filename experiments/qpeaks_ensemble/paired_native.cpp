// Experimental paired benchmark: unchanged native kernel and both cascade policies.
#include "../../src/qalgorithms_qpeaks.cpp"
#include <bit>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
using namespace qAlgorithms;
template<class T> void get(std::istream& f,T* p,size_t n=1){f.read(reinterpret_cast<char*>(p),sizeof(T)*n);if(!f)throw std::runtime_error("Truncated paired input.");}
template<class T> void put(std::ostream& f,const T* p,size_t n=1){f.write(reinterpret_cast<const char*>(p),sizeof(T)*n);if(!f)throw std::runtime_error("Paired field write failed.");}
void quoted(std::ostream& f,const std::string& s){f<<'"';for(char c:s){if(c=='"')f<<'"';f<<c;}f<<'"';}
int main(int argc,char** argv){try{
    if((argc!=3&&argc!=4)||std::endian::native!=std::endian::little)throw std::runtime_error("Usage: paired_native INPUT_QENS OUTPUT_DIRECTORY [--cascade-domain] (little endian)");
    const bool cascadeDomain=argc==4;
    if(cascadeDomain&&std::string(argv[3])!="--cascade-domain")throw std::runtime_error("Unknown paired benchmark mode.");
    const double domainFloor=cascadeDomain?minIntensity_global:0.;
    std::ifstream in(argv[1],std::ios::binary);char magic[8];get(in,magic,8);
    if(std::string(magic,8)!="QENS0001")throw std::runtime_error("Wrong paired format.");
    std::filesystem::path folder(argv[2]);std::filesystem::create_directories(folder);
    std::ofstream policy(folder/"field-policy.json");
    policy<<std::setprecision(17)<<"{\"domain\":\""<<(cascadeDomain?"cascade_floor":"positive")<<"\",\"minimum_intensity\":"<<domainFloor<<",\"existing_cascade_floor\":"<<minIntensity_global<<"}\n";
    policy.close();if(!policy)throw std::runtime_error("Cannot write field policy.");
    std::ofstream binary(folder/"fields.bin",std::ios::binary),runs(folder/"cascade-runs.csv"),peaks(folder/"cascade-peaks.csv"),timing(folder/"field-timing.csv");
    if(!binary||!runs||!peaks||!timing)throw std::runtime_error("Cannot open paired outputs.");
    binary.write("QFLD0001",8);
    runs<<"scenario,trial,policy,status,description_count,collective_count,discarded_partial_count,seconds,error\n";
    peaks<<"scenario,trial,policy,output,position_axis,start,length,height,is_collective,collective_reason\n";
    timing<<"scenario,trial,samples,domain_blocks,short_domain_samples,domain_excluded_samples,regressions,kernel_seconds\n";
    for(auto* f:{&runs,&peaks,&timing})*f<<std::setprecision(17);
    size_t records=0;
    while(in.peek()!=std::char_traits<char>::eof()){
        uint32_t s,t,n,truthN;get(in,&s);get(in,&t);get(in,&n);get(in,&truthN);
        if(n<1||n>UINT16_MAX||truthN>n)throw std::runtime_error("Invalid paired sizes.");
        std::vector<float> truth(truthN),x(n),y(n);std::vector<uint16_t> df(n);
        get(in,truth.data(),truthN);get(in,x.data(),n);get(in,y.data(),n);get(in,df.data(),n);
        for(size_t i=0;i<n;++i)if(!std::isfinite(y[i])||y[i]<0||!std::isfinite(x[i])||(i&&x[i]<=x[i-1])||df[i]!=i+1)throw std::runtime_error("Invalid paired observation.");
        using Clock=std::chrono::steady_clock;
        for(bool collective:{false,true}){
            const auto began=Clock::now();std::vector<RegressionGauss> result;std::vector<UnresolvedPeak> groups;std::string error;int status=0;
            try{status=qpeaks_find(y.data(),x.data(),df.data(),n,&result,collective?&groups:nullptr);if(status<0)error="negative native processing status";}
            catch(const std::exception& e){status=-1;error=e.what();}
            const double elapsed=std::chrono::duration<double>(Clock::now()-began).count();
            const std::string policy=collective?"collective":"strict";
            runs<<s<<','<<t<<','<<policy<<','<<status<<','<<(error.empty()?int(result.size()):-1)<<','<<(error.empty()?int(groups.size()):-1)<<','<<(error.empty()?0:result.size())<<','<<elapsed<<',';quoted(runs,error);runs<<'\n';
            if(error.empty())for(size_t k=0;k<result.size();++k){
                const auto& p=result[k];const auto g=std::find_if(groups.begin(),groups.end(),[k](const UnresolvedPeak& v){return v.resultIndex==k;});
                peaks<<s<<','<<t<<','<<policy<<','<<k<<','<<p.position<<','<<p.span.startIdx<<','<<p.span.length<<','<<p.height<<','<<(g!=groups.end())<<','<<(g==groups.end()?0:int(g->reason))<<'\n';
            }
        }
        struct Block{uint32_t start,length;std::vector<RegCoeffs> coeffs;};
        std::vector<Block> blocks;size_t shortSamples=0,nonpositive=0,count=0;
        const auto began=Clock::now();std::vector<float> z(n,std::numeric_limits<float>::quiet_NaN());
        for(size_t i=0;i<n;++i)if(y[i]>0)z[i]=std::log(y[i]);
        for(size_t a=0;a<n;){
            while(a<n&&!(y[a]>domainFloor)){++a;++nonpositive;}
            size_t b=a;while(b<n&&y[b]>domainFloor)++b;
            if(b-a>=MINLENGTH){Block v{uint32_t(a),uint32_t(b-a),{}};findCoefficients(z.data()+a,b-a,&v.coeffs);count+=v.coeffs.size();blocks.push_back(std::move(v));}
            else shortSamples+=b-a;
            a=b;
        }
        const double elapsed=std::chrono::duration<double>(Clock::now()-began).count();
        uint32_t header[]={s,t,n,uint32_t(blocks.size())};put(binary,header,4);put(binary,z.data(),n);
        for(const auto& b:blocks){uint32_t h[]={b.start,b.length,uint32_t(b.coeffs.size())};put(binary,h,3);
            for(const auto& c:b.coeffs){uint16_t h2[]={uint16_t(b.start+c.x0),c.scale};const double vals[]={c.b0,c.b1,c.b2,c.b3};for(double v:vals)if(!std::isfinite(v))throw std::runtime_error("Nonfinite native coefficient.");put(binary,h2,2);put(binary,vals,4);}}
        timing<<s<<','<<t<<','<<n<<','<<blocks.size()<<','<<shortSamples<<','<<nonpositive<<','<<count<<','<<elapsed<<'\n';++records;
        std::cout<<"Complete paired case="<<records<<" scenario="<<s<<" trial="<<t<<std::endl;
    }
    for(auto* f:{&binary,&runs,&peaks,&timing}){f->close();if(!*f)throw std::runtime_error("Paired output finalization failed.");}
    std::cout<<"Complete paired records="<<records<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
