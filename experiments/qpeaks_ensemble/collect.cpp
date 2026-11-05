// Exploratory collector: no production acceptance, grouping or output changes.
// The source include exposes the exact existing moving-regression kernels.
#include "../../src/qalgorithms_qpeaks.cpp"
#include "centered_kernel.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <map>
#include <set>
#include <tuple>
using namespace qAlgorithms;

namespace {
struct Input {
    uint32_t scenario=0,trial=0,n=0,t=0;
    std::vector<float> truth,x,y;
    std::vector<uint16_t> df;
};
struct Node {
    size_t id,block,center,scale;
    RegCoeffs c;
    double apex=0,rss=0;
    bool maximum=false,constant=false,aboveFloor=false,currentAccepted=false;
    size_t root=0;
};
struct DSU {
    std::vector<size_t> p;
    explicit DSU(size_t n):p(n){std::iota(p.begin(),p.end(),0);}
    size_t find(size_t a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;}
    void join(size_t a,size_t b){a=find(a);b=find(b);if(a!=b)p[std::max(a,b)]=std::min(a,b);}
};
struct Component {
    size_t count=0,minScale=SIZE_MAX,maxScale=0;
    double minApex=INFINITY,maxApex=-INFINITY;
    std::set<size_t> scales;
};
template<class T> void read(std::istream& in,T* p,size_t n=1){
    in.read(reinterpret_cast<char*>(p),std::streamsize(sizeof(T)*n));
    if(!in)throw std::runtime_error("Truncated input record.");
}
void quote(std::ostream& out,const std::string& s){out<<'"';for(char c:s){if(c=='"')out<<'"';out<<c;}out<<'"';}
std::pair<bool,double> interiorMaximum(const RegCoeffs& c){
    double root=0;
    if(c.b1<0 && c.b2<0){root=-c.b1/(2*c.b2);return {std::isfinite(root)&&root>-double(c.scale)&&root<0,root};}
    if(c.b1>0 && c.b3<0){root=-c.b1/(2*c.b3);return {std::isfinite(root)&&root>0&&root<double(c.scale),root};}
    if(c.b1==0 && c.b2<0 && c.b3<0)return {true,0};
    return {false,0};
}
std::set<std::pair<size_t,size_t>> currentAccepted(const Input& in){
    std::set<std::pair<size_t,size_t>> result;
    for(size_t start=0;start<in.n;){
        while(start<in.n && !(in.y[start]>minIntensity_global))++start;
        size_t end=start;while(end<in.n && in.y[end]>minIntensity_global)++end;
        if(end-start>=MINLENGTH){
            const size_t n=end-start;
            std::vector<float> z(n);
            std::vector<uint16_t> df(n);
            const uint16_t previous=start?in.df[start-1]:0;
            for(size_t i=0;i<n;++i){z[i]=std::log(in.y[start+i]);df[i]=in.df[start+i]-previous;}
            std::vector<RegCoeffs> coefficients;findCoefficients(z.data(),n,&coefficients);
            std::vector<RegressionGauss> accepted;
            validateRegressions(in.y.data()+start,z.data(),in.x.data()+start,df.data(),&coefficients,n,&accepted);
            for(const auto& r:accepted)result.emplace(start+r.coeffs.x0,r.coeffs.scale);
        }
        if(end==in.n)break;
        start=end+1;
    }
    return result;
}
}

int main(int argc,char** argv){
 try{
    if(argc!=3 && argc!=4){std::cerr<<"Usage: collect INPUT.bin OUTPUT_DIRECTORY [--centered]\n";return 2;}
    const bool centered=argc==4;
    if(centered && std::string(argv[3])!="--centered")throw std::invalid_argument("Unknown collector mode.");
    std::ifstream source(argv[1],std::ios::binary);
    if(!source)throw std::runtime_error("Cannot open input.");
    char magic[8];read(source,magic,8);
    if(std::string(magic,8)!="QENS0001")throw std::runtime_error("Unknown input format.");
    const std::filesystem::path folder(argv[2]);std::filesystem::create_directories(folder);
    std::ofstream summary(folder/"traces.csv"),components(folder/"components.csv"),candidates(folder/"candidates.csv"),edges(folder/"edges.csv"),observations(folder/"observations.csv");
    if(!summary||!components||!candidates||!edges||!observations)throw std::runtime_error("Cannot open output files.");
    for(auto* out:{&summary,&components,&candidates,&edges,&observations})*out<<std::setprecision(17);
    summary<<"scenario,trial,n,truth_count,positive_blocks,short_positive_samples,nonpositive_samples,all_regressions,interior_maxima,constant_log_windows,current_accepted_windows,graph_components,graph_edges,strict_output_count,strict_error,detail\n";
    components<<"scenario,trial,component,nodes,scales,min_scale,max_scale,min_apex,max_apex\n";
    candidates<<"scenario,trial,id,block,center,scale,b0,b1,b2,b3,full_rss,has_interior_maximum,apex,constant_log_window,all_above_production_floor,current_accepted,component\n";
    edges<<"scenario,trial,from_id,to_id,axis\n";
    observations<<"scenario,trial,index,x,y,cumulative_df,log_y_float\n";
    size_t records=0;
    while(source.peek()!=std::char_traits<char>::eof()){
        Input in;read(source,&in.scenario);read(source,&in.trial);read(source,&in.n);read(source,&in.t);
        if(in.n>UINT16_MAX||in.n<1||in.t>in.n)throw std::runtime_error("Unsupported input sizes.");
        in.truth.resize(in.t);in.x.resize(in.n);in.y.resize(in.n);in.df.resize(in.n);
        read(source,in.truth.data(),in.t);read(source,in.x.data(),in.n);read(source,in.y.data(),in.n);read(source,in.df.data(),in.n);
        for(size_t i=0;i<in.n;++i){
            if(!std::isfinite(in.y[i])||in.y[i]<0||!std::isfinite(in.x[i])||(i&&!(in.x[i]>in.x[i-1]))||in.df[i]!=i+1)
                throw std::runtime_error("Input requires finite nonnegative y, increasing x, cumulative measured-point counts.");
        }
        const bool detail=in.trial==0;
        const auto accepted=currentAccepted(in);
        std::vector<Node> nodes;size_t blocks=0,shortSamples=0,nonpositive=0;
        for(size_t start=0;start<in.n;){
            while(start<in.n&&!(in.y[start]>0)){++start;++nonpositive;}
            size_t end=start;while(end<in.n&&in.y[end]>0)++end;
            if(end-start>=MINLENGTH){
                const size_t n=end-start;std::vector<float> z(n);
                for(size_t i=0;i<n;++i)z[i]=std::log(in.y[start+i]);
                std::vector<RegCoeffs> coeffs;
                if(centered)qpeaks_ensemble::findCenteredCoefficients(z.data(),n,&coeffs);
                else findCoefficients(z.data(),n,&coeffs);
                for(const auto& c:coeffs){
                    const size_t first=c.x0-c.scale,last=c.x0+c.scale;
                    double rss=0;bool constant=true,above=true;
                    for(size_t i=first;i<=last;++i){const double d=double(z[i])-regAt(&c,double(i)-c.x0);rss+=d*d;constant=constant&&z[i]==z[first];above=above&&in.y[start+i]>minIntensity_global;}
                    if(!std::isfinite(rss)||!std::isfinite(c.b0)||!std::isfinite(c.b1)||!std::isfinite(c.b2)||!std::isfinite(c.b3))throw std::runtime_error("Nonfinite candidate; collection cannot silently skip it.");
                    const auto maximum=interiorMaximum(c);
                    const size_t center=start+c.x0;
                    nodes.push_back({nodes.size(),blocks,center,c.scale,c,maximum.first?center+maximum.second:0,rss,maximum.first,constant,above,accepted.count({center,c.scale})!=0,0});
                }
                ++blocks;
            }else shortSamples+=end-start;
            if(end==in.n)break;
            start=end; // next iteration explicitly counts the zero separator
        }
        DSU dsu(nodes.size());std::map<std::pair<size_t,size_t>,size_t> grid;
        size_t maxima=0,constant=0,edgeCount=0;
        for(const auto& node:nodes){constant+=node.constant;if(node.maximum){++maxima;grid.emplace(std::make_pair(node.center,node.scale),node.id);}}
        for(const auto& node:nodes)if(node.maximum){
            for(int axis=0;axis<2;++axis){
                const auto key=axis==0?std::make_pair(node.center+1,node.scale):std::make_pair(node.center,node.scale+1);
                const auto it=grid.find(key);if(it==grid.end())continue;
                dsu.join(node.id,it->second);++edgeCount;
                if(detail)edges<<in.scenario<<','<<in.trial<<','<<node.id<<','<<it->second<<','<<(axis==0?"center":"scale")<<'\n';
            }
        }
        std::map<size_t,Component> groups;
        for(auto& node:nodes)if(node.maximum){
            node.root=dsu.find(node.id);auto& g=groups[node.root];++g.count;g.scales.insert(node.scale);g.minScale=std::min(g.minScale,node.scale);g.maxScale=std::max(g.maxScale,node.scale);g.minApex=std::min(g.minApex,node.apex);g.maxApex=std::max(g.maxApex,node.apex);
        }
        for(const auto& item:groups){const auto& g=item.second;components<<in.scenario<<','<<in.trial<<','<<item.first<<','<<g.count<<','<<g.scales.size()<<','<<g.minScale<<','<<g.maxScale<<','<<g.minApex<<','<<g.maxApex<<'\n';}
        int count=-1;std::string error;
        try{std::vector<RegressionGauss> peaks;const int status=qpeaks_find(in.y.data(),in.x.data(),in.df.data(),in.n,&peaks);if(status<0)error="negative processing status";else count=int(peaks.size());}
        catch(const std::exception& e){error=e.what();}
        summary<<in.scenario<<','<<in.trial<<','<<in.n<<','<<in.t<<','<<blocks<<','<<shortSamples<<','<<nonpositive<<','<<nodes.size()<<','<<maxima<<','<<constant<<','<<accepted.size()<<','<<groups.size()<<','<<edgeCount<<','<<count<<',';quote(summary,error);summary<<','<<detail<<'\n';
        if(detail){
            for(const auto& node:nodes){candidates<<in.scenario<<','<<in.trial<<','<<node.id<<','<<node.block<<','<<node.center<<','<<node.scale<<','<<node.c.b0<<','<<node.c.b1<<','<<node.c.b2<<','<<node.c.b3<<','<<node.rss<<','<<node.maximum<<',';if(node.maximum)candidates<<node.apex;candidates<<','<<node.constant<<','<<node.aboveFloor<<','<<node.currentAccepted<<',';if(node.maximum)candidates<<node.root;candidates<<'\n';}
            for(size_t i=0;i<in.n;++i){observations<<in.scenario<<','<<in.trial<<','<<i<<','<<in.x[i]<<','<<in.y[i]<<','<<in.df[i]<<',';if(in.y[i]>0)observations<<std::log(in.y[i]);observations<<'\n';}
        }
        ++records;
    }
    for(auto* out:{&summary,&components,&candidates,&edges,&observations}){out->close();if(!*out)throw std::runtime_error("Output write failed.");}
    std::cout<<"Complete records="<<records<<"; kernel="<<(centered?"centered-reference":"production-kernel")<<"; no ensemble peak decision was made.\n";
    return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
