#pragma once
// Experimental arithmetic reference: same equidistant design and existing
// inverse matrices, with a symmetric log offset and double running sums.
namespace qpeaks_ensemble {
inline void findCenteredCoefficients(const float* z,size_t length,
                                    std::vector<qAlgorithms::RegCoeffs>* output)
{
    using namespace qAlgorithms;
    if(length<MINLENGTH)throw std::invalid_argument("Centered fit needs at least five points.");
    const size_t maxscale=std::min(maxscale_global,(length-1)/2);
    const size_t total=(maxscale-GLOBAL_MINSCALE+1)*length-maxscale*(maxscale+1)+2*GLOBAL_MINSCALE-2;
    output->resize(total);
    // An observed middle value (or the mean of the two middle values) makes
    // this arithmetic reference symmetric under reflection. It is an offset,
    // not a background estimate or a new statistical parameter.
    const double anchor=length%2?double(z[length/2]):(double(z[length/2-1])+double(z[length/2]))/2;
    std::vector<double> centered(length);
    for(size_t i=0;i<length;++i)centered[i]=double(z[i])-anchor;
    for(size_t center=GLOBAL_MINSCALE;center<length-GLOBAL_MINSCALE;++center){
        const double* c=centered.data()+center;
        double p0=c[-1]+c[0]+c[1],p1=c[1]-c[-1],p2=c[-1],p3=c[1];
        const size_t limit=std::min(maxscale,std::min(center,length-center-1));
        size_t offset=0;
        for(size_t scale=GLOBAL_MINSCALE;scale<=limit;++scale){
            const double left=c[-scale],right=c[scale],s=double(scale),sq=s*s;
            p0+=left+right;p1+=s*(right-left);p2+=sq*left;p3+=sq*right;
            const auto& m=qalgo_matInverse[scale];
            auto& result=output->at(center-scale+offset);
            result.b0=anchor+m.A*p0+m.B*(p2+p3);
            result.b1=m.C*p1+m.D*(p2-p3);
            result.b2=m.B*p0+m.D*p1+m.E*p2+m.F*p3;
            result.b3=m.B*p0-m.D*p1+m.F*p2+m.E*p3;
            result.scale=scale;result.x0=center;
            offset+=length-2*scale;
        }
    }
}
}
