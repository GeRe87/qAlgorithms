"""Package the preceding paired benchmark, reflections and optional local EICs.
Measurement values/provenance are written only to the requested local output.
"""
import argparse,csv,hashlib,json,math,struct
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paired_data',type=Path,help='Directory containing the preceding inputs.bin and scenarios.json')
    p.add_argument('output',type=Path);p.add_argument('--eic',type=Path,action='append',default=[])
    args=p.parse_args();blob=(args.paired_data/'inputs.bin').read_bytes();offset=0;records=[]
    meta=json.loads((args.paired_data/'scenarios.json').read_text(encoding='utf-8'))
    while offset<len(blob):
        scenario,trial,n,t=struct.unpack_from('<4I',blob,offset);offset+=16
        truth=struct.unpack_from('<'+'f'*t,blob,offset);offset+=4*t
        y=struct.unpack_from('<'+'f'*n,blob,offset);offset+=4*n
        records.append((scenario,trial,list(truth),list(map(float,range(n))),list(y)))
    if offset!=len(blob) or len(records)!=meta['records']:raise ValueError('Paired input record count mismatch.')
    original=list(meta['scenarios'])
    for s in original:
        r=next(r for r in records if r[0]==s['id'] and r[1]==0)
        scenario,trial,truth,x,y=r;reflected=sorted(len(y)-1-v for v in truth)
        records.append((1000+scenario,0,reflected,x,y[::-1]))
        meta['scenarios'].append(dict(id=1000+scenario,name=s['name']+'_reflected',samples=len(y),true_apices=reflected,trials=1,reflection_of=scenario))
    for j,path in enumerate(args.eic):
        with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        n=len(rows)
        if [int(r['df']) for r in rows]!=list(range(1,n+1)):raise ValueError('Only cumulative measured-point counts are supported.')
        records.append((2000+j,0,[],[float(r['x']) for r in rows],[float(r['y']) for r in rows]))
        meta['scenarios'].append(dict(id=2000+j,name='real_'+path.stem,samples=n,true_apices=[],trials=1,truth_unknown=True,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_file=path.name))
    for scenario,name,y,truth in [(3000,'constant_positive',[240.0]*99,[]),(3001,'exact_gaussian_float',[1000*math.exp(-(i-20)**2/18) for i in range(41)],[20.0])]:
        records.append((scenario,0,truth,list(map(float,range(len(y)))),y))
        meta['scenarios'].append(dict(id=scenario,name=name,samples=len(y),true_apices=truth,trials=1))
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'inputs.bin').open('wb') as f:
        f.write(b'QENS0001')
        for scenario,trial,truth,x,y in records:
            n=len(y);f.write(struct.pack('<4I',scenario,trial,n,len(truth)))
            for values in (truth,x,y):f.write(struct.pack('<'+'f'*len(values),*values))
            f.write(struct.pack('<'+'H'*n,*range(1,n+1)))
    meta.update(records=len(records),base_input_sha256=hashlib.sha256(blob).hexdigest(),input_sha256=hashlib.sha256((args.output/'inputs.bin').read_bytes()).hexdigest(),production_policy_changed=False,format='QENS0001; little endian uint32 scenario,trial,n,t; float32 truth[t],x[n],y[n]; uint16 cumulative_df[n]',graph='Axis-adjacent (center,scale) maximum-bearing models; diagnostic topology, not a same-peak relation.')
    if 'methods' in meta:meta['preceding_benchmark_methods']=meta.pop('methods')
    (args.output/'scenarios.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print('Prepared',len(records),'records.')

if __name__=='__main__':main()
