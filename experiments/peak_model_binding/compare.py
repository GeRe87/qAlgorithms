"""Independent Stage-A/Stage-B binary differential and alternating paired timing."""
import json
import statistics
import subprocess
import sys
from pathlib import Path
old_probe,new_probe,old_bench,new_bench,extra_bench=sys.argv[1:6]
a=subprocess.run([old_probe],check=True,capture_output=True).stdout
b=subprocess.run([new_probe],check=True,capture_output=True).stdout
assert a==b,'Legacy output mismatch (zero tolerance)'
data=[{},{}]
counts=[{},{}]
for trial in range(9):
    for index in ([0,1] if trial%2==0 else [1,0]):
        output=subprocess.check_output([[old_bench,new_bench][index]],text=True)
        for line in output.splitlines():
            name,m,ns,count=line.split(',')
            key=name+'_'+m
            data[index].setdefault(key,[]).append(float(ns))
            counts[index][key]=int(count)
    # Binding-only costs come from a separate executable, so common timing
    # binaries have identical translation units and compiler settings.
    for line in subprocess.check_output([extra_bench],text=True).splitlines():
        name,m,ns,count=line.split(',')
        if name.startswith('binding_'):
            key=name+'_'+m
            data[1].setdefault(key,[]).append(float(ns))
            counts[1][key]=int(count)
for key in data[0]:
    assert counts[0][key]==counts[1][key],(key,'unequal windows/traces')
def summary(values):
    return dict(median=statistics.median(values),min=min(values),max=max(values),stdev=statistics.stdev(values))
result=dict(exact_legacy_bytes=len(a),paired_trials=9,baseline={k:summary(v) for k,v in data[0].items()},stage_b={k:summary(v) for k,v in data[1].items()},counts=counts[1])
result['ratios']={k:result['stage_b'][k]['median']/v['median'] for k,v in result['baseline'].items()}
print(json.dumps(result,indent=2))
