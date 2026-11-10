"""Compare separately linked immutable-baseline and Stage-A probes exactly."""
import json
import statistics
import subprocess
import sys
from pathlib import Path
old, new = sys.argv[1:3]
a = subprocess.run([old], check=True, capture_output=True).stdout
b = subprocess.run([new], check=True, capture_output=True).stdout
assert a == b, 'Legacy differential mismatch (zero tolerance)'
measurements = [[], []]
for trial in range(5):
    for index in ([0, 1] if trial % 2 == 0 else [1, 0]):
        p = subprocess.run([str(Path([old, new][index]).with_name('test_stage_a_benchmark' + Path([old, new][index]).suffix))], check=True, capture_output=True)
        measurements[index].extend(float(line.split()[0]) for line in p.stderr.decode().splitlines())
summary = {'exact_bytes': len(a), 'trials': 55, 'ns_per_window': []}
for samples in measurements:
    summary['ns_per_window'].append({'median': statistics.median(samples), 'min': min(samples), 'max': max(samples), 'stdev': statistics.stdev(samples)})
summary['ratio'] = summary['ns_per_window'][1]['median']/summary['ns_per_window'][0]['median']
print(json.dumps(summary, indent=2))
