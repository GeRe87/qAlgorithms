"""Match the prior audit by original scan, center and half-width."""
import collections
import csv
import json
from pathlib import Path
import statistics
import sys

scans_path, windows_path, audit_path = map(Path, sys.argv[1:])
with scans_path.open(newline="") as stream:
    scans = list(csv.DictReader(stream))
with windows_path.open(newline="") as stream:
    fitted = {(int(r["scan"]), int(r["center"]), int(r["scale"])): r
              for r in csv.DictReader(stream)}
totals = {}
timing = {}
for name in scans[0]:
    if name.endswith("_us"):
        values = [float(row[name]) for row in scans]
        timing[name] = dict(median=statistics.median(values), min=min(values),
                            max=max(values), stdev=statistics.stdev(values))
    elif name != "scan":
        totals[name] = sum(int(row[name]) for row in scans)
comparison = collections.defaultdict(collections.Counter)
with audit_path.open() as stream:
    for line in stream:
        row = json.loads(line)
        if row["type"] != "window":
            continue
        key = f'{row["cohort"]}_n{2 * row["scale"] + 1}'
        counts = comparison[key]
        entry = fitted.get((row["scan"], row["center"], row["scale"]))
        counts["attempted"] += 1
        counts["old_positive_eligible"] += row["reason"] == "eligible"
        counts["segmented_evaluated"] += entry is not None
        if row["reason"] == "eligible" and entry is not None:
            counts["matched"] += 1
            counts["integrability_agrees"] += row["integrable"] == int(entry["integrable"])
            counts["availability_agrees"] += (row["inference_status"] == 1) == bool(int(entry["available"]))
        if row["reason"] == "eligible" and entry is None:
            counts["lost_with_positive_at_or_below_e"] += any(
                0 < x <= 2.7182818284590452 for x in row["raw"])
print(json.dumps(dict(scans=len(scans), totals=totals, timing_us_per_scan=timing,
                      matched_prior_audit=comparison), indent=2))
