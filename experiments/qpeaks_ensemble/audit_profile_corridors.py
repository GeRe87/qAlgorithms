"""Attribute profile maxima to native corridors and quantify adjacent valleys."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from bundle_method import topology
from compare_bundle_cascade import fields
from compare_profile_ensemble import sha256
from profile_ensemble import arrays, evaluate, influence


def corridor_kind(apex, own, block_corridors):
    if any(left <= apex <= right for left, right in own):
        return "own_corridor"
    if any(left <= apex <= right for left, right in block_corridors):
        return "other_corridor"
    if min(left for left, _ in own) <= apex <= max(right for _, right in own):
        return "between_own_corridors"
    return "outside_corridors"


def union_corridors(corridors):
    union = []
    for left, right in sorted(corridors):
        if union and left <= union[-1][1]:
            union[-1] = (union[-1][0], max(right, union[-1][1]))
        else:
            union.append((left, right))
    return union


def valley(models, left, right, sigma=None):
    """Contrast adjacent maxima with their shared sampled separator."""
    if not left["right"] == right["left"]:
        raise ValueError("Adjacent maxima do not share a separator.")
    xs = [left["apex"], left["right"], right["apex"]]
    profile = arrays(models)
    levels, _, _ = evaluate(profile, xs)
    drops = levels[[0, 2]] - levels[1]
    result = dict(separator=xs[1], log_drops=drops.tolist())
    if sigma is not None:
        centers, scales = profile.T[:2]
        first = int(min(centers-scales))
        length = int(max(centers+scales))-first+1
        jacobians = [influence(profile, x, first, length)[0] for x in xs]
        contrasts = np.array(jacobians)[[0, 2]] - jacobians[1]
        result["conditional_log_drop_se"] = (sigma*np.linalg.norm(contrasts, axis=1)).tolist()
    return result


def run(native, pilot, output):
    report = json.loads((pilot/"report.json").read_text())
    field_path = native/"fields.bin"
    if sha256(field_path) != report["fields_sha256"]:
        raise ValueError("The profile results belong to a different native field bank.")
    cases = {(r["scenario"], r["trial"]): r for r in json.loads((pilot/"cases.json").read_text())}
    results = []
    for key, logs, models, _ in fields(field_path):
        case = cases.pop(key)
        basins, _, groups, _ = topology(models, model_boundaries=True, log_values=logs,
                                         boundary_scope="singleton_blocks")
        roots = {g["id"]: [basins[i]["root"] for i in g["basins"]] for g in groups}
        # Only inspect corridors on the same positive block.
        block = {g["id"]: basins[g["basins"][0]]["block"] for g in groups}
        grouped = defaultdict(list)
        for g in groups:
            grouped[block[g["id"]]].extend(roots[g["id"]])
        grouped = {b: union_corridors(c) for b, c in grouped.items()}
        saved = {b["bundle"]: b for b in case["bundles"]}
        assert set(saved) == set(roots)
        diagnostics = []
        for g in groups:
            ident = g["id"]
            old = saved[ident]
            peaks = old.get("peaks", [])
            assignments = [corridor_kind(p["apex"], roots[ident], grouped[block[ident]]) for p in peaks]
            valleys = []
            if len(peaks) > 1:
                members = [models[i] for b in g["basins"] for i in basins[b]["candidate_ids"]]
                valleys = [valley(members, a, b) for a, b in zip(peaks, peaks[1:])]
            diagnostics.append(dict(bundle=ident, junction=old.get("junction", False),
                                    profile_status=old["status"], apex=[p["apex"] for p in peaks],
                                    corridors=roots[ident], classifications=assignments,
                                    valleys=valleys))
        results.append(dict(scenario=key[0], trial=key[1], name=case["name"],
                            reference_kind=case["reference_kind"], bundles=diagnostics))
    if cases:
        raise ValueError("Some profile results have no native field record.")
    output.mkdir(parents=True, exist_ok=True)
    (output/"cases.json").write_text(json.dumps(results, allow_nan=False), encoding="utf-8")
    summary = {}
    strata = defaultdict(list)
    for case in results:
        strata["real" if case["reference_kind"] == "unknown_real" else case["name"]].append(case)
    for name, cases in strata.items():
        bundles = [b for c in cases for b in c["bundles"]]
        summary[name] = dict(
            cases=len(cases), bundles=len(bundles),
            apex_classifications=Counter(k for b in bundles for k in b["classifications"]),
            bundles_with_own_corridor_maxima=Counter(sum(k == "own_corridor" for k in b["classifications"])
                                                    for b in bundles),
            per_case_all_maxima=Counter(sum(len(b["apex"]) for b in c["bundles"]) for c in cases),
            per_case_own_corridor_maxima=Counter(sum(b["classifications"].count("own_corridor")
                                                     for b in c["bundles"]) for c in cases),
            within_bundle_valleys=sum(len(b["valleys"]) for b in bundles),
            valleys_below_both_apices=sum(min(v["log_drops"]) > 0 for b in bundles for v in b["valleys"]))
    (output/"report.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("native", "pilot", "output"):
        parser.add_argument(name, type=Path)
    run(**vars(parser.parse_args()))
