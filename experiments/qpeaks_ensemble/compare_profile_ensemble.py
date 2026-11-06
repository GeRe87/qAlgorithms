"""Compare the profile ensemble with archived paired outputs on identical fields."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from analyse import read_inputs
from bundle_method import topology
from compare_bundle_cascade import fields, match_positions
from profile_ensemble import describe


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(inputs, native, output, refit, representative):
    cases = {(r["scenario"], r["trial"]): r for r in json.loads((inputs/"cases.json").read_text())}
    source = read_inputs(inputs/"inputs.bin")
    assert cases.keys() == source.keys()
    field_path = native/"fields.bin"
    field_digest = sha256(field_path)
    input_digest = sha256(inputs/"inputs.bin")
    baselines = []
    for path in (refit, representative):
        report = json.loads(path.with_name("report.json").read_text())
        if (report["native_fields_sha256"], report["inputs_sha256"], report["topology_version"]) != (
                field_digest, input_digest, "model_singleton_v3"):
            raise ValueError(f"Baseline does not match the input, native fields and v3 topology: {path}")
        baselines.append({(r["scenario"], r["trial"]): r for r in json.loads(path.read_text())})
    output.mkdir(parents=True, exist_ok=True)
    results = []
    seen = set()
    for key, z, candidates, _ in fields(field_path):
        seen.add(key)
        started = time.perf_counter()
        basins, _, groups, membership = topology(candidates, model_boundaries=True,
                                                 log_values=z, boundary_scope="singleton_blocks")
        assert len(membership) == len(candidates)
        bundles = []
        for g in groups:
            models = [candidates[i] for b in g["basins"] for i in basins[b]["candidate_ids"]]
            try:
                status, peaks = describe(models)
            except (ArithmeticError, ValueError) as error:
                status, peaks = "explicit_error", []
                bundles.append(dict(bundle=g["id"], status=status, fits=len(models), error=str(error)))
                continue
            bundles.append(dict(bundle=g["id"], status=status, fits=len(models), peaks=peaks,
                                junction=g["contains_merge"] or g["contains_split"]))
        elapsed = time.perf_counter() - started
        peaks = [p for b in bundles for p in b.get("peaks", [])]
        old, single = [data[key] for data in baselines]
        assert len(groups) == old["atomic_bundles"] == single["atomic_bundles"]
        result = dict(cases[key], bundles=bundles, native_fits=len(candidates),
                      profile_apices=[p["apex"] for p in peaks], profile_seconds=elapsed,
                      refit_apices=old["ensemble_proposal_apices"],
                      representative_apices=single["ensemble_proposal_apices"],
                      strict_status=old["strict_status"], strict_count=old["strict_descriptions"],
                      strict_positions_axis=old["strict_positions_axis"],
                      collective_status=old["collective_status"], collective_count=old["collective_descriptions"],
                      collective_positions_axis=old["collective_positions_axis"])
        if cases[key]["reference_kind"] != "unknown_real":
            result["localization"] = match_positions(source[key][0], result["profile_apices"])
        results.append(result)
        print(f"{len(results)} {key} {len(peaks)} {elapsed:.3f}", flush=True)
    assert seen == cases.keys() == baselines[0].keys() == baselines[1].keys()
    with (output/"cases.json").open("w", encoding="utf-8") as stream:
        json.dump(results, stream, allow_nan=False)
    report = dict(cases=len(results), fits=sum(r["native_fits"] for r in results),
                  bundle_status=Counter(b["status"] for r in results for b in r["bundles"]),
                  shapes=Counter(p["shape_status"] for r in results for b in r["bundles"] for p in b.get("peaks", [])),
                  real_cases=[{k: r[k] for k in ("name", "strict_status", "strict_count",
                                                 "collective_status", "collective_count")}
                              | dict(profile=len(r["profile_apices"]), refit=len(r["refit_apices"]),
                                     representative=len(r["representative_apices"]))
                              for r in results if r["reference_kind"] == "unknown_real"],
                  stage_seconds=sum(r["profile_seconds"] for r in results),
                  fields_sha256=field_digest, inputs_sha256=input_digest,
                  statistical_acceptance=False, rt_area=False)
    (output/"report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "native", "output", "refit", "representative"):
        parser.add_argument(name, type=Path)
    run(**vars(parser.parse_args()))
