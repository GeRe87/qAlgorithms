"""Compare qPeaks kernels with independent null calibration and paired evaluation."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import random
import re
from statistics import mean, median
import subprocess

from benchmark import records, save, score_features, simulate, validate_tolerances


KERNELS = ("quadratic", "hyperbolic")
STRENGTHS = (.03, .1, .3, 1.)


def resolved_features(executable, source, directory, model):
    directory.mkdir()
    log_path = directory / "qalgorithms.log"
    with log_path.open("wb") as log:
        subprocess.run([executable, "-i", source, "-o", directory, "--peak-model", model, "--benchmark-mode"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    diagnostics = log_path.read_text(errors="replace")
    if "Comparison skipped EIC" in diagnostics:
        raise ValueError(f"Incomplete run. An entire chromatogram was skipped in {log_path}.")
    suffix = "_hyperbolic" if model == "hyperbolic" else ""
    path = directory / f"{source.stem}{suffix}_features.csv"
    return {
        "features": [{"id": row["ID"], "mz": float(row["mz"]), "rt_s": float(row["retentionTime"]),
                      "area": float(row["area"]), "quality": float(row["DQSF"])}
                     for row in records(path) if row["isUnresolved"] == "0"],
        "omitted_centroid_components": sum(map(int, re.findall(r"omitted (\d+) unresolved centroid", diagnostics))),
        "omitted_feature_components": sum(map(int, re.findall(r"recovered \d+ rows and omitted (\d+) components", diagnostics))),
    }


def fits_budget(null_runs, target):
    """Least restrictive threshold giving the requested calibration budget."""
    scores = sorted((r["quality"] for run in null_runs for r in run
                     if math.isfinite(r["quality"])), reverse=True)
    limit = math.floor(target * len(null_runs))
    return math.nextafter(scores[limit], math.inf) if len(scores) > limit else None


def accepted(rows, threshold):
    return [r for r in rows if math.isfinite(r["quality"]) and
            (threshold is None or r["quality"] >= threshold)]


def aggregate(cases, model, threshold, ppm, seconds):
    null_counts = [len(accepted(case["null"][model]["features"], threshold)) for case in cases]
    scores = []
    levels = {level: [0, 0] for level in STRENGTHS}
    families = {family: [0, 0] for family in ("gaussian", "sech")}
    for case in cases:
        truth = case["truth"]
        signal = case["signal"][model]
        result = score_features(truth, accepted(signal["features"], threshold), ppm, seconds)
        found = {pair["truth"] for pair in result["matches"]}
        for peak in truth:
            levels[peak["strength"]][0] += peak["id"] in found
            levels[peak["strength"]][1] += 1
            families[peak["family"]][0] += peak["id"] in found
            families[peak["family"]][1] += 1
        scores.append(result)
    pairs = [pair for result in scores for pair in result["matches"]]
    errors = {}
    for field in ("mz_error_ppm", "rt_error_s", "area_ratio"):
        values = [pair[field] for pair in pairs if pair[field] is not None]
        errors[field] = median(abs(v - (1 if field == "area_ratio" else 0)) for v in values) if values else None
    tp, fp, fn = (sum(result[key] for result in scores) for key in ("tp", "fp", "fn"))
    return {"held_out_null_counts": null_counts, "held_out_fp_per_null_run": mean(null_counts),
            "omitted_centroid_components": sum(case[phase][model]["omitted_centroid_components"]
                                               for phase in ("null", "signal") for case in cases),
            "omitted_feature_components": sum(case[phase][model]["omitted_feature_components"]
                                               for phase in ("null", "signal") for case in cases),
            "tp": tp, "fp": fp, "fn": fn, "recall": tp / (tp + fn),
            "by_strength": {str(level): {"tp": a, "peaks": b, "recall": a / b}
                            for level, (a, b) in levels.items()},
            "by_family": {family: {"tp": a, "peaks": b, "recall": a / b if b else None}
                          for family, (a, b) in families.items()},
            "median_absolute_error": errors, "per_seed": scores}


def compare(calibration, executable, output, *, train=12, evaluate=24, peaks=12,
            target=.5, ppm=5., seconds=5., first_seed=41000, jobs=4):
    if min(train, evaluate, jobs) < 1 or peaks < len(STRENGTHS):
        raise ValueError("Provide positive run and worker counts and at least four peaks per run.")
    if not math.isfinite(target) or target < 0:
        raise ValueError("The false alarm budget must be finite and nonnegative.")
    validate_tolerances(ppm, seconds)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    executable = Path(executable).resolve()

    def run(seed):
        directory = output / f"seed-{seed}"
        simulate(calibration, directory, seed=seed, peaks=peaks, replicates=2, dropout=1,
                 strengths=STRENGTHS, point_truth=False)
        truth = json.loads((directory / "synthetic-001-truth.json").read_text())
        null = directory / "synthetic-002.mzML"
        signal = directory / "synthetic-001.mzML"
        result = {"truth": truth, "null": {}, "signal": {}}
        for model in KERNELS:
            result["null"][model] = resolved_features(executable, null, directory / f"null-{model}", model)
            if seed >= first_seed + train:
                result["signal"][model] = resolved_features(executable, signal, directory / f"signal-{model}", model)
        print(f"Seed {seed} complete", flush=True)
        return result

    seeds = list(range(first_seed, first_seed + train + evaluate))
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        cases = list(pool.map(run, seeds))
    models = {}
    for model in KERNELS:
        null_runs = [case["null"][model]["features"] for case in cases[:train]]
        threshold = fits_budget(null_runs, target)
        models[model] = {"dqsf_threshold": threshold,
                         "calibration_null_counts_before_threshold": [len(x) for x in null_runs],
                         "calibration_fp_per_null_run": mean(len(accepted(x, threshold)) for x in null_runs),
                         **aggregate(cases[train:], model, threshold, ppm, seconds)}
    left = models["quadratic"]["per_seed"]
    right = models["hyperbolic"]["per_seed"]
    differences = [(a["tp"] - b["tp"]) / peaks for a, b in zip(left, right)]
    rng = random.Random(first_seed)
    bootstrap = sorted(mean(rng.choices(differences, k=evaluate)) for _ in range(2000))
    profile_path = Path(calibration)
    if profile_path.is_dir():
        profile_path /= "calibration.json"
    summary = {"calibration": str(profile_path.resolve()),
               "calibration_sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
               "executable": str(executable), "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
               "training_seeds": seeds[:train],
               "held_out_seeds": seeds[train:], "strengths": STRENGTHS, "peaks_per_signal_run": peaks,
               "target_fp_per_null_run": target, "matching": {"ppm": ppm, "seconds": seconds},
               "recall_difference_quadratic_minus_hyperbolic": mean(differences),
               "paired_seed_bootstrap_95pct": [bootstrap[50], bootstrap[-51]], "models": models}
    save(output / "comparison.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration")
    parser.add_argument("--qalgo", dest="executable", required=True)
    parser.add_argument("-o", "--output", required=True)
    for key, kind, default in (("train", int, 12), ("evaluate", int, 24), ("peaks", int, 12),
                               ("target", float, .5), ("ppm", float, 5.), ("seconds", float, 5.),
                               ("first-seed", int, 41000), ("jobs", int, 4)):
        parser.add_argument("--" + key, type=kind, default=default)
    args = vars(parser.parse_args())
    try:
        compare(**args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")
    print(Path(args["output"]).resolve() / "comparison.json")


if __name__ == "__main__":
    main()
