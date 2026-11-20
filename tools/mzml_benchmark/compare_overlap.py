"""Compare kernels on matched isolated and overlapping signal placements."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import random
from statistics import mean, median
import subprocess

from benchmark import save, score_features, simulate, validate_tolerances
from compare_kernels import KERNELS, resolved_features


def evaluate(truth, rows, ppm, seconds):
    score = score_features(truth, rows, ppm, seconds)
    truth_by_id = {row["id"]: row for row in truth}
    found_by_id = {row["id"]: row for row in rows}
    matched = {row["truth"]: row for row in score["matches"]}
    groups = []
    for group in score["groups"].values():
        if len(group["present"]) != 2:
            continue
        ids = group["present"]
        selected = [matched[i] for i in ids if i in matched]
        summary = {"mode": group["mode"], "peaks": ids, "found": len(selected)}
        if len(selected) == 1:
            match = selected[0]
            measured = found_by_id[match["detected"]]["area"]
            individual = truth_by_id[match["truth"]]["area"]
            combined = sum(truth_by_id[i]["area"] for i in ids)
            summary.update(area_ratio_to_match=measured / individual, area_ratio_to_pair=measured / combined,
                           area_closer_to_pair=abs(measured - combined) < abs(measured - individual))
        groups.append(summary)
    return {"score": score, "pairs": groups}


def metrics(results):
    counts = Counter(pair["found"] for row in results for pair in row["pairs"])
    modes = {mode: Counter(pair["found"] for row in results for pair in row["pairs"] if pair["mode"] == mode)
             for mode in ("mz", "rt")}
    single = [pair for row in results for pair in row["pairs"] if pair["found"] == 1]
    errors = {}
    for field in ("mz_error_ppm", "rt_error_s", "area_ratio"):
        values = [abs(m[field] - (1 if field == "area_ratio" else 0))
                  for row in results for m in row["score"]["matches"] if m[field] is not None]
        errors[field] = median(values) if values else None
    return {"pairs_by_number_found": dict(counts), "by_mode": modes,
            "tp": sum(row["score"]["tp"] for row in results),
            "fn": sum(row["score"]["fn"] for row in results),
            "single_detection_median_area_ratio_to_member": median(p["area_ratio_to_match"] for p in single) if single else None,
            "single_detection_median_area_ratio_to_pair": median(p["area_ratio_to_pair"] for p in single) if single else None,
            "single_detection_closer_to_pair": sum(p["area_closer_to_pair"] for p in single),
            "mean_fp_per_run": mean(row["score"]["fp"] for row in results),
            "median_absolute_error": errors}


def shared_match_errors(results, model):
    examples = []
    for row in results:
        isolated, overlap = ({p["truth"]: p for p in row[condition][model]["score"]["matches"]}
                             for condition in ("isolated", "overlap"))
        paired_ids = {i for pair in row["overlap"][model]["pairs"] for i in pair["peaks"]}
        examples.extend((isolated[i], overlap[i])
                        for i in sorted(isolated.keys() & overlap.keys() & paired_ids))
    changes = {}
    for field in ("mz_error_ppm", "rt_error_s", "area_ratio"):
        reference = 1 if field == "area_ratio" else 0
        increments = [abs(b[field] - reference) - abs(a[field] - reference)
                      for a, b in examples if a[field] is not None and b[field] is not None]
        changes[field] = median(increments) if increments else None
    return {"count": len(examples), "median_paired_change_in_absolute_error": changes}


def compare(calibration, density, executable, output, *, seeds=12, peaks=24, first_seed=81000,
            jobs=4, ppm=5., seconds=5.):
    if min(seeds, jobs) < 1 or peaks < 2:
        raise ValueError("Provide positive run and worker counts and at least two peaks per run.")
    validate_tolerances(ppm, seconds)
    exe, out = Path(executable).resolve(), Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    qseeds = list(range(first_seed, first_seed + seeds))

    def run(seed):
        directory = out / f"seed-{seed}"
        manifest = simulate(calibration, directory, seed=seed, peaks=peaks, density=density,
                            paired_control=True, point_truth=False)
        files = manifest["runs"][0]
        result = {"seed": seed, "pair_count": manifest["pair_count"]}
        for condition, source, expected in (("overlap", files["mzml"], files["truth"]),
                                            ("isolated", files["control_mzml"], files["control_truth"])):
            truth = json.loads((directory / expected).read_text())
            result[condition] = {}
            for model in KERNELS:
                detections = resolved_features(exe, directory / source, directory / f"{condition}-{model}", model)
                score = evaluate(truth, detections.pop("features"), ppm, seconds)
                result[condition][model] = dict(score, **detections)
        save(directory / "scores.json", result)
        print(f"Seed {seed} complete", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=jobs) as workers:
        results = list(workers.map(run, qseeds))
    total_pairs = sum(row["pair_count"] for row in results)
    if not total_pairs:
        raise ValueError("No overlap pairs were sampled from this calibration.")
    models = {}
    rng = random.Random(first_seed)
    for model in KERNELS:
        paired = []
        scores = {}
        for condition in ("isolated", "overlap"):
            subset = [row[condition][model] for row in results]
            scores[condition] = metrics(subset)
            scores[condition].update({key: sum(row[key] for row in subset) for key in
                                     ("omitted_centroid_components", "omitted_feature_components")})
        for row in results:
            n = row["pair_count"]
            counts = [sum(g["found"] == 2 for g in row[c][model]["pairs"]) for c in ("isolated", "overlap")]
            paired.append((counts[0] - counts[1], n))
        bootstrap = []
        for _ in range(2000):
            sample = rng.choices(paired, k=seeds)
            n = sum(b for _, b in sample)
            if n:
                bootstrap.append(sum(a for a, _ in sample) / n)
        bootstrap.sort()
        scores["pair_resolution_loss"] = sum(a for a, _ in paired) / total_pairs
        scores["paired_seed_bootstrap_95pct"] = [bootstrap[int(p * (len(bootstrap) - 1))] for p in (.025, .975)]
        scores["shared_match_errors"] = shared_match_errors(results, model)
        models[model] = scores
    report = {"seeds": qseeds, "peaks_per_seed": peaks, "pair_count": total_pairs,
              "ppm": ppm, "seconds": seconds, "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
              "models": models}
    save(out / "comparison.json", report)
    return report


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("calibration")
    cli.add_argument("density")
    cli.add_argument("--qalgo", dest="executable", required=True)
    cli.add_argument("-o", "--output", required=True)
    cli.add_argument("--seeds", type=int, default=12)
    cli.add_argument("--peaks", type=int, default=24)
    cli.add_argument("--first-seed", type=int, default=81000)
    cli.add_argument("--jobs", type=int, default=4)
    cli.add_argument("--ppm", type=float, default=5.)
    cli.add_argument("--seconds", type=float, default=5.)
    options = vars(cli.parse_args())
    try:
        compare(**options)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        cli.exit(1, f"{error}\n")
    print(Path(options["output"]).resolve() / "comparison.json")


if __name__ == "__main__":
    main()
