"""Calibrate from qAlgorithms profiles, simulate mzML, and score feature lists."""

import argparse
import base64
from collections import Counter, defaultdict
import csv
from contextlib import nullcontext
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
from statistics import median
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET
import zlib

from density_simulation import contributions, describe_group, layout, load_density, mz_fwhm


VERSION = 1
NS = "http://psi.hupo.org/ms/mzml"
CLASSES = ("outside_centroid", "centroid_only", "bin_only",
           "non_reproducible_feature", "reproducible_feature", "reproducibility_unknown")


def records(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        yield from csv.DictReader(stream)


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def matching(left, right, ppm, seconds):
    """Maximum-cardinality matching with nearest candidates tried first."""
    from bisect import bisect_left, bisect_right
    ordered = sorted((row["mz"], i) for i, row in enumerate(right))
    masses = [mass for mass, _ in ordered]
    neighbours = []
    for row in left:
        tolerance = row["mz"] * ppm * 1e-6
        lo = bisect_left(masses, row["mz"] - tolerance)
        hi = bisect_right(masses, row["mz"] + tolerance)
        candidates = [(math.hypot((mz - row["mz"]) / tolerance,
                                 (right[j]["rt_s"] - row["rt_s"]) / seconds), j)
                      for mz, j in ordered[lo:hi]
                      if abs(right[j]["rt_s"] - row["rt_s"]) <= seconds]
        neighbours.append([j for _, j in sorted(candidates)])
    chosen, owners = {}, {}
    for start in sorted(range(len(left)), key=lambda i: (len(neighbours[i]), i)):
        pending = [start]
        via = {}
        seen = {start}
        for i in pending:
            for j in neighbours[i]:
                if j in via:
                    continue
                via[j] = i
                if j not in owners:
                    while True:
                        i = via[j]
                        previous = chosen.get(i)
                        chosen[i], owners[j] = j, i
                        if previous is None:
                            break
                        j = previous
                    break
                owner = owners[j]
                if owner not in seen:
                    seen.add(owner)
                    pending.append(owner)
            else:
                continue
            break
    return sorted(chosen.items())


def classify(directory, feature_classes):
    from itertools import groupby
    associated = defaultdict(int)
    for row in records(directory / "feature_centroids.csv"):
        cid = int(row["centroid"])
        associated[cid] = max(associated[cid], feature_classes[int(row["feature"])])
    groups = iter(groupby(records(directory / "centroids.csv"), lambda row: int(row["scan"])))
    current = next(groups, None)
    counts = Counter()
    templates = {}
    feature_by_centroid = defaultdict(list)
    for row in records(directory / "feature_centroids.csv"):
        feature_by_centroid[int(row["centroid"])].append(int(row["feature"]))
    centroid_count = 0
    with (directory / "point_classes.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scan", "start", "end_exclusive", "classification"))
        for scan in records(directory / "scans.csv"):
            source_scan, size = int(scan["scan"]), int(scan["points"])
            events = defaultdict(list)
            if current and current[0] == source_scan:
                for row in current[1]:
                    centroid_count += 1
                    cid = int(row["id"])
                    group = associated.get(cid, 2 if int(row["bin"]) >= 0 else 1)
                    events[int(row["start"])].append((group, 1))
                    events[int(row["end"])].append((group, -1))
                    if cid in feature_by_centroid:
                        width, step, height = (float(row[key]) for key in ("fwhm", "mz_step", "height"))
                        skew = float(row["right_left_halfwidth"] or 1)
                        if all(math.isfinite(v) and v > 0 for v in (width, step, height, skew)):
                            for fid in feature_by_centroid[cid]:
                                if fid not in templates or templates[fid][2] < height:
                                    templates[fid] = (width, step, height, skew)
                current = next(groups, None)
            events[size]
            active = Counter()
            last = 0
            for end in sorted(events):
                if end > last:
                    group = max((g for g, count in active.items() if count > 0), default=0)
                    writer.writerow((source_scan, last, end, CLASSES[group]))
                    counts[CLASSES[group]] += end - last
                for group, delta in events[end]:
                    active[group] += delta
                last = end
    return counts, templates, centroid_count


def calibrate(inputs, output, profiler="qmzml-profile", kernel="quadratic",
              replicate_inputs=False, ppm=5, seconds=5):
    """Read profile MS1 files and create a calibration directory."""
    inputs = [Path(path).resolve() for path in inputs]
    validate_tolerances(ppm, seconds)
    if len(set(inputs)) != len(inputs):
        raise ValueError("An mzML input was selected more than once.")
    if not inputs or (replicate_inputs and len(inputs) < 2):
        raise ValueError("Select at least one input, or two distinct files for replicate matching.")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("Choose an empty calibration directory.")
    calibration = {"version": VERSION, "kernel": kernel, "benchmark_mode": True, "replicate_inputs": replicate_inputs,
                   "matching": {"ppm": ppm, "seconds": seconds}, "runs": [], "templates": [],
                   "background": {}, "relative_noise": {}, "rt_steps": {}}
    runs = []
    for number, path in enumerate(inputs, start=1):
        directory = output / f"input-{number:03}"
        directory.mkdir()
        with (directory / "profiler.log").open("wb") as log:
            result = subprocess.run([str(profiler), str(path), str(directory.resolve()), kernel, "--benchmark-mode"],
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise ValueError(f"qAlgorithms profiling failed. See {directory / 'profiler.log'}.")
        for polarity in ("positive", "negative"):
            side = directory / polarity
            if not (side / "scans.csv").is_file():
                continue
            features = [{"id": int(row["ID"]), "mz": float(row["mz"]),
                         "rt_s": float(row["retentionTime"]), "width_s": float(row["highestRetentionTime"])
                         - float(row["lowestRetentionTime"])}
                        for row in records(side / "profile_features.csv") if row["isUnresolved"] == "0"]
            runs.append((number, path, polarity, side, features))
    replicated = defaultdict(set)
    if replicate_inputs:
        for i, (number, _, polarity, _, left) in enumerate(runs):
            for j in range(i + 1, len(runs)):
                other_number, _, other_polarity, _, right = runs[j]
                if number == other_number or polarity != other_polarity:
                    continue
                for a, b in matching(left, right, ppm, seconds):
                    replicated[i].add(left[a]["id"])
                    replicated[j].add(right[b]["id"])

    for i, (number, path, polarity, side, features) in enumerate(runs):
        classes = {row["id"]: (4 if row["id"] in replicated[i] else 3) if replicate_inputs else 5
                   for row in features}
        counts, templates, centroid_count = classify(side, classes)
        for row in features:
            if row["id"] in templates and row["width_s"] > 0:
                width, step, height, skew = templates[row["id"]]
                calibration["templates"].append(dict(row, polarity=polarity, fwhm=width,
                                                     mz_step=step, height=height, skew=skew))
        scan_times = [float(row["rt_s"]) for row in records(side / "scans.csv")]
        calibration["rt_steps"].setdefault(polarity, []).extend(
            b - a for a, b in zip(scan_times, scan_times[1:]) if b > a)
        with (side / "feature_reproducibility.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("feature", "classification"))
            writer.writerows((fid, CLASSES[group]) for fid, group in classes.items())
        bins = {row["bin"] for row in records(side / "bins.csv")}
        unresolved = sum(row["isUnresolved"] != "0" for row in records(side / "profile_features.csv"))
        calibration["runs"].append({
            "input": str(path), "polarity": polarity, "scans": len(scan_times),
            "centroids": centroid_count, "bins": len(bins), "features": len(features), "unresolved": unresolved,
            "reproducible_features": len(replicated[i]) if replicate_inputs else None,
            "point_classes": {name: counts[name] for name in CLASSES},
        })
        calibration["background"].setdefault(polarity, []).extend(
            float(row["intensity"]) for row in records(side / "background.csv"))
        samples = [float(row["log_residual_proxy"]) for row in records(side / "relative_noise.csv")]
        if samples:
            rng = random.Random(number)
            calibration["relative_noise"].setdefault(polarity, []).extend(
                rng.sample(samples, min(16384, len(samples))))
    save(output / "calibration.json", calibration)
    return calibration


def validate_tolerances(ppm, seconds):
    if not all(math.isfinite(v) and v > 0 for v in (ppm, seconds)):
        raise ValueError("Matching tolerances must be finite and positive.")


def cv(parent, accession, name, value="", **attributes):
    return ET.SubElement(parent, "cvParam", cvRef="MS", accession=accession, name=name,
                         value=str(value), **attributes)


def mzml_document(scans):
    root = ET.Element("mzML", xmlns=NS, id="qAlgorithms-benchmark", version="1.1.0")
    vocabularies = ET.SubElement(root, "cvList", count="2")
    ET.SubElement(vocabularies, "cv", id="MS", fullName="PSI-MS", URI="https://purl.obolibrary.org/obo/ms.obo")
    ET.SubElement(vocabularies, "cv", id="UO", fullName="Unit Ontology", URI="https://purl.obolibrary.org/obo/uo.obo")
    file_contents = ET.SubElement(ET.SubElement(root, "fileDescription"), "fileContent")
    cv(file_contents, "MS:1000579", "MS1 spectrum")
    software = ET.SubElement(ET.SubElement(root, "softwareList", count="1"), "software",
                             id="benchmark", version=str(VERSION))
    cv(software, "MS:1000799", "custom unreleased software tool", "qAlgorithms mzML benchmark")
    instruments = ET.SubElement(root, "instrumentConfigurationList", count="1")
    ET.SubElement(instruments, "instrumentConfiguration", id="IC")
    processing = ET.SubElement(ET.SubElement(root, "dataProcessingList", count="1"), "dataProcessing", id="DP")
    method = ET.SubElement(processing, "processingMethod", order="0", softwareRef="benchmark")
    cv(method, "MS:1000544", "Conversion to mzML")
    run = ET.SubElement(root, "run", id="run", defaultInstrumentConfigurationRef="IC")
    spectra = ET.SubElement(run, "spectrumList", count=str(scans), defaultDataProcessingRef="DP")
    return root, spectra


def add_scan(spectra, index, time, polarity, masses, intensities):
    spectrum = ET.SubElement(spectra, "spectrum", index=str(index), id=f"scan={index+1}",
                             defaultArrayLength=str(len(masses)))
    cv(spectrum, "MS:1000511", "ms level", "1")
    cv(spectrum, "MS:1000128", "profile spectrum")
    cv(spectrum, "MS:1000130" if polarity == "positive" else "MS:1000129", f"{polarity} scan")
    scan = ET.SubElement(ET.SubElement(spectrum, "scanList", count="1"), "scan")
    cv(scan, "MS:1000016", "scan start time", time, unitCvRef="UO",
       unitAccession="UO:0000010", unitName="second")
    arrays = ET.SubElement(spectrum, "binaryDataArrayList", count="2")
    for values, accession, name, vocab, unit, unit_name in (
        (masses, "MS:1000514", "m/z array", "MS", "MS:1000040", "m/z"),
        (intensities, "MS:1000515", "intensity array", "MS", "MS:1000131", "number of detector counts"),
    ):
        binary = base64.b64encode(zlib.compress(struct.pack(f"<{len(values)}d", *values))).decode("ascii")
        array = ET.SubElement(arrays, "binaryDataArray", encodedLength=str(len(binary)))
        cv(array, "MS:1000523", "64-bit float")
        cv(array, "MS:1000574", "zlib compression")
        cv(array, accession, name, unitCvRef=vocab, unitAccession=unit, unitName=unit_name)
        ET.SubElement(array, "binary").text = binary


def shape(offset, family):
    if family == "gaussian":
        return math.exp(-math.log(2) * offset * offset)
    tail = math.exp(-math.acosh(2) * abs(offset))
    return 2 * tail / (1 + tail * tail)


def integrate(x, y):
    return sum((b - a) * (u + v) / 2 for a, b, u, v in zip(x, x[1:], y, y[1:]))


def add_noise(signal, baseline, additive, multipliers):
    return [max(0., baseline + v + (n + v * e))
            for v, n, e in zip(signal, additive, multipliers)]


def simulate(calibration, output, seed=1, peaks=12, replicates=1, polarity=None, dropout=0.,
             strengths=(1.,), point_truth=True, density=None, overlap=None, overlap_mode="mixed",
             paired_control=False):
    """Generate profile mzML and its sampled noise-free truth."""
    if peaks < 0 or replicates < 1 or not 0 <= dropout <= 1:
        raise ValueError("Choose a nonnegative peak count, at least one run and dropout between 0 and 1.")
    if not strengths or not all(math.isfinite(value) and value > 0 for value in strengths):
        raise ValueError("Peak strengths must be finite and positive.")
    if (overlap is not None and not 0 <= overlap <= 1) or overlap_mode not in ("mixed", "mz", "rt"):
        raise ValueError("Choose an overlap fraction between 0 and 1 and a valid overlap mode.")
    if overlap and not density:
        raise ValueError("Overlap placement needs a density.json from measure_density.py.")
    if paired_control and not density:
        raise ValueError("A paired control requires the regional density profile.")
    path = Path(calibration)
    if path.is_dir():
        path /= "calibration.json"
    profile = json.loads(path.read_text(encoding="utf-8"))
    if profile["version"] != VERSION:
        raise ValueError("Unsupported calibration version.")
    available = profile["background"].keys()
    polarity = polarity or next((p for p in available if any(t["polarity"] == p for t in profile["templates"])),
                                next(iter(available), None))
    templates = [t for t in profile["templates"] if t["polarity"] == polarity]
    if polarity not in available or (peaks and not templates):
        raise ValueError("The requested polarity has no usable calibration templates for peaks.")
    steps = profile["rt_steps"].get(polarity, [])
    if not steps:
        raise ValueError("Calibration needs at least two distinct retention times.")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Choose an empty simulation directory.")
    rng = random.Random(seed)
    dt = median(steps)
    real_background = profile["background"][polarity] or [0.]
    baseline = median(real_background)
    noise = [v - baseline for v in real_background]
    relative_samples = profile["relative_noise"].get(polarity) or [0.]
    center = median(relative_samples)
    relative = [min(2., max(-2., value - center)) for value in relative_samples]
    windows = []
    masses = [t["mz"] for t in templates]
    mass_min, mass_max = (min(masses), max(masses)) if masses else (100., 500.)
    background_masses = (mass_min + mass_max) / 2 + 10
    levels = [strengths[i % len(strengths)] for i in range(peaks)]
    if len(strengths) > 1:
        rng.shuffle(levels)
    def make_peak(t, new_mass, identifier, strength):
        width = t["fwhm"]
        # Bound the segment size and use a spacing that survives the float axis
        # used by the production centroiding code.
        step = max(t["mz_step"], width / 40, abs(new_mass) * 2e-7)
        extent = max(8, min(120, math.ceil(4 * width / step)))
        rt_halfwidth = (t["observed_rt_fwhm"] / 2 if "observed_rt_fwhm" in t
                        else max(2 * dt, t["width_s"] / 6))
        return {"id": identifier, "mz": new_mass, "extent": extent, "step": step,
                "fwhm": width, "height": t["height"] * (strength if identifier else 1),
                "strength": strength, "rt_halfwidth": rt_halfwidth,
                "rt_width_source": "observed" if "observed_rt_fwhm" in t else "fit_window",
                "skew": min(3, max(1/3, t["skew"])), "family": rng.choice(("gaussian", "sech"))}

    if density:
        if not templates:
            raise ValueError("Regional sampling needs peak templates, even for a background-only run.")
        measured = load_density(density, path, polarity)
        grids, times, dt, observed_groups = layout(rng, templates, levels, dt, measured,
                                                   overlap, overlap_mode, make_peak)
        windows = [w for _, group in grids for w in group]
    else:
        for i in range(peaks + 1):
            t = rng.choice(templates) if templates else {"fwhm": .01, "mz_step": .001, "height": 1,
                                                         "width_s": 10, "skew": 1}
            new_mass = rng.uniform(mass_min, mass_max) if i < peaks else background_masses
            peak = make_peak(t, new_mass, i + 1 if i < peaks else None, levels[i] if i < peaks else None)
            step, extent = peak["step"], peak["extent"]
            while any(abs(new_mass - w["mz"]) <= (extent * step + w["extent"] * w["step"])
                      for w in windows):
                new_mass += max(0.1, 2 * extent * step)
            peak["mz"] = new_mass
            windows.append(peak)
        duration = max(w["rt_halfwidth"] for w in windows) * 16
        dt = max(dt, duration / 1999)
        times = [i * dt for i in range(max(21, math.ceil(duration / dt) + 1))]
        for w in windows:
            margin = min(7 * w["rt_halfwidth"], times[-1] / 2)
            peak_scan = round(rng.uniform(margin, times[-1] - margin) / dt)
            w["rt_s"] = times[peak_scan]
        windows.sort(key=lambda w: w["mz"])
        grids = [([w["mz"] + j * w["step"] for j in range(-w["extent"], w["extent"] + 1)], [w])
                 for w in windows]
    manifest = {"version": VERSION, "seed": seed, "polarity": polarity, "source_kernel": profile["kernel"],
                "python_version": sys.version.split()[0], "calibration_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "background_samples": len(profile["background"][polarity]), "baseline": baseline,
                "rt_step_s": dt, "dropout": dropout,
                "features": [{k: w[k] for k in ("id", "mz", "rt_s", "family")} for w in windows if w["id"]], "runs": []}
    if density:
        manifest.update(density_sha256=hashlib.sha256(Path(density).read_bytes()).hexdigest(),
                        requested_overlap_fraction=overlap, overlap_mode=overlap_mode,
                        measured_placement_groups=observed_groups,
                        pair_count=sum(len(group) == 2 for _, group in grids),
                        groups=[describe_group(group) for _, group in grids if group[0]["id"] is not None])
    shifts = {}
    if paired_control:
        # Reserve disjoint mass intervals for the relocated partners. Use the
        # same sampled shape and noise draws in the original and shifted lanes.
        right = max(masses[-1] for masses, _ in grids)
        for masses, group in grids:
            if len(group) == 2:
                left = right + max(.1, right * 2e-5)
                shifts[group[1]["id"]] = left - masses[0]
                right = left + masses[-1] - masses[0]
        manifest["control_mass_shifts"] = shifts
    for run in range(1, replicates + 1):
        stem = f"synthetic-{run:03}"
        tree, spectra = mzml_document(len(times))
        if paired_control:
            control_tree, control_spectra = mzml_document(len(times))
        scan_areas = defaultdict(list)
        present = {w["id"] for w in windows if w["id"] and (run == 1 or dropout == 0 or rng.random() >= dropout)}
        target_path = out / (stem + "-points.jsonl.gz")
        with (target_path.open("wb") if point_truth else nullcontext()) as target:
            with (gzip.GzipFile(fileobj=target, mode="wb", filename="", mtime=0)
                  if point_truth else nullcontext()) as truth_points:
                for index, time in enumerate(times):
                    mz_all, intensity_all = [], []
                    control_mz, control_intensity, relocated = [], [], []
                    for masses, group in grids:
                        if density:
                            components = contributions(masses, group, present, time, shape)
                            signal = list(map(sum, zip(*components.values()))) if components else [0.] * len(masses)
                        else:
                            w = group[0]
                            offsets = range(-w["extent"], w["extent"] + 1)
                            amp = w["height"] * shape((time - w["rt_s"]) / w["rt_halfwidth"], w["family"])
                            signal = [amp * shape(2 * j * w["step"] / w["fwhm"] *
                                                  (math.sqrt(w["skew"]) if j < 0 else 1 / math.sqrt(w["skew"])),
                                                  w["family"]) if w["id"] in present else 0. for j in offsets]
                            components = {w["id"]: signal} if w["id"] in present else {}
                        additive = rng.choices(noise, k=len(signal))
                        multipliers = [math.expm1(e) for e in rng.choices(relative, k=len(signal))]

                        values = add_noise(signal, baseline, additive, multipliers)
                        if paired_control:
                            control_mz.extend(masses)
                            if len(group) == 2:
                                zeros = [0.] * len(masses)
                                control_intensity.extend(add_noise(components.get(group[0]["id"], zeros),
                                                                   baseline, additive, multipliers))
                                shift = shifts[group[1]["id"]]
                                relocated.append(([mass + shift for mass in masses],
                                                  add_noise(components.get(group[1]["id"], zeros),
                                                            baseline, additive, multipliers)))
                            else:
                                control_intensity.extend(values)
                        noise_draw = [v - (baseline + clean) for v, clean in zip(values, signal)]
                        if truth_points:
                            row = {"scan": index, "start": len(mz_all),
                                   "feature": next(iter(components)) if len(components) == 1 else None,
                                   "signal": signal, "noise": noise_draw}
                            if density:
                                row["components"] = components
                            truth_points.write((json.dumps(row, separators=(",", ":")) + "\n").encode())
                        for identifier, points in components.items():
                            scan_areas[identifier].append(integrate(masses, points))
                        mz_all.extend(masses)
                        intensity_all.extend(values)
                    add_scan(spectra, index, time, polarity, mz_all, intensity_all)
                    if paired_control:
                        for masses, values in relocated:
                            control_mz.extend(masses)
                            control_intensity.extend(values)
                        add_scan(control_spectra, index, time, polarity, control_mz, control_intensity)
        ET.ElementTree(tree).write(out / (stem + ".mzML"), encoding="utf-8", xml_declaration=True)
        truth = [{"id": w["id"], "mz": w["mz"], "rt_s": w["rt_s"], "family": w["family"],
                  "mz_fwhm": mz_fwhm(w),
                  "rt_fwhm_s": 2 * w["rt_halfwidth"], "mz_right_left_halfwidth": w["skew"],
                  "height": w["height"], "strength": w["strength"],
                  "area": integrate(times, scan_areas[w["id"]])}
                 for w in windows if w["id"] in present]
        if density:
            by_id = {w["id"]: w for w in windows}
            for peak in truth:
                w = by_id[peak["id"]]
                peak.update(group=w["group"], overlap_mode=w["overlap_mode"],
                            rt_width_source=w["rt_width_source"])
        save(out / (stem + "-truth.json"), truth)
        entry = {"mzml": stem + ".mzML", "truth": stem + "-truth.json",
                 "points": target_path.name if point_truth else None}
        if paired_control:
            control = stem + "-isolated"
            ET.ElementTree(control_tree).write(out / (control + ".mzML"), encoding="utf-8", xml_declaration=True)
            save(out / (control + "-truth.json"),
                 [dict(peak, mz=peak["mz"] + shifts.get(peak["id"], 0.)) for peak in truth])
            entry.update(control_mzml=control + ".mzML", control_truth=control + "-truth.json")
        manifest["runs"].append(entry)
    save(out / "manifest.json", manifest)
    return manifest


def score(truth, detected, ppm=5, seconds=5):
    """Score resolved CSV features with m/z, retention time in seconds, and optional area."""
    validate_tolerances(ppm, seconds)
    expected = json.loads(Path(truth).read_text(encoding="utf-8"))
    found = []
    for i, row in enumerate(records(Path(detected)), start=1):
        if row.get("isUnresolved", "0") != "0":
            continue
        mass, rt = float(row["mz"]), float(row.get("retentionTime") or row["rt_s"])
        if not math.isfinite(mass) or mass <= 0 or not math.isfinite(rt):
            raise ValueError(f"Invalid coordinates in detection row {i}.")
        found.append({"id": row.get("ID") or row.get("id") or i, "mz": mass, "rt_s": rt,
                      "area": float(row["area"]) if row.get("area") else None})
    return score_features(expected, found, ppm, seconds)


def score_features(expected, found, ppm, seconds):
    pairs = []
    matches = matching(expected, found, ppm, seconds)
    for i, j in matches:
        a, b = expected[i], found[j]
        area_ratio = b["area"] / a["area"] if b["area"] is not None and math.isfinite(b["area"]) else None
        pairs.append({"truth": a["id"], "detected": b["id"], "mz_error_ppm": (b["mz"] - a["mz"]) / a["mz"] * 1e6,
                      "rt_error_s": b["rt_s"] - a["rt_s"], "area_ratio": area_ratio})
    selected_left, selected_right = {i for i, _ in matches}, {j for _, j in matches}
    result = {"tp": len(matches), "fp": len(found) - len(matches), "fn": len(expected) - len(matches),
              "matching": {"ppm": ppm, "seconds": seconds}, "matches": pairs,
              "false_positives": [row["id"] for j, row in enumerate(found) if j not in selected_right],
              "false_negatives": [row["id"] for i, row in enumerate(expected) if i not in selected_left]}
    if any("group" in peak for peak in expected):
        groups = {}
        for i, peak in enumerate(expected):
            group = groups.setdefault(peak["group"], {"mode": peak["overlap_mode"],
                                                      "present": [], "detected": []})
            group["present"].append(peak["id"])
            if i in selected_left:
                group["detected"].append(peak["id"])
        result["groups"] = groups
        result["pairs_with_both_members_detected"] = sum(
            len(group["present"]) == len(group["detected"]) == 2 for group in groups.values())
        result["pairs_present"] = sum(len(group["present"]) == 2 for group in groups.values())
    return result


def build(inputs, output, profiler="qmzml-profile", kernel="quadratic",
          replicate_inputs=False, ppm=5, seconds=5, seed=1, peaks=12, replicates=1, polarity=None,
          dropout=0., regional=True, overlap=None, overlap_mode="mixed"):
    """Calibrate and generate synthetic runs in one operation."""
    output = Path(output)
    if peaks < 0 or replicates < 1 or not 0 <= dropout <= 1:
        raise ValueError("Choose a nonnegative peak count, at least one run and dropout between 0 and 1.")
    if (overlap is not None and not 0 <= overlap <= 1) or overlap_mode not in ("mixed", "mz", "rt"):
        raise ValueError("Choose an overlap fraction between 0 and 1 and a valid overlap mode.")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose an empty benchmark directory.")
    calibrate(inputs, output / "calibration", profiler, kernel, replicate_inputs, ppm, seconds)
    density = None
    if (regional or overlap) and peaks:
        from measure_density import profile_density
        density = profile_density(output / "calibration", output / "density")
    return simulate(output / "calibration", output / "synthetic", seed, peaks, replicates, polarity, dropout,
                    density=density, overlap=overlap, overlap_mode=overlap_mode)


def default_profiler():
    name = "qmzml-profile.exe" if os.name == "nt" else "qmzml-profile"
    path = Path(__file__).resolve().parents[2] / "build"
    return (os.environ.get("QALGORITHMS_PROFILER") or shutil.which(name)
            or next((str(x / name) for x in (path, path / "Release", path / "Debug")
                     if (x / name).is_file()), name))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    c = commands.add_parser("calibrate")
    b = commands.add_parser("build")
    for command in (c, b):
        command.add_argument("inputs", nargs="+")
        command.add_argument("-o", "--output", required=True)
        command.add_argument("--profiler", default=default_profiler())
        command.add_argument("--kernel", choices=("quadratic", "hyperbolic"), default="quadratic")
        command.add_argument("--replicate-inputs", action="store_true")
    g = commands.add_parser("simulate")
    g.add_argument("calibration")
    g.add_argument("-o", "--output", required=True)
    g.add_argument("--density", help="density.json from measure_density.py")
    b.add_argument("--regional", action=argparse.BooleanOptionalAction, default=True,
                   help="Measure and use regional peak density (default)")
    for command in (g, b):
        command.add_argument("--seed", type=int, default=1)
        command.add_argument("--peaks", type=int, default=12)
        command.add_argument("--replicates", type=int, default=1)
        command.add_argument("--polarity", choices=("positive", "negative"))
        command.add_argument("--dropout", type=float, default=0.)
        command.add_argument("--overlap", type=float,
                             help="Override the measured overlap with an experimental paired fraction")
        command.add_argument("--overlap-mode", choices=("mixed", "mz", "rt"), default="mixed")
    s = commands.add_parser("score")
    s.add_argument("truth")
    s.add_argument("detected")
    s.add_argument("-o", "--output")
    for command in (c, b, s):
        command.add_argument("--ppm", type=float, default=5)
        command.add_argument("--seconds", type=float, default=5)
    commands.add_parser("gui")
    args = vars(parser.parse_args())
    command = args.pop("command")
    try:
        if command == "gui":
            from gui import launch
            launch()
        elif command == "score":
            output = args.pop("output")
            result = score(**args)
            if output:
                save(Path(output), result)
            print(json.dumps({k: result[k] for k in ("tp", "fp", "fn")} if output else result))
        else:
            {"calibrate": calibrate, "simulate": simulate, "build": build}[command](**args)
            print(Path(args["output"]).resolve())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
