"""Measure local peak density and overlap from qAlgorithms calibration exports."""

import argparse
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
from heapq import heappop, heappush
import hashlib
from itertools import groupby
import json
import math
from pathlib import Path
from statistics import median

from benchmark import records, save


def quantiles(values):
    if not values:
        return {}
    values = sorted(values)
    return {f"p{p}": values[round((len(values) - 1) * p / 100)] for p in (10, 25, 50, 75, 90, 99)}


def describe_separations(values):
    counts = Counter(min(19, int(v * 20)) for v in values)
    return {"count": len(values), "separation_in_mean_fwhm": quantiles(values),
            "histogram_20": [counts[i] for i in range(20)]}


def overlaps(intervals):
    """Positive-length intersections of half-open intervals."""
    active, hits, pairs = [], set(), 0
    for start, end, identifier in sorted(intervals):
        if not end > start:
            continue
        while active and active[0][0] <= start:
            heappop(active)
        if active:
            hits.add(identifier)
            hits.update(partner for _, partner in active)
            pairs += len(active)
        heappush(active, (end, identifier))
    return pairs, hits


def observed_fwhm(times, areas, apex, low, high):
    first, end = bisect_left(times, low), bisect_right(times, high)
    nearest = bisect_left(times, apex, first, end)
    candidates = range(max(first, nearest - 1), min(nearest + 1, end))
    if not candidates:
        return None
    center = min(candidates, key=lambda i: abs(times[i] - apex))
    half = areas[center] / 2
    if half <= 0:
        return None
    left, right = center, center
    while left > first and areas[left] > half:
        left -= 1
    while right + 1 < end and areas[right] > half:
        right += 1
    if left == center or right == center or areas[left] > half or areas[right] > half:
        return None

    def crossing(outside, inside):
        fraction = (half - areas[outside]) / (areas[inside] - areas[outside])
        return times[outside] + fraction * (times[inside] - times[outside])

    width = crossing(right, right - 1) - crossing(left, left + 1)
    return width if width > 0 else None


def measure(folder, polarity, output, mass_width, rt_width, rt_width_method):
    scans = {int(row["scan"]): float(row["rt_s"]) for row in records(folder / "scans.csv")}
    exposure = Counter(math.floor(rt / rt_width) for rt in scans.values())
    regions = defaultdict(Counter)
    by_bin = defaultdict(lambda: array("d"))
    traces = defaultdict(lambda: (array("d"), array("d")))
    scan_counts, spacing, scaled_spacing = array("I"), array("d"), array("d")
    close_mass_spacing = array("d")
    totals = Counter()

    def cell(mass, rt):
        return math.floor(mass / mass_width), math.floor(rt / rt_width)

    with (folder / "centroids.csv").open(newline="", encoding="utf-8-sig") as stream:
        rows = csv.reader(stream)
        header = next(rows)
        expected = ("scan", "start", "end", "mz", "fwhm", "bin")
        columns = [header.index(name) for name in expected]
        scan_col, start_col, end_col, mass_col, width_col, bin_col = columns
        area_col = header.index("area")
        for scan, group in groupby(rows, key=lambda row: int(row[scan_col])):
            rt = scans[scan]
            centroids = []
            for row in group:
                mass, width, bin_id = float(row[mass_col]), float(row[width_col]), int(row[bin_col])
                centroids.append((mass, width, int(row[start_col]), int(row[end_col])))
                if bin_id >= 0:
                    regions[cell(mass, rt)]["binned_centroids"] += 1
                    if math.isfinite(width) and width > 0:
                        by_bin[bin_id].append(width)
                    if rt_width_method == "observed":
                        times, areas = traces[bin_id]
                        times.append(rt)
                        areas.append(float(row[area_col]))
            centroids.sort()
            scan_counts.append(len(centroids))
            totals["centroids"] += len(centroids)
            point_pairs, point_hits = overlaps((a, b, i) for i, (_, _, a, b) in enumerate(centroids))
            mass_pairs, mass_hits = overlaps((mz - width / 2, mz + width / 2, i)
                                            for i, (mz, width, _, _) in enumerate(centroids)
                                            if math.isfinite(width) and width > 0)
            totals["pair_profile_span_overlap"] += point_pairs
            totals["pair_mz_fwhm_overlap"] += mass_pairs
            totals["centroids_with_profile_span_overlap"] += len(point_hits)
            totals["centroids_with_mz_fwhm_overlap"] += len(mass_hits)
            totals["centroids_with_valid_fwhm"] += sum(math.isfinite(x[1]) and x[1] > 0 for x in centroids)
            for i, (mz, _, _, _) in enumerate(centroids):
                region = regions[cell(mz, rt)]
                region["centroids"] += 1
                region["centroids_profile_overlap"] += i in point_hits
                region["centroids_fwhm_overlap"] += i in mass_hits
            for a, b in zip(centroids, centroids[1:]):
                distance = b[0] - a[0]
                if a[0] > 0:
                    spacing.append(distance / a[0] * 1e6)
                if min(a[1], b[1]) > 0 and math.isfinite(a[1] + b[1]):
                    ratio = 2 * distance / (a[1] + b[1])
                    scaled_spacing.append(ratio)
                    if 0 < ratio < 1:
                        close_mass_spacing.append(ratio)
    n_scans_with_centroids = len(scan_counts)
    scan_counts.extend([0] * (len(scans) - n_scans_with_centroids))
    mz_widths = {key: median(widths) for key, widths in by_bin.items()}
    for key, (times, areas) in traces.items():
        if any(a > b for a, b in zip(times, times[1:])):
            ordered = sorted(zip(times, areas))
            traces[key] = tuple(array("d", values) for values in zip(*ordered))
    by_trace = defaultdict(list)
    features = []
    rt_fwhms = []
    observed_width_rows = []
    for row in records(folder / "profile_features.csv"):
        mz, rt = float(row["mz"]), float(row["retentionTime"])
        if not math.isfinite(mz) or not math.isfinite(rt):
            totals["invalid_feature_coordinates"] += 1
            continue
        region = regions[cell(mz, rt)]
        if row["isUnresolved"] != "0":
            totals["unresolved_descriptions"] += 1
            totals["unresolved_observed_maxima"] += int(row["unresolvedMaximaCount"])
            region["unresolved_descriptions"] += 1
            continue
        region["resolved_features"] += 1
        bin_id = int(row["binID"])
        low, high = float(row["lowestRetentionTime"]), float(row["highestRetentionTime"])
        width_rt = None
        if rt_width_method == "observed" and bin_id in traces:
            times, areas = traces[bin_id]
            # Feature CSV stores window bounds to four decimal places.
            width_rt = observed_fwhm(times, areas, rt, low - .000051, high + .000051)
            if width_rt is not None:
                observed_width_rows.append((row["ID"], mz, rt, width_rt))
        width_rt = width_rt if width_rt is not None else (high - low) / 3
        rt_fwhms.append(width_rt)
        by_trace[bin_id].append((low, high, len(features), rt))
        features.append((mz, rt, low, high, mz_widths.get(bin_id, 0), bin_id))
    same_trace_hits, same_trace_pairs = set(), 0
    rt_spacing = array("d")
    close_rt_spacing = array("d")
    candidates = []
    for trace in by_trace.values():
        pairs, hits = overlaps((low, high, i) for low, high, i, _ in trace)
        same_trace_hits.update(hits)
        same_trace_pairs += pairs
        by_apex = sorted((rt, rt_fwhms[i], i) for _, _, i, rt in trace)
        for (a, width_a, i), (b, width_b, j) in zip(by_apex, by_apex[1:]):
            rt_spacing.append(b - a)
            if width_a > 0 and width_b > 0:
                ratio = 2 * (b - a) / (width_a + width_b)
                if 0 < ratio < 1:
                    close_rt_spacing.append(ratio)
                    candidates.append((ratio, "rt", i, j))

    # m/z uses the bin's median centroid FWHM, whereas RT uses the fit window.
    # Exclude same-bin pairs since those already share a mass trace.
    active, simultaneous, cross_trace_hits = [], {}, set()
    cross_trace_pairs = 0
    intervals = sorted((mz - width / 2, mz + width / 2, i)
                       for i, (mz, _, low, high, width, _) in enumerate(features)
                       if width > 0 and high > low)
    for start, end, i in intervals:
        while active and active[0][0] <= start:
            _, other = heappop(active)
            del simultaneous[other]
        mz, apex, low, high, width, trace = features[i]
        for other, (other_low, other_high, other_trace) in simultaneous.items():
            if trace != other_trace and max(low, other_low) < min(high, other_high):
                cross_trace_pairs += 1
                cross_trace_hits.update((i, other))
                other_mz, other_rt, _, _, other_width, _ = features[other]
                ratio_mz = 2 * abs(mz - other_mz) / (width + other_width)
                ratio_rt = 2 * abs(apex - other_rt) / (rt_fwhms[i] + rt_fwhms[other])
                if 0 < ratio_mz < 1 and ratio_rt < 1:
                    candidates.append((ratio_mz, "mz", i, other))
        simultaneous[i] = low, high, trace
        heappush(active, (end, i))
    for name, hits in (("feature_same_bin_overlap", same_trace_hits),
                       ("feature_cross_bin_overlap", cross_trace_hits)):
        for i in hits:
            regions[cell(*features[i][:2])][name] += 1

    # The simulator accepts disjoint pairs. Count a feature at most once,
    # preferring pairs with the smallest separation relative to peak width.
    assigned = set()
    selected_spacings = {"mz": [], "rt": []}
    for ratio, mode, i, j in sorted(candidates):
        if i in assigned or j in assigned:
            continue
        assigned.update((i, j))
        selected_spacings[mode].append(ratio)
        anchor = min(i, j)
        regions[cell(*features[anchor][:2])]["groups_" + mode] += 1
    for i, feature in enumerate(features):
        if i not in assigned:
            regions[cell(*feature[:2])]["groups_isolated"] += 1

    keys = ("centroids", "binned_centroids", "centroids_profile_overlap", "centroids_fwhm_overlap",
            "resolved_features", "unresolved_descriptions", "feature_same_bin_overlap", "feature_cross_bin_overlap",
            "groups_isolated", "groups_mz", "groups_rt")
    tiles = []
    for (m, t), counts in sorted(regions.items()):
        tiles.append({"mz_start": m * mass_width, "rt_start_s": t * rt_width,
                      "selected_scans": exposure[t], **{key: counts[key] for key in keys}})
    with (output / "tiles.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("mz_start", "rt_start_s", "selected_scans", *keys))
        writer.writeheader()
        writer.writerows(tiles)
    if rt_width_method == "observed":
        with (output / "rt_fwhm.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("feature", "mz", "rt_s", "rt_fwhm_s"))
            writer.writerows(observed_width_rows)
    return {
        "polarity": polarity, "scans": len(scans), "rt_range_s": [min(scans.values()), max(scans.values())],
        "rt_width_method": rt_width_method, "features_with_measured_rt_fwhm": len(observed_width_rows),
        "tile_width": {"mz": mass_width, "rt_s": rt_width},
        "totals": dict(totals, resolved_features=len(features), pair_same_bin_rt_overlap=same_trace_pairs,
                       features_with_same_bin_rt_overlap=len(same_trace_hits),
                       features_with_mz_width_and_rt_window=len(intervals),
                       pair_cross_bin_mz_and_rt_overlap=cross_trace_pairs,
                       features_with_cross_bin_mz_and_rt_overlap=len(cross_trace_hits)),
        "centroids_per_scan": quantiles(scan_counts), "adjacent_centroid_spacing_ppm": quantiles(spacing),
        "adjacent_spacing_in_mean_fwhm": quantiles(scaled_spacing), "same_bin_apex_spacing_s": quantiles(rt_spacing),
        "close_pairs": {
            "mz": {"count": len(close_mass_spacing), "separation_in_mean_fwhm": quantiles(close_mass_spacing)},
            "rt": {"count": len(close_rt_spacing), "separation_in_mean_fwhm": quantiles(close_rt_spacing)},
        },
        "simulation_groups": {
            "isolated": len(features) - len(assigned),
            "paired_feature_fraction": len(assigned) / len(features) if features else 0.,
            "candidate_pairs": len(candidates),
            "pairs": {mode: describe_separations(values) for mode, values in selected_spacings.items()},
        },
        "densest_feature_tiles": sorted(tiles, key=lambda r: -r["resolved_features"])[:10],
        "densest_centroid_tiles": sorted((row for row in tiles if row["selected_scans"] >= 5),
                                        key=lambda r: -r["centroids"] / r["selected_scans"])[:10],
    }


def profile_density(calibration, output, mass_width=20, rt_width=30, rt_width_method="window"):
    calibration, output = Path(calibration), Path(output)
    if not all(math.isfinite(step) and step > 0 for step in (mass_width, rt_width)):
        raise ValueError("Tile widths must be finite and positive.")
    if rt_width_method not in ("window", "observed"):
        raise ValueError("Choose window or observed for the RT width method.")
    raw_profile = (calibration / "calibration.json").read_bytes()
    profile = json.loads(raw_profile)
    output.mkdir(parents=True, exist_ok=False)
    directories = sorted(calibration.glob("input-*"))
    report = []
    run = 0
    for directory in directories:
        for polarity in ("positive", "negative"):
            folder = directory / polarity
            if not (folder / "scans.csv").is_file():
                continue
            run_dir = output / (directory.name + "-" + polarity)
            run_dir.mkdir()
            summary = measure(folder, polarity, run_dir, mass_width, rt_width, rt_width_method)
            summary.update(input=profile["runs"][run]["input"], tile_file=str(run_dir.relative_to(output) / "tiles.csv"))
            if rt_width_method == "observed":
                summary["rt_width_file"] = str(run_dir.relative_to(output) / "rt_fwhm.csv")
            report.append(summary)
            print(f"{directory.name} {polarity} complete", flush=True)
            run += 1
    if not report or run != len(profile["runs"]):
        raise ValueError("Calibration profiles do not match the manifest.")
    path = output / "density.json"
    save(path, {"calibration": str(calibration.resolve()),
                "calibration_sha256": hashlib.sha256(raw_profile).hexdigest(), "runs": report})
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calibration", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--mass-width", type=float, default=20)
    parser.add_argument("--rt-width", type=float, default=30)
    parser.add_argument("--rt-width-method", choices=("window", "observed"), default="window")
    args = parser.parse_args()
    try:
        profile_density(args.calibration, args.output, args.mass_width, args.rt_width, args.rt_width_method)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
