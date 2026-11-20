"""Place calibrated peak groups on regional m/z and retention-time grids."""

import csv
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
from statistics import median


def load_density(path, calibration, polarity):
    path = Path(path)
    profile = json.loads(path.read_text(encoding="utf-8"))
    digest = hashlib.sha256(calibration.read_bytes()).hexdigest()
    if profile.get("calibration_sha256") != digest:
        raise ValueError("Density measurements must come from this calibration. Regenerate density.json.")
    runs = [run for run in profile["runs"] if run["polarity"] == polarity]
    if not runs:
        raise ValueError(f"The density profile has no {polarity} runs.")
    mass_width, rt_width = (runs[0]["tile_width"][name] for name in ("mz", "rt_s"))
    counts = Counter()
    spacings = defaultdict(list)
    regional_groups = {mode: Counter() for mode in ("isolated", "mz", "rt")}
    histograms = {mode: Counter() for mode in ("mz", "rt")}
    observed_widths = defaultdict(list)
    for run in runs:
        if run["tile_width"] != runs[0]["tile_width"]:
            raise ValueError("Density runs must use the same grid.")
        with (path.parent / run["tile_file"]).open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                tile = (float(row["mz_start"]), float(row["rt_start_s"]))
                counts[tile] += int(row["resolved_features"])
                for mode, groups in regional_groups.items():
                    groups[tile] += int(row.get("groups_" + mode, 0))
        for mode, data in run["close_pairs"].items():
            q = data["separation_in_mean_fwhm"]
            if q:
                spacings[mode].append((data["count"], [q[key] for key in ("p10", "p50", "p90")]))
        for mode, data in run.get("simulation_groups", {}).get("pairs", {}).items():
            histograms[mode].update({i: n for i, n in enumerate(data["histogram_20"]) if n})
        if run.get("rt_width_file"):
            with (path.parent / run["rt_width_file"]).open(newline="", encoding="utf-8-sig") as stream:
                for row in csv.DictReader(stream):
                    key = int(row["feature"]), float(row["mz"]), float(row["rt_s"])
                    observed_widths[key].append(float(row["rt_fwhm_s"]))
    if not sum(counts.values()):
        raise ValueError("The density profile has no resolved features.")
    start = min(run["rt_range_s"][0] for run in runs)
    stop = max(run["rt_range_s"][1] for run in runs)
    return counts, (mass_width, rt_width), (start, stop), spacings, regional_groups, histograms, observed_widths


def mz_fwhm(window):
    skew = window["skew"]
    return (skew ** .5 + skew ** -.5) * window["fwhm"] / 2


def geometry(windows):
    start = min(w["mz"] - w["extent"] * w["step"] for w in windows)
    end = max(w["mz"] + w["extent"] * w["step"] for w in windows)
    return start, end


def layout(rng, templates, levels, dt, density, overlap, overlap_mode, make_peak):
    counts, (mass_width, rt_width), (rt_start, rt_end), spacings, observed, histograms, observed_widths = density
    by_tile = defaultdict(list)
    for row in templates:
        widths = observed_widths.get((row["id"], row["mz"], row["rt_s"]))
        if widths:
            row = dict(row, observed_rt_fwhm=median(widths))
        by_tile[(math.floor(row["mz"] / mass_width) * mass_width,
                 math.floor(row["rt_s"] / rt_width) * rt_width)].append(row)
    tiles = sorted(tile for tile in counts if counts[tile] > 0 and tile in by_tile
                   and max(tile[1], rt_start) < min(tile[1] + rt_width, rt_end))
    if not tiles:
        raise ValueError("No density tiles contain usable peak templates.")
    weights = [counts[tile] for tile in tiles]
    choices = {mode: [(tile, observed[mode][tile]) for tile in tiles if observed[mode][tile] > 0]
               for mode in observed}
    group_totals = {mode: sum(weight for _, weight in choices[mode]) for mode in choices}
    plan = []
    if overlap is None:
        if not any(group_totals.values()):
            raise ValueError("Density profile lacks measured feature groups. Regenerate density.json.")
        remaining = len(levels)
        while remaining:
            modes = ["isolated"]
            if remaining >= 2:
                modes += ["mz", "rt"] if overlap_mode == "mixed" else [overlap_mode]
            mode_weights = [group_totals[mode] for mode in modes]
            mode = rng.choices(modes, mode_weights)[0] if any(mode_weights) else "isolated"
            plan.append(mode)
            remaining -= 1 if mode == "isolated" else 2
    else:
        pairs = round(len(levels) * overlap / 2)
        plan = [overlap_mode] * pairs
        if overlap_mode == "mixed":
            plan = ["mz" if i % 2 == 0 else "rt" for i in range(pairs)]
        plan += ["isolated"] * (len(levels) - 2 * pairs)
        rng.shuffle(plan)
    groups = []
    next_id = 1
    for group_id, mode in enumerate(plan, start=1):
        regional_tiles, regional_weights = (zip(*choices[mode]) if overlap is None and choices[mode]
                                            else (tiles, weights))
        for _ in range(512):
            mass, time = rng.choices(regional_tiles, regional_weights)[0]
            center_mz = rng.uniform(mass, mass + mass_width)
            center_rt = rng.uniform(max(time, rt_start), min(time + rt_width, rt_end))
            first = make_peak(rng.choice(by_tile[(mass, time)]), center_mz, next_id, levels[next_id - 1])
            first.update(rt_s=center_rt, group=group_id, overlap_mode=mode, source_tile=(mass, time))
            windows = [first]
            if mode != "isolated":
                second = make_peak(rng.choice(by_tile[(mass, time)]), center_mz, next_id + 1, levels[next_id])
                second.update(rt_s=center_rt, group=group_id, overlap_mode=mode, source_tile=(mass, time))
                bins = histograms[mode]
                samples = spacings.get(mode)
                if bins:
                    bucket = rng.choices(list(bins), list(bins.values()))[0]
                    distance = (bucket + rng.random()) / 20
                elif samples:
                    options = rng.choices(samples, weights=[n for n, _ in samples])[0][1]
                    distance = rng.triangular(options[0], options[-1], options[1])
                else:
                    distance = .7
                distance = max(1e-6, min(.999999, distance))
                direction = rng.choice((-1, 1))
                if mode == "mz":
                    second["mz"] += direction * distance * (mz_fwhm(first) + mz_fwhm(second)) / 2
                else:
                    second["rt_s"] += direction * distance * (first["rt_halfwidth"] + second["rt_halfwidth"])
                    if not rt_start <= second["rt_s"] <= rt_end:
                        continue
                windows.append(second)
            start, end = geometry(windows)
            if start > 0 and all(end < a or b < start for a, b, _ in groups):
                groups.append((start, end, windows))
                next_id += len(windows)
                break
        else:
            raise ValueError("Could not place disjoint m/z groups. Try fewer peaks or a wider calibration range.")

    # Keep one signal-free segment, with sample spacing from the same device.
    background = make_peak(rng.choice(templates), max((b for _, b, _ in groups), default=100.) + .5,
                           None, None)
    if groups:
        right = max(b for _, b, _ in groups)
        background["mz"] = right + max(.5, background["extent"] * background["step"] * 2)
    start, end = geometry([background])
    groups.append((start, end, [background]))
    dt = max(dt, (rt_end - rt_start) / 2999)
    scans = max(2, math.ceil((rt_end - rt_start) / dt) + 1)
    times = [rt_start + i * dt for i in range(scans)]
    grids = []
    for start, end, windows in sorted(groups):
        step = max(min(w["step"] for w in windows), (end - start) / 240, end * 2e-7)
        # Keep the grid inside the interval used by the placement check.
        masses = [start + j * step for j in range(math.floor((end - start) / step) + 1)]
        if len(masses) < 5:
            raise ValueError("A mass group has too few representable sample points.")
        grids.append((masses, windows))
    return grids, times, dt, group_totals


def describe_group(group):
    first = group[0]
    data = {"id": first["group"], "mode": first["overlap_mode"], "source_tile": first["source_tile"],
            "features": [w["id"] for w in group]}
    if len(group) == 2:
        second = group[1]
        data.update(mz_distance_in_mean_fwhm=2 * abs(first["mz"] - second["mz"]) /
                    (mz_fwhm(first) + mz_fwhm(second)),
                    rt_distance_in_mean_fwhm=abs(first["rt_s"] - second["rt_s"]) /
                    (first["rt_halfwidth"] + second["rt_halfwidth"]))
    return data


def contributions(masses, windows, present, rt, shape):
    result = {}
    for w in windows:
        if w["id"] not in present:
            continue
        amplitude = w["height"] * shape((rt - w["rt_s"]) / w["rt_halfwidth"], w["family"])
        skew = math.sqrt(w["skew"])
        result[w["id"]] = [amplitude * shape(2 * (mz - w["mz"]) / w["fwhm"] *
                                             (skew if mz < w["mz"] else 1 / skew), w["family"])
                           for mz in masses]
    return result
