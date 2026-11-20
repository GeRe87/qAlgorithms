"""Verify profile calibration, reproducible mzML and independent scoring."""

import base64
import csv
import gzip
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET
import zlib

from test_mzml_cli import write_fixture

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "mzml_benchmark"))
from benchmark import build, matching, records, score, simulate
from compare_kernels import accepted, fits_budget
from measure_density import observed_fwhm, profile_density


def test(profiler, qalgo):
    with TemporaryDirectory(prefix="qAlgorithms benchmark ") as temp:
        root = Path(temp)
        inputs = [root / "real 1.mzML", root / "real 2.mzML"]
        for path in inputs:
            write_fixture(path, mz_compressed=True)
        out = root / "benchmark"
        build(inputs, out, profiler, replicate_inputs=True, seed=11, peaks=5, replicates=2, dropout=1,
              regional=False)
        calibration = out / "calibration"
        synthetic = out / "synthetic"
        report = json.loads((calibration / "calibration.json").read_text())
        assert len(report["runs"]) == 2
        for run in report["runs"]:
            assert run["reproducible_features"] == run["features"] > 0
            assert sum(run["point_classes"].values()) == 65 * 242
        assert all(int(r["end"]) > int(r["start"]) for r in
                   records(calibration / "input-001" / "positive" / "centroids.csv"))
        repeated = root / "same seed"
        simulate(calibration, repeated, 11, 5, 2, dropout=1)
        assert {p.name: p.read_bytes() for p in repeated.iterdir()} == {
            p.name: p.read_bytes() for p in synthetic.iterdir()}
        truth = synthetic / "synthetic-001-truth.json"
        assert len(json.loads(truth.read_text())) == 5
        assert json.loads((synthetic / "synthetic-002-truth.json").read_text()) == []

        mzml = synthetic / "synthetic-001.mzML"
        baseline = json.loads((synthetic / "manifest.json").read_text())["baseline"]
        spectrum = ET.parse(mzml).find(".//{*}spectrum")
        intensities = spectrum.findall(".//{*}binary")[1].text
        measured = tuple(x[0] for x in struct.iter_unpack("<d", zlib.decompress(base64.b64decode(intensities))))
        with gzip.open(synthetic / "synthetic-001-points.jsonl.gz", "rt") as stream:
            checked = 0
            for line in stream:
                row = json.loads(line)
                if row["scan"] != 0:
                    break
                assert row["start"] == checked
                for signal, noise in zip(row["signal"], row["noise"]):
                    assert math.isclose(signal + baseline + noise, measured[checked], rel_tol=1e-12, abs_tol=1e-12)
                    checked += 1
            assert checked == len(measured)
        for kernel, suffix in (("quadratic", ""), ("hyperbolic", "_hyperbolic")):
            found = root / kernel
            subprocess.run([qalgo, "-i", mzml, "-o", found, "--peak-model", kernel, "--benchmark-mode"],
                           check=True, capture_output=True)
            result = score(truth, found / ("synthetic-001" + suffix + "_features.csv"))
            assert result["tp"] + result["fn"] == 5
        found = root / "known.csv"
        first = json.loads(truth.read_text())[0]
        with found.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("id", "mz", "rt_s", "area"))
            writer.writerow((1, first["mz"], first["rt_s"], first["area"]))
            writer.writerow((2, 100000, 0, 1))
        result = score(truth, found)
        assert (result["tp"], result["fp"], result["fn"]) == (1, 1, 4)

        density = profile_density(calibration, root / "density")
        automatic = simulate(calibration, root / "automatic", 11, 4, density=density)
        assert automatic["requested_overlap_fraction"] is None
        assert automatic["pair_count"] == 0
        regional = root / "regional"
        manifest = simulate(calibration, regional, 11, 4, density=density, overlap=1, paired_control=True)
        assert manifest["pair_count"] == 2
        pairs = {group["mode"]: group for group in manifest["groups"]}
        assert 0 < pairs["mz"]["mz_distance_in_mean_fwhm"] < 1
        assert pairs["mz"]["rt_distance_in_mean_fwhm"] == 0
        assert 0 < pairs["rt"]["rt_distance_in_mean_fwhm"] < 1
        assert pairs["rt"]["mz_distance_in_mean_fwhm"] == 0
        seen = set()
        with gzip.open(regional / "synthetic-001-points.jsonl.gz", "rt") as stream:
            for line in stream:
                row = json.loads(line)
                if row["scan"] > 0:
                    break
                if len(row["components"]) == 2:
                    assert [sum(x) for x in zip(*row["components"].values())] == row["signal"]
                    seen.update(map(int, row["components"]))
        assert seen == set(range(1, 5))
        run = manifest["runs"][0]
        original = json.loads((regional / run["truth"]).read_text())
        control = json.loads((regional / run["control_truth"]).read_text())
        shifts = manifest["control_mass_shifts"]
        assert len(shifts) == 2
        for a, b in zip(original, control, strict=True):
            assert b == dict(a, mz=a["mz"] + shifts.get(a["id"], 0.))
        assert (regional / run["control_mzml"]).is_file()
        repeated = root / "regional repeat"
        simulate(calibration, repeated, 11, 4, density=density, overlap=1)
        assert (regional / "synthetic-001.mzML").read_bytes() == (repeated / "synthetic-001.mzML").read_bytes()
        regional_truth = regional / "synthetic-001-truth.json"
        perfect = root / "perfect.csv"
        with perfect.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("id", "mz", "rt_s", "area"), extrasaction="ignore")
            writer.writeheader()
            writer.writerows(json.loads(regional_truth.read_text()))
        paired_result = score(regional_truth, perfect)
        assert paired_result["pairs_present"] == paired_result["pairs_with_both_members_detected"] == 2

    # The nearest detection of the first peak is the only option of the third.
    rows = lambda coordinates: [{"mz": 100., "rt_s": rt} for rt in coordinates]
    assert len(matching(rows([1.1, 1.8, 3.7]), rows([0., 2., 4.]), 5, 2)) == 3
    null_runs = [[{"quality": .9}, {"quality": .7}], [{"quality": .7}]]
    threshold = fits_budget(null_runs, .5)
    assert [len(accepted(run, threshold)) for run in null_runs] == [1, 0]
    assert fits_budget([[], []], .5) is None
    assert observed_fwhm([0, 1, 3, 6, 9], [0, 2, 4, 2, 0], 3, 0, 9) == 5
    assert observed_fwhm([0, 1, 3, 6, 9], [0, 2, 4, 2, 0], 3, 2, 9) is None


if __name__ == "__main__":
    test(*map(Path, sys.argv[1:]))
