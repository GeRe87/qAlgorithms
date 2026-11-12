"""Exercise the opt-in calibrated scan workflow on affine mzML profiles."""

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

from test_mzml_cli import NS, binary_array, element, write_fixture


def write_calibrated_fixture(path, hyperbolic=False):
    write_fixture(path, blank=True)
    tree = ET.parse(path)
    spectra = tree.find(f".//{{{NS}}}spectrumList")
    spectra.set("count", "3")
    for spectrum in list(spectra)[3:]:
        spectra.remove(spectrum)
    for i, spectrum in enumerate(spectra):
        mz = [100 + j / 64 for j in range(17)]
        raw = []
        for x in range(-8, 9):
            u = .22 * x
            shape = .9 * (1 - math.cosh(u)) + .06 * u - .09 * math.sinh(u) if hyperbolic else -.055 * x * x
            raw.append(800 * math.exp(shape + .003 * math.sin(2 * x)))
        if i == 1:
            mz[5] += .001
        if i == 2:
            raw = [0] * len(raw)
        spectrum.set("defaultArrayLength", str(len(mz)))
        spectrum.remove(spectrum.find(f"{{{NS}}}binaryDataArrayList"))
        arrays = element(spectrum, "binaryDataArrayList", count="2")
        binary_array(arrays, mz, 64, False, "mz")
        binary_array(arrays, raw, 32, True, "intensity")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def run_calibrated(host, path, output, model, *, step="1", area=".5", rotations="399", max_scale="8",
                   extra=(), expected=True):
    command = [str(host), "-i", str(path), "-o", str(output), "--calibrated-peaks",
               model, step, area, ".65", "7", max_scale, ".05", rotations, *extra]
    result = subprocess.run(command, text=True, capture_output=True)
    assert (result.returncode == 0) == expected, (result.stdout, result.stderr)
    return result


def check_output(path, model):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert all(None not in row and all(value is not None for value in row.values()) for row in rows), rows
    assert {row["model"] for row in rows} == {model}, rows
    scans = [row for row in rows if row["record"] == "scan"]
    assert [row["status"] for row in scans] == ["evaluated", "unsupportedAxis", "evaluated"], scans
    assert [int(row["familySize"]) for row in scans] == [4, 0, 4], scans
    assert float(scans[0]["backgroundAlphaPerWindow"]) == .05 / 4
    assert scans[2]["candidateCount"] == "0"
    measured = [row for row in rows if row["status"] == "measured"]
    assert len(measured) == 1, rows
    peak = measured[0]
    assert peak["spectrumIndex"] == scans[0]["spectrumIndex"] and peak["retentionTimeSeconds"] == scans[0]["retentionTimeSeconds"]
    assert abs(float(peak["mz"]) - 100.125) < .02 and float(peak["areaIntensityMz"]) > .5, peak
    variance_mz, covariance, variance_area = (float(peak[key]) for key in ("varianceMz", "covarianceMzArea", "varianceArea"))
    assert variance_mz > 0 and variance_area > 0 and covariance * covariance <= variance_mz * variance_area, peak
    return peak


def main():
    host = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="calibrated qPeaks ") as directory:
        root = Path(directory)
        for model, step in (("quadratic", "1"), ("hyperbolic", ".22")):
            path = root / f"{model}.mzML"
            write_calibrated_fixture(path, hyperbolic=model == "hyperbolic")
            output = root / model
            run_calibrated(host, path, output, model, step=step)
            csv_path = output / f"{model}_calibrated_peaks.csv"
            check_output(csv_path, f"builtin:{model}")
            assert list(output.iterdir()) == [csv_path]
            original = csv_path.read_bytes()
            run_calibrated(host, path, output, model, step=step, area="1e9", extra=("--skip-existing",))
            assert csv_path.read_bytes() == original
            run_calibrated(host, path, output, model, step=step, area="1e9")
            assert "measured" not in csv_path.read_text()
            prior = csv_path.read_bytes()
            insufficient = run_calibrated(host, path, output, model, step=step, rotations="10", expected=False)
            assert "Rotation budget cannot resolve" in insufficient.stderr
            staging = csv_path.with_suffix(".csv.partial")
            assert csv_path.read_bytes() == prior and not staging.exists()
            staging.mkdir()
            reserved = run_calibrated(host, path, output, model, step=step, expected=False)
            assert "partial output already exists" in reserved.stderr and csv_path.read_bytes() == prior
            assert staging.is_dir()
            staging.rmdir()
            run_calibrated(host, path, output, model, step=step, extra=("--printcentroids",), expected=False)
            if model == "quadratic":
                run_calibrated(host, path, output, model, max_scale="65536", expected=False)
    print("Calibrated affine mzML workflow passed for both built-in models")


if __name__ == "__main__":
    main()
