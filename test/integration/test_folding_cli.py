#!/usr/bin/env python3
"""End-to-end CLI smoke test for experimental profile-folding output.

Uses the repository's mzML fixture writer and supplies analytically known,
non-centroided sampled TOF-like profile peaks. This is an engineering test, not
proof of real-world mass accuracy, isolation or calibrated uncertainty.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import tempfile

import test_mzml_cli as fixtures


def folding_profile(
    polarity: str, rt: float, scan: int, blank: bool
) -> tuple[list[float], list[float]]:
    masses: list[float] = []
    intensities: list[float] = []
    # Constant TOF native bin spacing in sqrt(m/z), distinct local anchors,
    # unequal heights and shifted true apexes, with a fixed shared asymmetry.
    for j, mass in enumerate((109.21, 116.37, 125.19)):
        anchor = math.sqrt(mass)
        shift = 0.12 + j * 0.08
        log_height = 8.5 + j * 0.35 + 0.01 * math.sin(0.1 * scan)
        for k in range(-7, 8):
            coordinate = anchor + 0.00015 * k
            local = float(k) - shift
            log_signal = log_height + (-0.18 if local < 0 else -0.13) * local**2
            masses.append(coordinate**2)
            intensities.append(0.0 if blank else math.exp(log_signal))
    return masses, intensities


def read_rows(path: Path) -> list[dict[str, str]]:
    fixtures.require(path.is_file(), f"Missing folding CSV: {path}")
    with path.open(newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        required = {
            "spectrumIndex", "retentionTimeSeconds", "polarity",
            "axis", "model", "samplingStep", "profilesInFit", "anchorMz",
            "validApex", "apexMz", "apexUncertaintyMz", "b0", "b1",
            "b2", "b3", "pooledLogMSE", "degreesOfFreedom"
        }
        fixtures.require(required <= set(reader.fieldnames or []),
                         f"Missing folding columns: {reader.fieldnames}")
        return list(reader)


def run_cases(executable: Path, root: Path) -> None:
    original_profile = fixtures.profile
    fixtures.profile = folding_profile
    try:
        positive = root / "positive profile.mzML"
        fixtures.write_fixture(positive, indexed=True, mz_compressed=True,
                               intensity_compressed=True)
        output = root / "asymmetric"
        output.mkdir()
        cmd = ["-i", str(positive), "-o", str(output), "--peak-folding",
               "--folding-axis", "tof", "--folding-model", "asymmetric",
               "--folding-max-scans", "3", "--folding-min-intensity", "100"]
        fixtures.run_command(executable, cmd, root, expect_success=True,
                             log_prefix="asymmetric-")
        filename = output / "positive profile_folding.csv"
        rows = read_rows(filename)
        fixtures.require(len(rows) == 9, f"Expected 3 peaks x 3 scans, got {len(rows)}")
        for j, row in enumerate(rows):
            fixtures.require(row["axis"] == "tof" and row["model"] == "asymmetric",
                             "Incorrect folding settings in output")
            fixtures.require(row["validApex"] == "1", "Analytical peak lost fitted apex")
            fixtures.require(int(row["profilesInFit"]) == 3, "Peak pooling count incorrect")
            fixtures.require(int(row["spectrumIndex"]) == j // 3,
                             "Unexpected spectrum index selection")
            fixtures.require(float(row["apexUncertaintyMz"]) >= 0,
                             "Invalid analytical uncertainty")
            fixtures.require(abs(float(row["apexMz"]) -
                                 (109.21, 116.37, 125.19)[j % 3]) < 0.01,
                             "Incorrect fitted mass coordinate")

        # Existing folding CSV is not overwritten with --skip-existing.
        initial = filename.read_bytes()
        fixtures.run_command(executable, cmd + ["--skip-existing"], root,
                             expect_success=True, log_prefix="existing-")
        fixtures.require(filename.read_bytes() == initial,
                         "--skip-existing changed the folding CSV")

        symmetric = root / "symmetric"
        symmetric.mkdir()
        fixtures.run_command(executable, [
            "-i", str(positive), "-o", str(symmetric), "--peak-folding",
            "--folding-axis", "tof", "--folding-model", "symmetric",
            "--folding-max-scans", "2"
        ], root, expect_success=True, log_prefix="symmetric-")
        sym_rows = read_rows(symmetric / "positive profile_folding.csv")
        fixtures.require(len(sym_rows) == 6, "Symmetric mode lost candidates")
        for row in sym_rows:
            fixtures.require(abs(float(row["b2"]) - float(row["b3"])) < 1e-12,
                             "Symmetric mode has non-identical curvatures")

        orbitrap = root / "orbitrap"
        orbitrap.mkdir()
        # This tests coordinate-axis selection and code execution. The fixture
        # is physically TOF-like: the Orbitrap setting is intentionally not
        # presented as scientifically correct for the generated data.
        fixtures.run_command(executable, [
            "-i", str(positive), "-o", str(orbitrap), "--peak-folding",
            "--folding-axis", "orbitrap", "--folding-max-scans", "2"
        ], root, expect_success=True, log_prefix="orbitrap-")
        fixtures.require(bool(read_rows(orbitrap / "positive profile_folding.csv")),
                         "Orbitrap coordinate option produced no output")

        mixed = root / "mixed profile.mzML"
        fixtures.write_fixture(mixed, polarity="mixed")
        mixed_output = root / "mixed"
        mixed_output.mkdir()
        fixtures.run_command(executable, [
            "-i", str(mixed), "-o", str(mixed_output), "--peak-folding",
            "--folding-axis", "tof", "--folding-max-scans", "2"
        ], root, expect_success=True, log_prefix="mixed-")
        for polarity in ("positive", "negative"):
            rows = read_rows(mixed_output / f"mixed profile_{polarity}_folding.csv")
            fixtures.require(len(rows) == 6 and
                             all(r["polarity"] == polarity for r in rows),
                             "Mixed-polarity folding output is incorrectly separated")

        blank = root / "blank profile.mzML"
        fixtures.write_fixture(blank, blank=True)
        blank_output = root / "blank"
        blank_output.mkdir()
        fixtures.run_command(executable, [
            "-i", str(blank), "-o", str(blank_output), "--peak-folding",
            "--folding-axis", "tof"
        ], root, expect_success=True, log_prefix="blank-")
        fixtures.require(read_rows(blank_output / "blank profile_folding.csv") == [],
                         "Blank profile created false folding peaks")

        missing = root / "missing-axis"
        missing.mkdir()
        process = fixtures.run_command(executable, [
            "-i", str(positive), "-o", str(missing), "--peak-folding"
        ], root, expect_success=False, log_prefix="missing-")
        fixtures.require("--folding-axis" in process.stderr,
                         "Missing instrument coordinate option was not rejected")
        fixtures.require(not list(missing.glob("*_folding.csv")),
                         "Invalid flag created a folding output")

        incompatible = root / "incompatible"
        incompatible.mkdir()
        fixtures.run_command(executable, [
            "-i", str(positive), "-o", str(incompatible), "--peak-folding",
            "--folding-axis", "tof", "--printfeatures"
        ], root, expect_success=False, log_prefix="incompatible-")
        fixtures.require(not list(incompatible.glob("*_folding.csv")),
                         "Incompatible flags created a misleading output")

        # A folding-related option without --peak-folding must not be
        # silently ignored in favour of the standard feature workflow.
        standalone = root / "standalone"
        standalone.mkdir()
        for flag, value in (("--folding-model", "asymmetric"),
                            ("--folding-min-intensity", "100"),
                            ("--folding-axis", "tof")):
            process = fixtures.run_command(executable, [
                "-i", str(positive), "-o", str(standalone), flag, value
            ], root, expect_success=False, log_prefix="standalone-" + flag[2:] + "-")
            fixtures.require("--peak-folding" in process.stderr,
                             "Unused folding-specific option was not rejected")
        fixtures.require(not list(standalone.iterdir()),
                         "Standalone folding option unexpectedly created outputs")

        centroid = root / "centroid profile.mzML"
        fixtures.write_fixture(centroid, centroided=True)
        cent_output = root / "centroid"
        cent_output.mkdir()
        fixtures.run_command(executable, [
            "-i", str(centroid), "-o", str(cent_output), "--peak-folding",
            "--folding-axis", "tof"
        ], root, expect_success=False, log_prefix="centroid-")
        fixtures.require(not list(cent_output.glob("*_folding.csv")),
                         "Centroid-only mzML should not produce a profile folding CSV")

        for bad in ("0", "-1", "1.5"):
            fixtures.run_command(executable, [
                "-i", str(positive), "-o", str(missing), "--peak-folding",
                "--folding-axis", "tof", "--folding-max-scans", bad
            ], root, expect_success=False,
                log_prefix="bad-scans-" + bad.replace("-", "neg") + "-")
        print("PASS experimental folding CLI for synthetic mzML inputs")
    finally:
        fixtures.profile = original_profile


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    if not executable.is_file():
        parser.error(f"Binary missing: {executable}")
    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=True)
        work = args.work_dir
        run_cases(executable, work)
    else:
        with tempfile.TemporaryDirectory(prefix="qalgorithms-folding-cli-") as directory:
            run_cases(executable, Path(directory))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
