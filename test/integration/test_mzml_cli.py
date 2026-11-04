#!/usr/bin/env python3
"""Small deterministic mzML -> feature CSV integration test (Python stdlib only).

The generated profiles contain two isolated synthetic Gaussian signals with
known mass and retention-time centres. This checks the complete executable,
binary-array decoding, polarity selection, coordinates, and CSV output. It is
an engineering smoke test, not validation of analytical accuracy on real data.

Usage: python3 test/integration/test_mzml_cli.py /path/to/qAlgorithms
Use --work-dir DIR to retain generated inputs, outputs, and process logs.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import math
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import zlib


NS = "http://psi.hupo.org/ms/mzml"
ET.register_namespace("", NS)
EXPECTED = {
    "positive": ((100.125, 22.25), (150.375, 31.4)),
    "negative": ((101.25, 23.25), (151.5, 32.4)),
}
REQUIRED_COLUMNS = {
    "ID", "binID", "mz", "mzUncertainty", "retentionTime",
    "RT_Uncertainty", "area", "height", "DQSC", "DQSB", "DQSF",
}
PROFILE_CASES = ("positive", "negative", "mixed", "blank", "malformed", "centroided", "ms2_only")
CLI_CASES = ("missing_output_value", "output_collision", "skip_error", "default_features",
             "unresolved_trace_write_failure", "orphan_unresolved_trace", "orphan_feature_output", "collective_output")


class TestFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise TestFailure(message)


def element(parent: ET.Element, tag: str, **attributes: str) -> ET.Element:
    return ET.SubElement(parent, f"{{{NS}}}{tag}", attributes)


def cv(parent: ET.Element, accession: str, name: str,
       value: str = "", **attributes: str) -> None:
    element(parent, "cvParam", cvRef="MS", accession=accession,
            name=name, value=value, **attributes)


def binary_array(parent: ET.Element, values: list[float], bits: int,
                 compressed: bool, kind: str) -> None:
    fmt = "d" if bits == 64 else "f"
    data = struct.pack(f"<{len(values)}{fmt}", *values)
    if compressed:
        data = zlib.compress(data)
    encoded = base64.b64encode(data).decode("ascii")
    array = element(parent, "binaryDataArray", encodedLength=str(len(encoded)))
    cv(array, "MS:1000523" if bits == 64 else "MS:1000521",
       f"{bits}-bit float")
    cv(array, "MS:1000574" if compressed else "MS:1000576",
       "zlib compression" if compressed else "no compression")
    if kind == "mz":
        cv(array, "MS:1000514", "m/z array", unitCvRef="MS",
           unitAccession="MS:1000040", unitName="m/z")
    else:
        cv(array, "MS:1000515", "intensity array", unitCvRef="MS",
           unitAccession="MS:1000131", unitName="number of detector counts")
    element(array, "binary").text = encoded


def profile(polarity: str, rt: float, scan: int,
            blank: bool) -> tuple[list[float], list[float]]:
    masses: list[float] = []
    intensities: list[float] = []
    for peak_number, (mass, apex_rt) in enumerate(EXPECTED[polarity]):
        sigma_rt = 2.0 + 0.2 * peak_number
        height = (150000.0 - peak_number * 30000.0) * math.exp(
            -0.5 * ((rt - apex_rt) / sigma_rt) ** 2)
        # Seeded point-wise intensity noise gives a nonzero regression residual.
        # The very small mass drift stays below the uncertainty from that noise;
        # a drifting but almost exact curve would not form a common mass bin.
        noise = random.Random(123456 + scan * 101 + peak_number)
        mass_shift = 0.000005 * math.sin(0.55 * scan + peak_number)
        for point in range(121):
            mz = mass + (point - 60) * 0.001
            deviation = 1.0 + 0.02 * noise.uniform(-1.0, 1.0)
            intensity = height * math.exp(
                -0.5 * ((mz - mass - mass_shift) / 0.005) ** 2)
            masses.append(mz)
            intensities.append(0.0 if blank else intensity * deviation)
    return masses, intensities


def write_fixture(path: Path, *, polarity: str = "positive",
                  indexed: bool = False, mz_bits: int = 64,
                  intensity_bits: int = 32, mz_compressed: bool = False,
                  intensity_compressed: bool = False,
                  reverse_arrays: bool = False, rt_minutes: bool = False,
                  blank: bool = False, centroided: bool = False,
                  ms_level: int = 1, intensity_scale: float = 1.0,
                  collective: bool = False) -> None:
    root = ET.Element(f"{{{NS}}}mzML", {"id": "synthetic-smoke", "version": "1.1.0"})
    cv_list = element(root, "cvList", count="2")
    element(cv_list, "cv", id="MS", fullName="PSI-MS",
            URI="https://purl.obolibrary.org/obo/ms.obo")
    element(cv_list, "cv", id="UO", fullName="Unit Ontology",
            URI="https://purl.obolibrary.org/obo/uo.obo")
    description = element(root, "fileDescription")
    content = element(description, "fileContent")
    cv(content, "MS:1000579" if ms_level == 1 else "MS:1000580",
       "MS1 spectrum" if ms_level == 1 else "MSn spectrum")
    source_list = element(description, "sourceFileList", count="1")
    source = element(source_list, "sourceFile", id="source1",
                     name="synthetic", location="file:///synthetic")
    cv(source, "MS:1000584", "mzML format")
    samples = element(root, "sampleList", count="1")
    element(samples, "sample", id="sample1", name="synthetic engineering test")
    software_list = element(root, "softwareList", count="1")
    software = element(software_list, "software", id="fixture", version="1")
    cv(software, "MS:1000799", "custom unreleased software tool", "qAlgorithms smoke fixture")
    instruments = element(root, "instrumentConfigurationList", count="1")
    instrument = element(instruments, "instrumentConfiguration", id="IC1")
    components = element(instrument, "componentList", count="3")
    source = element(components, "source", order="1")
    cv(source, "MS:1000073", "electrospray ionization")
    analyzer = element(components, "analyzer", order="2")
    cv(analyzer, "MS:1000484", "orbitrap")
    detector = element(components, "detector", order="3")
    cv(detector, "MS:1000624", "inductive detector")
    processing_list = element(root, "dataProcessingList", count="1")
    processing = element(processing_list, "dataProcessing", id="DP1")
    method = element(processing, "processingMethod", order="0", softwareRef="fixture")
    cv(method, "MS:1000544", "Conversion to mzML")
    run = element(root, "run", id="run1", defaultInstrumentConfigurationRef="IC1",
                  sampleRef="sample1")
    scans = 61 if collective else (8 if blank else 65)
    modes = ("positive", "negative") if polarity == "mixed" else (polarity,)
    spectrum_list = element(run, "spectrumList", count=str(scans * len(modes)),
                            defaultDataProcessingRef="DP1")
    spectrum_ids: list[str] = []
    for scan_number in range(scans):
        for mode in modes:
            rt = 10.0 + 0.5 * scan_number + (0.25 if mode == "negative" else 0.0)
            mz, intensity = profile(mode, rt, scan_number, blank)
            if collective:
                mz, intensity = profile(mode, EXPECTED[mode][0][1], scan_number, False)
                z = (scan_number-30)/12
                factor = max(0.0001, math.exp(-0.5*z**4)*(1+0.1*math.sin(2*math.pi*scan_number/3)))
                # Normalise the analytical mass-Gaussian integral so the EIC
                # has the same 10,000-unit plateau as the selector fixture.
                area_scale = 10000/(150000*0.005*math.sqrt(2*math.pi))
                mz, intensity = mz[:121], [value*factor*area_scale for value in intensity[:121]]
            if intensity_scale != 1.0:
                intensity = [value * intensity_scale for value in intensity]
            spectrum_id = f"scan={len(spectrum_ids) + 1}"
            spectrum = element(spectrum_list, "spectrum", index=str(len(spectrum_ids)),
                               id=spectrum_id, defaultArrayLength=str(len(mz)))
            spectrum_ids.append(spectrum_id)
            cv(spectrum, "MS:1000511", "ms level", str(ms_level))
            cv(spectrum, "MS:1000127" if centroided else "MS:1000128",
               "centroid spectrum" if centroided else "profile spectrum")
            cv(spectrum, "MS:1000130" if mode == "positive" else "MS:1000129",
               f"{mode} scan")
            scan_list = element(spectrum, "scanList", count="1")
            cv(scan_list, "MS:1000795", "no combination")
            scan = element(scan_list, "scan")
            cv(scan, "MS:1000016", "scan start time",
               str(rt / 60.0 if rt_minutes else rt), unitCvRef="UO",
               unitAccession="UO:0000031" if rt_minutes else "UO:0000010",
               unitName="minute" if rt_minutes else "second")
            arrays = element(spectrum, "binaryDataArrayList", count="2")
            array_data = ((mz, mz_bits, mz_compressed, "mz"),
                          (intensity, intensity_bits, intensity_compressed, "intensity"))
            for values, bits, compressed, kind in reversed(array_data) if reverse_arrays else array_data:
                binary_array(arrays, values, bits, compressed, kind)

    xml = ET.tostring(root, encoding="utf-8")
    declaration = b'<?xml version="1.0" encoding="UTF-8"?>\n'
    if not indexed:
        path.write_bytes(declaration + xml)
        return

    # Include real spectrum offsets and a checksum instead of merely wrapping
    # the XML in an indexedmzML tag.
    data = declaration + f'<indexedmzML xmlns="{NS}">\n'.encode("ascii") + xml + b"\n"
    index_offset = len(data)
    entries = []
    previous = 0
    for spectrum_id in spectrum_ids:
        offset = data.index(b"<spectrum ", previous)
        previous = offset + 1
        entries.append(f'<offset idRef="{spectrum_id}">{offset}</offset>')
    index = '<indexList count="1"><index name="spectrum">' + "".join(entries) + '</index></indexList>\n'
    data += index.encode("ascii")
    data += f"<indexListOffset>{index_offset}</indexListOffset>\n<fileChecksum>".encode("ascii")
    data += hashlib.sha1(data).hexdigest().encode("ascii")
    data += b"</fileChecksum>\n</indexedmzML>\n"
    path.write_bytes(data)


def run_command(executable: Path, arguments: list[str], case_dir: Path,
                *, expect_success: bool, log_prefix: str = "") -> subprocess.CompletedProcess[str]:
    command = [str(executable), *arguments]
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=60,
                                 encoding="utf-8", errors="replace", cwd=case_dir)
    except subprocess.TimeoutExpired as error:
        raise TestFailure(f"CLI timed out for {case_dir.name}; possible processing loop") from error
    (case_dir / f"{log_prefix}stdout.log").write_text(process.stdout, encoding="utf-8")
    (case_dir / f"{log_prefix}stderr.log").write_text(process.stderr, encoding="utf-8")
    require(process.returncode >= 0,
            f"{case_dir.name}: process terminated by signal {-process.returncode}\n"
            f"stderr:\n{process.stderr}")
    success = process.returncode == 0
    require(success == expect_success,
            f"{case_dir.name}: exit {process.returncode}; expected "
            f"{'success' if expect_success else 'failure'}\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}")
    return process


def run_program(executable: Path, input_path: Path, output_dir: Path,
                *, expect_success: bool, option: str | None = "--printfeatures",
                extra_arguments: tuple[str, ...] = (),
                log_prefix: str = "") -> subprocess.CompletedProcess[str]:
    output_dir.mkdir(parents=True)
    arguments = ["-i", str(input_path), "-o", str(output_dir)]
    if option is not None:
        arguments.append(option)
    arguments.extend(extra_arguments)
    return run_command(executable, arguments, output_dir.parent,
                       expect_success=expect_success, log_prefix=log_prefix)


def read_features(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), f"Missing feature CSV: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        require(REQUIRED_COLUMNS <= set(reader.fieldnames or []),
                f"Missing required feature columns in {path.name}: {reader.fieldnames}")
        rows = list(reader)
    for row in rows:
        require(None not in row and all(value is not None for value in row.values()),
                f"Truncated/invalid CSV row in {path.name}: {row}")
        for key, value in row.items():
            try:
                number = float(value)
            except ValueError as error:
                raise TestFailure(f"Nonnumeric {key}={value!r} in {path.name}") from error
            require(math.isfinite(number), f"Nonfinite {key}={value} in {path.name}")
        require(float(row["area"]) > 0 and float(row["height"]) > 0,
                f"Nonpositive peak area or height in {path.name}")
    status_columns = {"isUnresolved", "unresolvedReason", "unresolvedMaximaCount", "unresolvedLowerRT", "unresolvedUpperRT"}
    require(status_columns <= set(reader.fieldnames or []), "Unresolved status columns must always be present")
    trace_path = path.with_name(path.name.removesuffix("_features.csv") + "_unresolved_groups.csv")
    require(trace_path.is_file(), f"Missing unresolved-group trace: {trace_path}")
    with trace_path.open(newline="", encoding="utf-8") as stream:
        trace_reader = csv.DictReader(stream)
        require({"featureID", "binID", "maximumNumber", "retentionTime", "reason"} <= set(trace_reader.fieldnames or []), "Invalid unresolved trace schema")
        trace = list(trace_reader)
    known_ids = {row["ID"] for row in rows}
    require(all(item["featureID"] in known_ids for item in trace), "Trace references an unknown feature")
    for row in rows:
        members = [item for item in trace if item["featureID"] == row["ID"]]
        count = int(row["unresolvedMaximaCount"])
        if row["isUnresolved"] == "0":
            require(count == 0 and row["unresolvedReason"] == "0" and not members and
                    float(row["unresolvedLowerRT"]) == 0 and float(row["unresolvedUpperRT"]) == 0,
                    "Resolved features must not acquire collective metadata")
        else:
            require(row["isUnresolved"] == "1" and row["unresolvedReason"] in {"1", "2"} and count >= 2,
                    "Collective status must identify its reason and supported maxima")
            require(len(members) == count, "Collective maxima were lost from its trace")
            lo, hi = float(row["unresolvedLowerRT"]), float(row["unresolvedUpperRT"])
            require(float(row["lowestRetentionTime"])-1e-4 <= lo <= float(row["retentionTime"])+1e-4 and
                    float(row["retentionTime"])-1e-4 <= hi <= float(row["highestRetentionTime"])+1e-4,
                    "Collective region must be covered by its unchanged representative")
            require(all(float(a["retentionTime"])<float(b["retentionTime"]) for a,b in zip(members,members[1:])),
                    "Source maxima must remain unique and ordered")
            require([int(item["maximumNumber"]) for item in members] == list(range(1,count+1)), "Trace maximum numbering is incomplete")
            require(all(item["binID"] == row["binID"] and item["reason"] == row["unresolvedReason"] and
                        float(row["unresolvedLowerRT"])-1e-4 <= float(item["retentionTime"]) <= float(row["unresolvedUpperRT"])+1e-4
                        for item in members), "Trace must preserve origin, reason and region coverage")
    return rows


def check_peaks(path: Path, polarity: str) -> None:
    rows = read_features(path)
    require(bool(rows), f"Feature CSV is empty for known synthetic peaks: {path}")
    require(all(row["isUnresolved"] == "0" for row in rows), "Known resolved synthetic peaks must stay resolved")
    for mass, rt in EXPECTED[polarity]:
        matching = [row for row in rows
                    if abs(float(row["mz"]) - mass) <= 0.002
                    and abs(float(row["retentionTime"]) - rt) <= 0.5]
        require(bool(matching),
                f"Expected peak near m/z={mass}, RT={rt}s not found in {path.name}; "
                f"observed: {[(r['mz'], r['retentionTime']) for r in rows]}")


def run_cases(executable: Path, work_dir: Path, selected: str) -> None:
    for name in PROFILE_CASES:
        if selected != "all" and selected != name:
            continue
        case_dir = work_dir / name
        case_dir.mkdir()
        # A space in both file and directory names exercises real CLI path handling.
        input_path = case_dir / f"{name} profile.mzML"
        output = case_dir / "feature output"
        if name == "malformed":
            input_path.write_text("<mzML><run>", encoding="utf-8")
        else:
            write_fixture(
                input_path,
                polarity=name if name in ("positive", "negative", "mixed") else "positive",
                indexed=name in ("negative", "mixed"),
                mz_compressed=name in ("negative", "mixed"),
                intensity_compressed=name == "negative",
                reverse_arrays=name == "negative",
                rt_minutes=name == "negative",
                blank=name == "blank", centroided=name == "centroided",
                ms_level=2 if name == "ms2_only" else 1,
            )
        success = name not in ("malformed", "centroided", "ms2_only")
        run_program(executable, input_path, output, expect_success=success,
                    option="-printfeatures" if name == "positive" else "--printfeatures")
        if name == "mixed":
            for mode in ("positive", "negative"):
                check_peaks(output / f"{input_path.stem}_{mode}_features.csv", mode)
            require(len(list(output.glob("*_features.csv"))) == 2,
                    "Mixed-polarity input should produce exactly two feature lists")
        elif name in ("positive", "negative"):
            check_peaks(output / f"{input_path.stem}_features.csv", name)
        elif name == "blank":
            require(read_features(output / f"{input_path.stem}_features.csv") == [],
                    "Blank profile input should produce a header-only feature CSV")
        else:
            require(not list(output.glob("*_features.csv")),
                    f"Rejected input {name} produced a misleading feature CSV")
        print(f"PASS {name}", flush=True)

    for name in CLI_CASES:
        if selected != "all" and selected != name:
            continue
        case_dir = work_dir / name
        case_dir.mkdir()
        if name == "collective_output":
            input_path = case_dir / "unresolved plateau.mzML"
            write_fixture(input_path, collective=True)
            output = case_dir / "feature output"
            process = run_program(executable, input_path, output, expect_success=True)
            rows = read_features(output / "unresolved plateau_features.csv")
            require(len(rows)==1 and rows[0]["isUnresolved"]=="1", "The ambiguous flat-top fixture must produce one marked collective")
            require(abs(float(rows[0]["mz"])-EXPECTED["positive"][0][0])<0.002, "Collective must preserve its mass bin")
            require("explicitly marked unresolved" in process.stderr, "A successful collective output must remain visible in the CLI")
        elif name in {"unresolved_trace_write_failure", "orphan_unresolved_trace"}:
            input_path = case_dir / "sample.mzML"
            write_fixture(input_path)
            output = case_dir / "feature output"
            output.mkdir()
            trace = output / "sample_unresolved_groups.csv"
            if name == "unresolved_trace_write_failure":
                trace.mkdir()
                process = run_command(executable, ["-i", str(input_path), "-o", str(output), "--printfeatures"], case_dir, expect_success=False)
                require("unresolved-group trace" in process.stderr, "A missing companion must fail visibly")
                require(trace.is_dir(), "The blocking directory must not be replaced")
            else:
                trace.write_text("preserve this prior trace\n", encoding="utf-8")
                run_command(executable, ["-i", str(input_path), "-o", str(output), "--printfeatures", "--skip-existing"], case_dir, expect_success=False)
                require(trace.read_text(encoding="utf-8") == "preserve this prior trace\n", "An orphan trace must not be overwritten")
                require(not (output / "sample_features.csv").exists(), "A new feature CSV must not be paired with an old orphan trace")
        elif name == "orphan_feature_output":
            input_path = case_dir / "sample.mzML"
            write_fixture(input_path)
            output = case_dir / "feature output"
            run_program(executable, input_path, output, expect_success=True)
            primary = output / "sample_features.csv"
            trace = output / "sample_unresolved_groups.csv"
            before = primary.read_bytes(), trace.read_bytes()
            run_command(executable, ["-i", str(input_path), "-o", str(output), "--skip-existing"], case_dir, expect_success=True)
            require(before == (primary.read_bytes(), trace.read_bytes()), "An existing output pair must be preserved")
            trace.unlink()
            process = run_command(executable, ["-i", str(input_path), "-o", str(output), "--skip-existing"], case_dir, expect_success=False)
            require("unresolved-group trace pair" in process.stderr, "An orphan primary must fail visibly")
            require(primary.read_bytes() == before[0] and not trace.exists(), "An orphan primary must not be overwritten or silently repaired")
        elif name == "missing_output_value":
            input_path = case_dir / "valid.mzML"
            write_fixture(input_path)
            run_command(executable, ["-i", str(input_path), "-o", "-printfeatures"],
                        case_dir, expect_success=False)
            require(not (case_dir / "-printfeatures").exists(),
                    "The switch following -o was incorrectly used as an output directory")
        elif name == "output_collision":
            input_dir = case_dir / "inputs"
            input_dir.mkdir()
            pure = input_dir / "sample_positive.mzML"
            mixed = input_dir / "sample.mzML"
            # The smaller pure file is processed first. Different intensities
            # make an accidental overwrite observable even when both inputs
            # contain peaks at the same m/z and retention times.
            write_fixture(pure, intensity_scale=1.5)
            write_fixture(mixed, polarity="mixed")
            reference = case_dir / "reference output"
            run_program(executable, pure, reference, expect_success=True,
                        log_prefix="reference-")
            expected_csv = reference / "sample_positive_features.csv"
            check_peaks(expected_csv, "positive")
            output = case_dir / "batch output"
            run_program(executable, input_dir, output, expect_success=False)
            actual_csv = output / "sample_positive_features.csv"
            check_peaks(actual_csv, "positive")
            require(actual_csv.read_bytes() == expected_csv.read_bytes(),
                    "The mixed-polarity input overwrote the first sample's feature CSV")
        elif name == "skip_error":
            input_dir = case_dir / "inputs"
            input_dir.mkdir()
            (input_dir / "malformed.mzML").write_text("<mzML><run>", encoding="utf-8")
            valid = input_dir / "valid.mzML"
            write_fixture(valid)
            output = case_dir / "batch output"
            run_program(executable, input_dir, output, expect_success=False,
                        extra_arguments=("--skip-error",))
            check_peaks(output / "valid_features.csv", "positive")
            require(not (output / "malformed_features.csv").exists(),
                    "Malformed input produced a misleading feature CSV during --skip-error")
        elif name == "default_features":
            input_path = case_dir / "default output.mzML"
            write_fixture(input_path)
            output = case_dir / "feature output"
            run_program(executable, input_path, output, expect_success=True, option=None)
            check_peaks(output / "default output_features.csv", "positive")
        print(f"PASS {name}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--work-dir", type=Path, help="retain generated files beneath this directory")
    parser.add_argument("--case", default="all",
                        choices=("all", *PROFILE_CASES, *CLI_CASES))
    args = parser.parse_args()
    executable = args.executable.resolve()
    if not executable.is_file():
        parser.error(f"Executable does not exist: {executable}")
    if args.work_dir:
        args.work_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix="mzml-cli-", dir=args.work_dir)).resolve()
    try:
        run_cases(executable, work_dir, args.case)
    except (TestFailure, OSError) as error:
        print(f"FAIL: {error}\nGenerated files and logs: {work_dir}", file=sys.stderr)
        return 1
    if args.work_dir:
        print(f"Generated files and logs: {work_dir}")
    else:
        shutil.rmtree(work_dir)
    print("Synthetic mzML CLI smoke tests passed (not analytical validation).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
