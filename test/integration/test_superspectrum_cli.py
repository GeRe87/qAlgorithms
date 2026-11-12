"""Fit unequal-length scans on one fixed grid and verify FIFO membership."""

import csv
import math
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET

from test_mzml_cli import NS, binary_array, element, run_command, run_program, write_fixture


def write_profiles(path, boundary=False):
    write_fixture(path, blank=True)
    tree = ET.parse(path)
    spectra = tree.find(f".//{{{NS}}}spectrumList")
    spectra.set("count", "6")
    for spectrum in list(spectra)[6:]:
        spectra.remove(spectrum)
    for scan, (spectrum, count) in enumerate(zip(spectra, (12, 13, 14, 0, 0, 15))):
        mz = [100.0 + j * .24 / (count - 1) for j in range(count)]
        if scan == 3:
            # A trace starting at the apex must not split the full traces.
            mz = [100.125 + j * .015 for j in range(8)]
        raw = [(scan + 1) * 100 * math.exp(-.5 * ((mass - 100.125) / .038)**2)
               for mass in mz]
        # A distant flat profile cannot bridge the gap or supply missing coverage.
        mz += [101 + j * .02 for j in range(12)]
        raw += [4] * 12
        if boundary:
            lower = 111.92 if scan == 0 else 111.964
            high_points = 21 if scan == 5 else 7
            mz = [lower + j * .001 for j in range(round((112 - lower) / .001))]
            mz += [112 + j * .004 for j in range(high_points)]
            raw = [max(.001, (scan + 1) * 100 * math.exp(-.5 * ((mass - 111.99) / .005)**2))
                   for mass in mz]
        spectrum.set("defaultArrayLength", str(len(mz)))
        spectrum.remove(spectrum.find(f"{{{NS}}}binaryDataArrayList"))
        arrays = element(spectrum, "binaryDataArrayList", count="2")
        binary_array(arrays, mz, 64, False, "mz")
        binary_array(arrays, raw, 32, True, "intensity")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def main(host):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "profile.mzML"
        write_profiles(source)
        settings = ("--superspectrum-fixed", ".01", ".025", "11")
        output = root / "out"
        run_program(host, source, output, expect_success=True, option=None,
                    extra_arguments=settings)
        path = output / "profile_superspectrum_centroids.csv"
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 6, rows
        for start, members in ((0, {0, 1, 2}), (1, {1, 2, 5})):
            group = [row for row in rows if int(row["firstScan"]) == start]
            assert {int(row["scan"]) for row in group} == members, group
            assert {int(row["lastScan"]) for row in group} == {start + 4}
            for key in ("centroid", "mz", "gridCenterMz", "b1", "b2", "b3"):
                assert len({row[key] for row in group}) == 1
            for row in group:
                assert abs(float(row["mz"]) - 100.125) < .004, row
                assert math.isclose(float(row["height"]), (int(row["scan"]) + 1) * 100, rel_tol=.10)
                grid_center = float(row["gridCenterMz"]) / .01
                assert abs(grid_center - round(grid_center)) < 1e-8
        previous = path.read_bytes()
        failed = run_command(host, ["-i", str(source), "-o", str(output),
                                   "--superspectrum-fixed", "1e-20", ".025", "11"], root,
                             expect_success=False)
        assert "step is too small" in failed.stderr
        assert path.read_bytes() == previous and list(output.iterdir()) == [path]
        auto_output = root / "auto"
        result = run_program(host, source, auto_output, expect_success=True, option=None,
                             extra_arguments=("--superspectrum",))
        assert "Superspectrum pooled settings:" in result.stdout
        with (auto_output / path.name).open(newline="") as stream:
            automatic = list(csv.DictReader(stream))
        assert {(row["firstScan"], row["scan"]) for row in automatic} == {
            (row["firstScan"], row["scan"]) for row in rows}
        assert len(automatic) == len(rows)
        assert all(abs(float(row["mz"]) - 100.125) < .004 for row in automatic)
        steps = {float(row["gridStepMz"]) for row in automatic}
        assert len(steps) == 1 and .007 < next(iter(steps)) < .012

        boundary = root / "boundary.mzML"
        write_profiles(boundary, boundary=True)
        run_program(host, boundary, root / "boundary", expect_success=True, option=None,
                    extra_arguments=("--superspectrum",))
        with (root / "boundary" / "boundary_superspectrum_centroids.csv").open(newline="") as stream:
            near_boundary = [row for row in csv.DictReader(stream) if abs(float(row["mz"]) - 111.99) < .004]
        assert {row["firstScan"] for row in near_boundary} == {"0", "1"}
        assert len({row["gridStepMz"] for row in near_boundary}) == 1
        assert all(math.isclose(float(row["gridStepMz"]), .0005, rel_tol=1e-8) for row in near_boundary)


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
