"""Fit unequal-length scans on one fixed grid and verify FIFO membership."""

import csv
import math
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET

from test_mzml_cli import NS, binary_array, element, run_command, run_program, write_fixture


def write_profiles(path):
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
        settings = ("--superspectrum", ".01", ".025", "11")
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
                                   "--superspectrum", "1e-20", ".025", "11"], root,
                             expect_success=False)
        assert "step is too small" in failed.stderr
        assert path.read_bytes() == previous and list(output.iterdir()) == [path]


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
