"""Exercise shared profile masses through the normal mzML feature pipeline."""

import csv
import math
from pathlib import Path
import sys
import tempfile

from test_mzml_cli import EXPECTED, run_program, write_fixture


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main(executable):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "profile.mzML"
        write_fixture(source)
        run_program(executable, source, root / "default", expect_success=True)
        run_program(executable, source, root / "fifo", expect_success=True, option=None,
                    extra_arguments=("--fifo-feature-mz",))
        original = read_rows(root / "default" / "profile_features.csv")
        measured = read_rows(root / "fifo" / "profile_features.csv")
        assert len(original) == len(measured) == len(EXPECTED["positive"])
        for old, new, (expected_mz, _) in zip(original, measured, EXPECTED["positive"]):
            assert all(new[key] == value for key, value in old.items())
            assert abs(float(new["mzJoint"]) - expected_mz) < 0.0001
            assert math.isfinite(float(new["mzJointStdError"])) and float(new["mzJointStdError"]) > 0
            assert int(new["mzJointScans"]) > 5
            assert int(new["mzJointPoints"]) >= 3 * int(new["mzJointScans"])
        run_program(executable, source, root / "incompatible", expect_success=False, option=None,
                    extra_arguments=("--fifo-feature-mz", "--printcentroids"))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
