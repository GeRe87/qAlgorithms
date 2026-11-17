"""Exercise both kernels through the regular mzML workflow."""
import csv
import math
from pathlib import Path
import sys
import tempfile

from test_mzml_cli import check_peaks, run_program, write_fixture


def main(host):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "profile.mzML"
        write_fixture(source)
        for mode, flags in (("default", ()), ("quadratic", ("--peak-model", "quadratic")),
                            ("hyperbolic", ("--peak-model", "hyperbolic"))):
            run_program(host, source, root / mode, expect_success=True,
                        extra_arguments=("--printcentroids", "--printbins", *flags))
        for suffix in ("centroids", "bins", "features", "unresolved_groups"):
            name = f"profile_{suffix}.csv"
            assert (root / "default" / name).read_bytes() == (root / "quadratic" / name).read_bytes()
        features = root / "hyperbolic" / "profile_hyperbolic_features.csv"
        check_peaks(features, "positive")
        with features.open(newline="") as stream:
            for row in csv.DictReader(stream):
                assert float(row["b2"]) > abs(float(row["b3"]))
                for field in ("area", "height", "areaUncertainty", "heightUncertainty", "RT_Uncertainty"):
                    assert math.isfinite(float(row[field])) and float(row[field]) >= 0
        run_program(host, source, root / "invalid", expect_success=False,
                    extra_arguments=("--peak-model", "unknown"))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
