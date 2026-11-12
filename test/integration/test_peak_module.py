"""Build a model after the host exists, then load it into that same host."""

import hashlib
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    host_arg, cmake, generator, compiler, config, platform, toolset = sys.argv[1:]
    host = Path(host_arg).resolve()
    original_hash = hashlib.sha256(host.read_bytes()).hexdigest()
    original_mtime = host.stat().st_mtime_ns
    test_dir = Path(__file__).resolve().parent
    samples = [str(800 * math.exp(-.055 * x * x + .003 * math.sin(2 * x))) for x in range(-8, 9)]

    with tempfile.TemporaryDirectory(prefix="qAlgorithms external model ") as temporary:
        build = Path(temporary)
        configure = [cmake, "-S", str(test_dir / "peak_module_fixture"), "-B", str(build), "-G", generator,
                     f"-DCMAKE_CXX_COMPILER={compiler}", f"-DCMAKE_BUILD_TYPE={config}",
                     f"-DQALGORITHMS_ABI_INCLUDE={test_dir.parents[1] / 'include'}"]
        if platform:
            configure.extend(["-A", platform])
        if toolset:
            configure.extend(["-T", toolset])
        subprocess.run(configure, check=True)
        build_command = [cmake, "--build", str(build)]
        if config:
            build_command.extend(["--config", config])
        subprocess.run(build_command, check=True)
        valid = Path((build / f"external_gaussian-{config}.txt").read_text())
        incompatible = Path((build / f"incompatible_gaussian-{config}.txt").read_text())

        def run(module, minimum="100", values=samples):
            return subprocess.run([str(host), "--peak-model-module", str(module), "1", minimum, ".7", *values],
                                  capture_output=True, text=True)

        def stages(process):
            assert process.returncode == 0, (process.stdout, process.stderr)
            return dict(line.split("=", 1) for line in process.stdout.splitlines())

        passed = stages(run(valid))
        assert passed.pop("model") == "external-gaussian-ols", passed
        assert passed == dict.fromkeys(("status", "form", "apex", "formConfidence", "area", "coverage", "background"),
                                       "passed"), passed
        failed = stages(run(valid, "1e9"))
        assert failed["status"] == failed["area"] == "failed", failed
        assert failed["coverage"] == failed["background"] == "notEvaluated", failed
        invalid = run(incompatible)
        assert invalid.returncode != 0 and "Invalid v1 peak model module descriptor" in invalid.stderr, invalid
        missing = run(build / "missing-module")
        assert missing.returncode != 0 and "Cannot load the peak model module" in missing.stderr, missing
        zero = stages(run(valid, values=["0", *samples[1:]]))
        assert zero["status"] == "unavailable" and zero["background"] == "notEvaluated", zero
        failed_fit = stages(run(valid, values=[".5", *samples[1:]]))
        assert failed_fit["status"] == "failed" and failed_fit["apex"] == "notEvaluated", failed_fit

    assert hashlib.sha256(host.read_bytes()).hexdigest() == original_hash
    assert host.stat().st_mtime_ns == original_mtime
    print(f"Separately built external Gaussian passed all stages on unchanged host SHA-256 {original_hash}")


if __name__ == "__main__":
    main()
