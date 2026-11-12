"""Build a model after the host exists, then load it into that same host."""

import hashlib
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    host_arg, cmake, generator, compiler, config, platform, toolset, ownership_test = sys.argv[1:]
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
        validation_only = Path((build / f"validation_only_gaussian-{config}.txt").read_text())

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

        def search(values, *, mixed=False, rotations=399, alpha=".05", minimum_area="100", axis=None, module=valid):
            builtin = ["--quadratic", "100", ".7"] if mixed else []
            axis_args = ["--axis", *map(str, axis)] if axis else []
            return subprocess.run([str(host), "--peak-model-module-search", str(module), "1", minimum_area, ".7",
                                   str(alpha), str(rotations), "7", "8", *builtin, *axis_args, *values],
                                  capture_output=True, text=True)

        def family(process):
            assert process.returncode == 0, (process.stdout, process.stderr)
            summary, windows = {}, []
            for line in process.stdout.splitlines():
                key, value = line.split("=", 1)
                if key == "window":
                    name, center, scale, status, background, adjusted = value.split(",")
                    windows.append((name, int(center), int(scale), status, background, float(adjusted)))
                else:
                    summary[key] = float(value)
            assert len(windows) == summary["familySize"], (summary, windows)
            assert sum(window[3] == "passed" for window in windows) == summary["candidates"], (summary, windows)
            return summary, windows

        native_summary, native_windows = family(search(samples))
        mixed_summary, mixed_windows = family(search(samples, mixed=True))
        assert native_summary["familySize"] == 4 and native_summary["windowAlpha"] == .05 / 4, native_summary
        assert mixed_summary["familySize"] == 8 and mixed_summary["windowAlpha"] == .05 / 8, mixed_summary
        assert {row[0] for row in mixed_windows if row[3] == "passed"} == {
            "builtin:piecewise-quadratic", "native:external-gaussian-ols"}, mixed_windows
        assert [row[1:3] for row in native_windows] == [(7, 7), (8, 7), (8, 8), (9, 7)], native_windows
        for solo, shared in zip(native_windows, mixed_windows[4:]):
            assert solo[:3] == shared[:3], (solo, shared)
            assert math.isclose(shared[-1], min(1, 2 * solo[-1])), (solo, shared)

        skipped_summary, skipped_windows = family(search(["0", *samples[1:]], mixed=True))
        assert skipped_summary["familySize"] == 8 and sum(row[3] == "unavailable" for row in skipped_windows) == 4
        assert all(math.isnan(row[-1]) for row in skipped_windows if row[3] == "unavailable"), skipped_windows
        budget = search(samples, mixed=True, rotations=99)
        assert budget.returncode != 0 and "Rotation budget cannot resolve" in budget.stderr, budget

        def measurement(process):
            assert process.returncode == 0, (process.stdout, process.stderr)
            peaks = [list(map(float, row[5:].split(","))) for row in process.stdout.splitlines() if row.startswith("peak=")]
            assert len(peaks) == 1, process.stdout
            return peaks[0]

        base_peak = measurement(search(samples, axis=(0, 1)))
        reversed_peak = measurement(search(samples, axis=(150, -.25)))
        assert reversed_peak[:3] == base_peak[:3] and base_peak[:2] == [8, 8]
        assert math.isclose(reversed_peak[3], 150 - .25 * base_peak[3])
        assert math.isclose(reversed_peak[4], .25 * base_peak[4])
        assert base_peak[5] > 0 and base_peak[8] > 0 and base_peak[6] != 0, base_peak
        for i, factor in enumerate((.0625, -.0625, -.0625, .0625)):
            assert math.isclose(reversed_peak[5 + i], factor * base_peak[5 + i]), reversed_peak
        old_abi = search(samples, axis=(0, 1), module=validation_only)
        assert old_abi.returncode == 0 and "unmeasured=8,8," in old_abi.stdout and "peak=" not in old_abi.stdout
        subprocess.run([ownership_test, str(valid)], check=True)

    assert hashlib.sha256(host.read_bytes()).hexdigest() == original_hash
    assert host.stat().st_mtime_ns == original_mtime
    print(f"Separately built external Gaussian passed all stages on unchanged host SHA-256 {original_hash}")


if __name__ == "__main__":
    main()
