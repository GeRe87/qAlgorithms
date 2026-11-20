# qAlgorithms mzML benchmark

Build the native profiler with the project's normal CMake toolchain. Python
3.10 or newer runs the benchmark. The desktop interface also needs Tkinter.

```sh
cmake -S . -B build
cmake --build build --target qmzml-profile qAlgorithms
python tools/mzml_benchmark/benchmark.py gui
```

In the desktop window choose the profile mzML files and an empty output
directory. Tick the replicate box when the input files measure the same sample.
On Windows, `benchmark.pyw` in this directory opens the same window by double-click
when Windows associates .pyw files with Python.
If the profiler lives outside the default build directory, select its executable
once in the profiler field for the current session. Select **Kalibrieren und
erzeugen** to create the benchmark. The second tab scores a feature list.

The equivalent command uses one input by default. Pass `--replicate-inputs`
with two or more files from the same sample, and choose `--replicates` for the
number of simulated runs. A separate background-only experiment uses
`--peaks 0`. Optional `--dropout 0.2` removes each known feature independently
with probability 0.2 from runs after the first.

```sh
python tools/mzml_benchmark/benchmark.py build measured.mzML -o benchmark --seed 53
build/qAlgorithms --benchmark-mode -i benchmark/synthetic/synthetic-001.mzML -o detected
python tools/mzml_benchmark/benchmark.py score benchmark/synthetic/synthetic-001-truth.json detected/synthetic-001_features.csv -o score.json
```

Use `--profiler PATH` or set `QALGORITHMS_PROFILER` for a custom build
directory. Use the individual `calibrate` and `simulate` commands to reuse
one calibration with several seeds. Both support `--help`.

The profiler and comparison scripts use qAlgorithms benchmark mode.
For direct CLI use, pass `--benchmark-mode`. It permits equal adjacent
float coordinates during hyperbolic fitting and reports unsupported
apex groups omitted during centroiding. If an EIC has no valid collective
description, it retries using the existing policy that retains supported
independent peaks and reports omitted components. If recovery fails,
the file fails explicitly. The default qAlgorithms CLI uses its usual
conflict policy.

A clean-checkout smoke test creates its own small profile mzML fixture,
builds a calibration, generates mzML and truth, and scores the outputs of
both kernels.

```sh
ctest --test-dir build -R 'test_mzml_benchmark|test_mzml_cli|test_single_trace_kernel_cli' --output-on-failure
```

Python callers can import the functions directly.

```python
from benchmark import build, score

build(["measured.mzML"], "benchmark", profiler="build/qmzml-profile", seed=53)
result = score("benchmark/synthetic/synthetic-001-truth.json", "detected_features.csv")
```

Set the Python import path to `tools/mzml_benchmark` when calling from another
directory. In C++ and R the thin interfaces invoke the same CLI and write the
same output files.

```cpp
#include "benchmark.hpp"

int main()
{
    return qAlgorithms::benchmark::run("python", "tools/mzml_benchmark/benchmark.py",
        {"simulate", "benchmark/calibration", "-o", "synthetic", "--seed", "53"});
}
```

```r
source("tools/mzml_benchmark/benchmark.R")
mzml_benchmark(c("simulate", "benchmark/calibration", "-o", "synthetic", "--seed", "53"),
               script = "tools/mzml_benchmark/benchmark.py")
```

The calibration directory contains input-level counts and, for each polarity,
original scan and point intervals, centroid properties and uncertainties, bin
membership and pointwise prediction intervals, resolved feature membership and
the standard feature CSV including unresolved descriptions. Point intervals are
zero-based and end-exclusive. A point covered by several fits counts once in
its deepest class. With replicate matching a feature is reproducible when it
has a one-to-one partner in at least one other input file of the same polarity
within the configured m/z and RT tolerances. Otherwise its reproducibility is
unknown when the inputs have not been declared replicates.

The synthetic directory contains the profile mzML files, the manifest with the
seed, calibration hash and version, per-run feature truth, and gzip-compressed
point truth. Each truth record locates a scan segment and records its noise-free
signal and effective additive noise. The manifest supplies the baseline.
Synthetic areas integrate the noise-free signal on the actual synthetic m/z
and RT grids. The scorer accepts qAlgorithms feature CSVs or tables containing
`mz`, `rt_s` (seconds), optional `id` and optional `area`. It reports
TP, FP and FN using one-to-one matching with configurable tolerances.

This implementation samples isolated Gaussian and hyperbolic-secant profiles
using joint templates from detected real features, local m/z spacing, typical RT
sampling and observed half-width asymmetry where available. It uses empirical
outside-centroid intensities for the background and local log differences near
centroid apices as a proxy for relative point noise. Outside-centroid
measurements may include unresolved signals. The background is sampled
independently across points. Features retain the same m/z and RT in simulated
replicates. The optional dropout and independent noise change their observations.
This first version covers isolated MS1 peaks. Overlapping ions, correlated
background, acquisition-specific drift and a calibrated false alarm rate
require additional models and validation.

To compare kernels on a common false alarm budget, use:

```sh
python tools/mzml_benchmark/compare_kernels.py benchmark/calibration --qalgo build/qAlgorithms -o comparison
```

Each seed creates a signal run and a background-only run on the same scan times
and m/z grid. The first 12 seeds calibrate a minimum DQSF threshold separately
for each kernel with a shared maximum of 0.5 resolved false features per null
run. The next 24 seeds evaluate false alarms and recovery with the threshold
held fixed. If a kernel already meets the budget, it keeps its original output.
Each signal run contains a balanced mix of the four relative amplitude factors
0.03, 0.1, 0.3 and 1. The factors scale the sampled template height and are
not signal-to-noise ratios. The report includes matched measurement errors
and a bootstrap interval that resamples paired runs by seed.

Measure density and overlap in an existing real-data calibration with:

```sh
python tools/mzml_benchmark/measure_density.py benchmark/calibration -o density
```

The tile tables count centroid apices and feature apices on a 20 m/z by 30 s
grid and include the number of selected MS1 scans per time tile. The summary
also reports neighbouring centroid distances, intersections of the original
profile fit spans, and intersections of centroid m/z FWHM intervals. At the
feature level it counts RT fit-window intersections within each bin and
simultaneous RT and m/z intersections between bins. The latter uses each
bin's median centroid FWHM as an approximate mass width. Unresolved
descriptions and their supported maxima are counted separately. These are
measurements of qAlgorithms output on that input. The similarity of these
measurements under different kernels has to be evaluated separately.

The GUI and `build` use measured density and overlap automatically:

```sh
python tools/mzml_benchmark/benchmark.py build measured.mzML -o benchmark
python tools/mzml_benchmark/benchmark.py simulate benchmark/calibration --density benchmark/density/density.json -o second-run
```

The density step links neighbouring resolved features into disjoint pairs.
A same-bin pair must have an RT separation below one mean estimated RT FWHM.
A cross-bin pair must intersect in RT and be within one mean FWHM in each
dimension. The nearest eligible pairs are selected first. Regional simulation
samples the measured mixture of isolated groups, m/z pairs and RT pairs, and
places them in the corresponding observed regions. Partner separations are
sampled from the measured histograms in units of mean FWHM. The RT
normalization uses the reported fit-window duration divided by three as a width
proxy. Larger clusters and unresolved descriptions do not enter the estimate.

Use `--overlap .5` to force half the features into pairs for a controlled run,
or `--overlap 0` for regional sampling with isolated peaks. Use
`--overlap-mode mz` or `rt` to focus on one type. If a forced pair type has
no measured neighbours, its fallback separation is 0.7 mean FWHM.
Use `build --no-regional` for the original isolated benchmark with no
density calibration.

Each pair is rendered as two additive signals on a shared profile grid. The
manifest stores pair members, source tiles and realized separations. The point
truth stores each member's contribution and the noise added to the combined
signal. Regional mode with zero overlap samples density with isolated peaks.

`compare_overlap.py` compares both kernels on a paired experiment:

```sh
python tools/mzml_benchmark/compare_overlap.py benchmark/calibration benchmark/density/density.json --qalgo build/qAlgorithms -o overlap-comparison
```

For each seed the generator renders the original placement and a control
with the second member of each pair on its own mass grid. Signal shape,
retention time, amplitude, sampled truth area and noise draws are shared.
The control uses extra m/z grid points for the relocated lanes. Its raw
false-positive counts therefore have a different exposure. The report gives
the number of pairs with zero, one and two recovered features and compares
the area of single detections to each matched member and to the pair total.

Repeat the experiment with RT widths estimated from the measured points by
running `measure_density.py --rt-width-method observed` into a new density
directory. This width uses the area trace, finds the observed point nearest
the fitted apex, then interpolates the first half-height crossing on each
side within the fit window. Features missing a crossing use the previous
window-based width. The density summary counts the directly estimated widths.
