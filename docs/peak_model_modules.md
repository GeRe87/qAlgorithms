# Native peak-model modules

The host loads a separately built shared library from an explicit path.
The library exports the C symbol `qalgorithms_peak_model_v1` and returns
the descriptor defined in `include/qalgorithms_peak_model_abi.h`.
The host checks the ABI version and descriptor size, then dispatches through
its callbacks. The descriptor supplies the model name.

## Source layout

- `qalgorithms_peak_module.cpp` owns the library handle and adapts the native fit to the qPeaks validation code.
- `detail/qalgorithms_peak_cascade.h` orders the shared validation stages.
- `qalgorithms_peak_background.cpp` compares raw backgrounds and replays the same model fit on rotated observations.
- `qalgorithms_peak_builtin_models.cpp` connects the enum-based models to the validation code.

The library owns its opaque fit and supplies its destructor. The host
destroys each fit before starting another one in the same validation.
Callbacks borrow input and output buffers for the duration of the call.
Concurrent calls require support from the module. Use modules you trust
to execute as native code in the qAlgorithms process.

## Run a fixed-window diagnostic

```text
qAlgorithms --peak-model-module PATH STEP MIN_AREA MIN_COVERAGE I0 I1 I2 I3 I4 [...]
```

Supply an odd-length prespecified window with at least five intensities.
The window center is coordinate zero. The module receives the declared
sample step and reports area in intensity times that coordinate's unit.
The command prints the model name and each stage's status. The C++ entry
point is `NativePeakModel::validateCandidate` in `include/qalgorithms_peak_module.h`.

The module defines its shape, fit and statistical callbacks. The host
checks the structure and the plausibility of a passing result. Validation
of the statistical assumptions belongs to that model. The external entry
point currently covers a single prespecified window.

The `test_peak_module` integration test builds the fixture in a separate
CMake project after the host exists. It loads the resulting library into
that executable, checks the stages and rejects an incompatible descriptor.
The test compares the host's hash and modification time before and after.
The fixture demonstrates the ABI and uses test-specific margins.

## Run a declared window family

```text
qAlgorithms --peak-model-module-search PATH STEP MIN_AREA MIN_COVERAGE ALPHA ROTATIONS MIN_SCALE MAX_SCALE [--quadratic MIN_AREA MIN_COVERAGE] I0 I1 I2 I3 I4 [...]
```

The command enumerates every fitting center and half-width in the supplied
range. Add `--quadratic` to test the same windows with the built-in
quadratic model as well. Its area threshold is in intensity times sample
index. The module threshold uses intensity times the declared coordinate unit.

In C++, pass `NativeCandidateModelWindows` entries to `CandidateSearchPlan`.
The plan owns shared references to the loaded models. The C++ API also
accepts both built-in models and explicit, model-specific window lists.
Each native descriptor name must be unique within the family. Results use
the descriptor name for native windows and an optional model enum for built-in
windows. The command prefixes window labels with `native:` or `builtin:`
and percent-encodes delimiters in native names.

Declare the complete model and window family, thresholds and random seeds
before inspecting the intensities. The raw-background stage uses the
same Bonferroni allocation for all model-window pairs, including skipped
windows. The rotation budget must resolve that threshold. The guarantee
uses the fixed-window raw Gaussian null and the rotation assumptions.
Module-specific statistical stages keep their own assumptions.

The integration test uses the separately built module in both a native-only
family and a mixed family. It checks that adding the built-in model doubles
the denominator and that windows containing zero remain counted.

## Select and measure accepted peaks

Add `--axis ORIGIN UNITS_PER_SAMPLE` after the optional quadratic settings
and before the intensities. The origin addresses sample zero. Supply a
signed affine scale valid over the full modeled peak.

The search assigns accepted windows to observed maxima using their fitted
apices. Windows assigned to the same maximum are compared on their common
sample interval using raw-intensity squared error from the original log-OLS
curves. Exact score ties and insufficient common support remain unresolved.
The background test continues to use its own raw amplitude and background fit.

A measured record includes the observed maximum, selected window index,
position, modeled total area and the four entries of their covariance.
Area uses intensity times the measurement coordinate. The covariance is
conditional on a fixed window, model and affine calibration. It excludes
cross-peak, calibration and selection uncertainty. The modeled area uses
the original log fit, with the legacy continuous valley truncation for
mixed quadratic flanks.

A native module may export `qalgorithms_peak_model_properties_v1` with
the signature in the ABI header. It supplies position, total area and
their joint covariance from the same log-OLS fit. A validation-only module
still participates in selection. A selected window without valid joint
properties is reported as `unmeasured`.

## Use the regular mzML workflow

```text
qAlgorithms -i profile.mzML -o output --calibrated-peaks MODEL STEP MIN_AREA_MZ MIN_COVERAGE MIN_SCALE MAX_SCALE ALPHA ROTATIONS
```

`MODEL` is `quadratic`, `hyperbolic` or `module:PATH`. Quadratic requires
`STEP=1` and a maximum scale within the compiled coefficient table.
The other models use their declared fit-coordinate step. Native windows
must have more observations than fitted parameters.
`MIN_AREA_MZ` is in intensity times m/z and is converted to fit units
before the area test.

The mode uses the normal mzML reader, MS1 profile selection and polarity
handling. It writes `<input>_calibrated_peaks.csv` per polarity. The CSV
contains one summary per scan, measured representatives and explicit
unresolved or unassigned records. Position is in m/z, area is in intensity
times m/z, and the covariance columns describe those two quantities.
This mode ends after scan-level peak output. The existing feature workflow
remains the default when the option is omitted.

Each scan declares its family from sample count and the scale range before
intensities are tested. The Bonferroni limit applies within that scan under
the raw Gaussian and rotation assumptions. A window with zero or other
nonpositive intensities remains in the denominator. A nonaffine m/z grid
is reported as `unsupportedAxis`. The grid check allows float storage
rounding of one ULP relative to the affine line through the scan endpoints.
The reported covariance is conditional on that calibration and the chosen
fit. The rotation budget must resolve the largest tested scan family.

## Joint m/z for chromatographic features

Use `qAlgorithms -i profile.mzML -o output --fifo-feature-mz` to add
`mzJoint`, `mzJointStdError`, `mzJointScans` and `mzJointPoints` to the
regular feature CSV. The ordinary qPeaks workflow first detects centroids,
builds EICs and identifies resolved features. The joint fit then rereads
only the spectra needed by those features. The existing feature columns
retain their original definitions, including the area along retention time.
Use a fresh output directory when adding this option with `--skip-existing`.

Each feature uses a fixed common m/z interval, the intersection of the
centroid m/z ± one FWHM intervals in its measured scans. The fit uses each
scan's recorded profile points at their original precision, including
nonuniform and scan-specific m/z grids. Scans with fewer than three usable
points, repeated coordinates in the segment or nonpositive intensity in
the segment are omitted. Interpolated EIC entries do not contribute.

The five-scan FIFO adds and subtracts centered regression cross-products.
The closing fit combines the expired and current contributions, using each
scan once. The local model uses a shared log-quadratic shape with a separate
log amplitude per scan. The reported m/z is its common apex. The standard
error is conditional on independent, equal-variance log residuals, the
selected scans and the fixed segment. It does not include m/z calibration,
drift or feature selection uncertainty. The joint columns are empty when
there is insufficient support for a concave fit with an apex covered by
all contributing scans, or when the chromatographic feature is unresolved.

## Centroiding on a five-scan superspectrum

```text
qAlgorithms -i profile.mzML -o output --superspectrum
```

Each polarity is processed separately. The first five MS1 profile scans
form the initial window, then each new scan advances it by one.

The automatic mode samples up to 24 evenly spaced scans across each
polarity. Adjacent m/z spacings from the original scans estimate local
sampling density. A spacing larger than twice its smallest adjacent
spacing is excluded from this estimate. Four mantissa bands per octave
give fixed m/z regions. In each region the grid step is half the median
spacing. The gap limit is the larger of three times the median and 1.5
times the 90th percentile. Regions with fewer than 24 spacing observations
use pooled estimates. If the initial sample has no spacings, the scan list
is checked for the first usable scan.

The regression half-width is 1.5 times the 90th percentile of measured
full widths at half maximum, rounded up in grid samples and clipped to
the compiled range. Widths come from maxima whose falling flanks reach
half height on both sides within a covered run. Fewer than 12 measured
widths use the compiled upper bound. The run reports the estimated
settings and includes the actual grid step in `gridStepMz`.

To reproduce a fixed-settings comparison, use
`--superspectrum-fixed GRID_STEP MAX_GAP MAX_HALF_WIDTH`.
`GRID_STEP` and `MAX_GAP` are in m/z, and the half-width is in grid samples.

Adjacent positive profile points within the gap limit at their midpoint
define covered runs.
Their union defines the superspectrum segments before interpolation.
Common intervals supported by at least two runs from distinct scans
are projected onto a fixed lattice `index * gridStep`. Fits near regional
boundaries can use support from either side. The apex region determines
which lattice owns the result. Interpolation is
linear in intensity and confined to each scan's covered runs. Resampled
runs stay cached until their scan leaves the FIFO. The fit uses the
precomputed piecewise quadratic qPeaks inverse and its running sums along
m/z. Each window estimates a shared shape and a separate height per scan.

The experimental `*_superspectrum_centroids.csv` has a row per contributing
scan, with the window bounds and centroid ID as the grouping key.
Concave fits must have an interior apex, a fitted height at least twice
both fitted edge heights, and explain more than half of the mean log
trace's variation in the fit window. Overlapping candidates are ranked by
explained variation times the contributing scan count and grid step.
These heuristics
operate on the shared shape. A centroid can recur in successive windows.
The window CSV includes the assigned track ID. Linking requires a common
contributing scan, mutual coverage of the fitted apices and an m/z
distance within the larger of two grid steps or one quarter of the
shorter fit interval. Each window contributes at most one centroid to
each track. Eligible links are assigned in order of increasing m/z
distance. A track can resume while the FIFO windows still overlap.

`*_superspectrum_tracks.csv` contains one row per track and contributing
scan. If several fits include the same scan, its height comes from the
window in which the scan is closest to the middle. Ties use the earlier
window. The mean and sample standard deviation use one fitted m/z
position per linked window. A single-window track has an empty
standard deviation. The standard deviation describes variation among
overlapping windows, which share profile data.
