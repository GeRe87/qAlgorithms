# Experimental qPeaks regression ensembles

This is an inspectable research path, **not a replacement peak detector**.
Production source, statistical thresholds, collective-output policy, feature
CSVs and default build are unchanged. No measurement data are committed here.

The [one-versus-two comparison](COMPARISON.md) adds matched component-family
fits, complete start accounting and hinge sensitivity. It makes no peak-count
release decision.

The [second checkpoint](REFINEMENT.md) adds optional centered arithmetic and a
shared-curve point estimator on prespecified supports. It does not yet make
ensemble peak decisions. The original-kernel experiment below remains the
first checkpoint; its default collection mode is retained.

## What is implemented

`collect.cpp` includes the existing qPeaks implementation in a separate
executable. It collects every precomputed moving-regression fit on complete
contiguous **positive** input blocks. Zero samples and positive blocks shorter
than the five points needed by the four-coefficient model are counted
explicitly; nothing is silently omitted. This intentionally studies fits below
the production intensity floor too. The original x coordinates and measured-
point counts are retained for comparison with the current strict API.
All ensemble topology uses equidistant sample indices; physical RT is not
interpolated, regularized or used to redefine the regression kernels.

Every fit is retained whether or not it has an interior maximum or passes
production validation. Recorded fields include the unmodified coefficients,
full-window log RSS, exact constant-log-window flag, geometric maximum/apex,
production-floor coverage and current local-cascade acceptance. No raw-height
correction is applied to the exploratory coefficients. Current acceptance is
annotated by running the existing validation separately on its own above-floor
blocks; it is not used to build the experimental graph.

The first graph connects maximum-bearing fits at adjacent enumerated centers
or adjacent enumerated scales. This four-neighbour lattice is an explicit,
inspectable diagnostic design choice, **not a justified same-peak relation**.
Diagonal motion without an intervening eligible model can break a continuation.
Connected components can merge physical peaks or remain extensive in noise.
No component, birth or long-lived branch is declared a detected peak.

`analyse.py` contracts the existing same-scale edges into contiguous center
runs. All vertical connections between runs remain in a directed acyclic graph,
including merges, splits, births and deaths. Shared-center edge counts are
bookkeeping, not independent confirmations. No minimum run length, minimum
number of fits, apex-distance cutoff, persistence threshold or new p-value is
introduced. Run mean apices are **display coordinates only**; their spans are
not uncertainty intervals or optimized peak estimates. Plot opacity, marker
sizes and repeating colour palettes have no role in analysis decisions.

## Reproduction

Build the normal repository library first. From the repository root with GCC
and zlib available (commands shown for a Release library):

```text
g++ -std=c++20 -O3 -DNDEBUG -Iinclude -Iexternal/libcerf experiments/qpeaks_ensemble/collect.cpp build-release/libqalgo_lib.a -lz -o ensemble-collect.exe
python experiments/qpeaks_ensemble/test_collect.py ./ensemble-collect.exe
./ensemble-collect.exe /path/to/inputs.bin /path/to/study/release
python experiments/qpeaks_ensemble/analyse.py /path/to/study /path/to/analysis
python experiments/qpeaks_ensemble/test_analysis.py
```

On Windows, use the corresponding executable path syntax. Analysis and plotting
require NumPy and Matplotlib; the collector and its six tests use only the
existing C++ dependencies and Python standard library. A separate NumPy install
location can be passed with `--numpy-root PATH`. `--no-plots` skips Matplotlib.
The checked reference environment used Python 3.14, NumPy 2.5.4 and Matplotlib
3.11.2. No experiment is added to the default CMake targets.

The local evidence package supplies `inputs.bin`, `scenarios.json`, collected
CSVs and the analysis outputs. To repackage the preceding paired benchmark
with optional real EIC observation CSVs, use:

```text
python experiments/qpeaks_ensemble/prepare_inputs.py /path/to/preceding-paired-data /path/to/study --eic /path/to/trace.csv
```

Repeat `--eic` for additional traces. Observation CSVs require `x,y,df` columns;
df must count measured points cumulatively. Provenance and measurement values
are written only to the requested local study directory. The preparation
script reproduces every byte of the checked input binary. The analysis expects
the scenario IDs of that preceding paired protocol, not arbitrary new IDs.

Input binary: ASCII `QENS0001`, followed by repeated little-endian records:
uint32 scenario, trial, n, truth_count; float32 truth[truth_count], x[n], y[n];
uint16 cumulative_measured_points[n]. Nonfinite/negative y, non-increasing x,
unsupported counts, truncation and nonfinite fits fail explicitly.

## Outputs and accounting

- `traces.csv`: one row per input, including every successful strict output or
  explicit strict processing error. Errors have count -1, never successful zero.
- `components.csv`: every diagnostic graph component for every trace.
- `candidates.csv`, `edges.csv`, `observations.csv`: full detail for **trial 0**
  of each scenario. This is deterministic output selection, not a peak filter.
- `runs.csv`, `lineage-edges.csv`: contracted lineage graphs for those detailed
  traces, retaining all original vertical graph connections.
- `analysis.json`: complete-input accounting, scenario summaries, reflection
  comparisons, lineage event counts and the correlated-ensemble reference.

The audit independently checks the candidate enumeration from all positive
blocks, record keys, every detailed node/edge/component, finite values and
accepted-fit bookkeeping. The six collector tests include an independent grid
BFS, reflection, unit conversion, point accounting, exact-flat recording and
explicit invalid-input errors. Five analysis tests check merges, splits,
reflection/empty fields, overlap covariance and the common-model GLS identity.

## Checked first experiment

The local study contains 5220 traces: the preceding 5202 paired inputs,
ten reflected trial-zero controls, six previously captured real EICs and two deterministic controls.
There are 6502840 fitted windows in total; 524349 windows have full CSV detail
across 28 selected traces. The synthetic protocol includes 4000 specified null
traces, 1200 noisy signal traces, and the unchanged pt_03/pt_04 fixtures.
The positive flat baselines and zero-clipped raw noise have the interpretation
limits documented in the preceding statistical audit. They are controlled
inputs, not an instrument-noise calibration or real-file FDR estimate.

The collector completes in Release and assertion-enabled -Og Debug with all
**five CSV outputs byte-identical**. Its six structural tests also pass in
standard unoptimized Debug. All five Python analysis tests pass. This does
not claim a new whole-file mzML run or a fully implemented ensemble estimator.

Important observations:

* At the smallest scale, pt_04 contains three separate center runs with fitted
  apices around 3, 7 and 11. The first run has no fit accepted by the current
  local cascade. The later two runs join at scale 3, so the complete lattice
  has only two components. Thus connected-component count is not peak count.
* The overlapping-pair example has two initial branches which join at larger
  scales; its complete graph has one component. The three-Gaussian scenario
  likewise has a median component count of one across its 300 inputs.
* Specified flat log/raw noise also produces extensive multiscale structures:
  all 1000 traces in either generator have graph components, with a median
  largest-component coverage of all 48 available scales. Persistence alone
  is therefore not adopted as evidence. These comparisons are descriptive;
  different trace lengths and positive-block extents affect available scales.
* For an exactly constant positive 99-point trace, mathematical least squares
  has zero slope/curvatures. The current finite-precision kernel nevertheless
  produces 2301 formal interior maxima among 2304 windows, forming two graph
  components. Exact flatness is recorded for every window, without a tuned
  near-zero tolerance. The current cascade emits zero peaks for this control.
* All ten reflected examples preserve maximum classifications and lattice
  topology. Raw fitted apex coordinates are not exactly reflection-invariant:
  the largest discrepancy in these examples is about 0.001086 sample in
  pt_03. Very weak curvature amplifies kernel arithmetic error. No numerical
  reflection correction is made in this experiment.
* The six real traces still raise their explicit **strict API** conflict
  errors. This is expected: the optional collective-output mode is not enabled
  by this collector. It does not revoke the complete production run with
  explicitly marked collective descriptions from the preceding checkpoint.

## Why correlated ensembles cannot create independent information

For a fixed local window, write the ideal-arithmetic log-domain estimator as

    beta_hat_j = A_j z.

With common independent error variance sigma^2,

    Cov(beta_hat_j,beta_hat_k) = sigma^2 A_j A_k'.

With general measurement covariance Sigma, replace this with A_j Sigma A_k'.
Disjoint supports have zero cross-covariance only under the independent-error
reference. Selected groups do not automatically inherit this unconditional
fixed-window law. Peak positions are nonlinear functions of coefficients;
their covariance additionally needs a justified transformation/approximation.

The Python QR/SVD calculations are diagnostic reference work, not new per-scan
inversions in the production kernel or a production speed claim.

The analysis implements a narrow, analytically checkable reference: a common
**global quadratic log mean**, i.e. one Gaussian on a predetermined support,
with all local windows included before any data selection. Both local quadratic
curvatures then equal the global quadratic curvature. This is not an assumption
that every observed bundle follows one Gaussian.

Stacking A_j gives a covariance sigma^2 A A', often singular. GLS is computed
through the observable row-space projector P=A^+ A, avoiding inversion of this
singular covariance. For common design X, the estimator is

    theta_hat = (X' P X)^(-1) X' P z.

For the checked supports, A has full column rank and P=I. The correctly
correlated ensemble therefore equals the direct original-data OLS estimator,
to about 3.4e-16 in its linear operator. This is a mathematical reference,
not an automatically selected peak fit and not an uncertainty-coverage study.

| Original points | Local fits | Stacked coefficients | Covariance rank | Nullity |
| --- | ---: | ---: | ---: | ---: |
| 9 | 9 | 36 | 9 | 27 |
| 21 | 81 | 324 | 21 | 303 |
| 41 | 361 | 1444 | 41 | 1403 |

Ignoring inter-fit covariance changes the estimator to a repeated-data
weighting and understates its uncertainty. In the 41-point reference, the
actual slope standard error is about **11.12 times** the naively reported
one. The naive estimator also has larger actual variance than direct OLS.
An ensemble can help with support choice, bias and competing descriptions;
this reference shows why counting overlapping fits as new observations is
not that improvement.

## What remains before an ensemble peak decision

The lattice/lineage is an exploratory baseline. A defensible association
criterion must address diagonal continuations, numerical degeneracy, noise
branches and the meeting of neighboring peaks. A final common model must be
specified; coefficients with different hinge positions cannot simply be
averaged into another member of the current piecewise-quadratic peak family.
Prediction/model comparison must use a common original-data domain and a
stated noise model, counting each measurement once. If multiple physical peaks
are fitted together, their intensities add in the raw domain; their individual
log profiles do not add linearly.

Calibration must include support/search/branch selection, not merely the
winning fit. No persistence cutoff, ensemble vote, best residual or small
within-bundle spread is currently promoted to a statistical test. Areas,
interpolation and physical RT corrections remain deferred. The production
pt_03/pt_04 expectations remain unchanged and unresolved.

Related primary-method precedents, not drop-in guarantees for this prototype:
[SiZer](https://doi.org/10.1080/01621459.1999.10474186) and
[multiscale mass-spectrum peak detection](https://doi.org/10.1093/bioinformatics/btl355).

The [real mzML extension](REAL_DATA.md) applies the conditional component
comparison to freshly exported real EIC windows, with complete input accounting
and explicit unconverged fits. It does not turn diagnostic fits into peak counts.
