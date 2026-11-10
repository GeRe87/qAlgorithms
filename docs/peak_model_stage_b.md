# Stage B: opt-in model fitting diagnostics

## Scope and API

Stage B starts at `6a1a4c39b213fa1d0ff84a611f3e1a8da9620b7f`, a clean tree.
The only advance from verified Stage A `1039ad0` was its Stage-B specification.
The old detection API, CLI and persisted types remain unchanged. In particular,
`qpeaks_find(..., PeakFitOptions{PeakModelKind::Hyperbolic}, ...)` still throws
before processing; it does not perform statistical peak detection.

Include `qalgorithms_peak_fit.h` for the new explicit diagnostics API:

```cpp
using namespace qAlgorithms;
PeakFitRequest request{PeakModelKind::Hyperbolic, {logTrace}, length, 24, 0.1};
forEachPeakFit(request, [&](const PeakFitView &view) {
    const auto &fit = *std::get<const peak_model::HyperbolicModel::Fit*>(view);
    // Borrowed fit and locals are overwritten by the next callback.
    // Copy fit if this is an externally selected window to retain.
});
```

`forEachPeakFit` validates all finite log samples, dimensions, visitor, model and
geometry before invoking any callback. It dispatches once per call. The
Hyperbolic branch prepares a cache once and uses reusable kernel scratch.
`Model<Hyperbolic>::runningFit` accepts a prepared typed cache for repeated fits.
The callbacks are center-major/increasing-scale, including unconstrained fits;
they do not imply accepted peaks. The variant holds the appropriate typed pointer.

The PiecewiseQuadratic branch exposes the existing shared-shape regression for
m=1 or multiple traces with step=1 in sample coordinates, filtering its existing
full-scale stream to maxHalfWidth. It intentionally does not change that kernel
or introduce a new early-stop optimization. Existing production paths never use
this diagnostics filter.

There is one mathematical implementation in
`include/detail/qalgorithms_hyperbolic.h`. The old research header is a namespace
alias/include shim. Apart from namespace spelling, the moved kernel is identical
(mechanically checked in CI). All independent research tests and tools remain.
No hyperbolic equations were copied or rewritten in the binding.

## Selected-window properties and units

```cpp
peak_model::HyperbolicModel::Geometry geometry(24, 0.1);
// 'saved' is an owned copy from the same fixed-step request.
auto props = peak_model::evaluateHyperbolicWindow(saved, traceIndex, geometry,
                                                rawTrace, length);
```

The caller must use the SAME step/cache geometry and aligned raw samples as the
original fit. Sample-array lengths are caller-provided pointer contracts; the API
cannot discover allocation sizes or trace provenance. Passing nullptr raw
explicitly requests properties of the uncorrected log fit. Invalid supplied raw
data does not silently fall back to uncorrected height/area. The fit is immutable;
corrected intercept and amplitude-refit result are separate fields.

The external fixed positive step defines dimensionless u relative to each window
center. Apex/width use u, area uses intensity*u, height inherits raw-intensity
units. The sample center remains separate. No instrument transform, inferred
calibration, fitted step, physical-axis conversion or Jacobian is implemented.

Finite coefficients, identifiable geometry, integrable tails, apex-in-support,
raw-amplitude validity and numerical property availability are independent flags.
Numerically representable diagnostic properties may exist for an unidentifiable
fit; always inspect identifiable separately. Nonintegrable tails and Bessel
failure have explicit statuses. An invalid requested refit leaves corrected
height/area unavailable. Width, area uncertainty and statistical acceptance are
unavailable (NaN plus status), never valid zeroes. No legacy F/t gates are reused.

Covariance factors are exposed with order `[b,a0,C0,D0,a1,C1,D1,...]` and include
shared cross-trace D covariance. They assume fixed aligned grid, independent
homoskedastic log observations and require an externally justified residual
variance before being interpreted as covariance. The refit influence variance
is not corrected-area variance. Pooled noise estimation and acceptance are Stage C.

## Validation and reproduction

Windows Intel i7-11850H / MSYS2 GCC 15.1.0; native arch OFF. Before production
edits, Debug and Release builds with research ON succeed: 13/16 tests pass.
After Stage B, 14/17 pass in both. Exactly the old failures remain:
`test_qpeaks`, `test_qpeaks_pt03`, `test_qpeaks_pt04`. Research OFF has one fewer
passing research test. No existing test expectations were changed. The research
test main now explicitly returns 0 so its independent helpers can also be reused
in the binding test without falling off a renamed non-main function.

Independent full dense OLS and full covariance checks exercise the binding for
m=1,3,10, steps .06/.12/.25, scales 8/16/24. Relative-plus-absolute tolerances
are 1e-7 for coefficients and 1e-6 for covariance, inherited from the independent
research reference; derivative-root apex 1e-10, Simpson area 1e-7, area correction
scaling 1e-8. Legacy differential is exact (zero tolerance). Existing research
quadrature, scalar raw LS, extreme amplitude, and cross-trace finite-difference
influence tests are reused and run independently too. Additional guards cover
finite-input prescanning, invalid enums, zero step, visitor lifetime/storage,
zero/negative/nonfinite raw input, invalid and near-boundary tails, unsupported
area numeric range, poor identifiability and apex outside support.

Local separate-library differential: 5,646,223 serialized bytes identical;
23 synthetic mzML feature/group CSVs identical. This exercises centroid/EIC/feature
orchestration, not an exhaustive real-instrument accuracy corpus.

Configure with research ON and optional binding benchmark:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQALGORITHMS_ENABLE_HYPERBOLIC_RESEARCH=ON -DQALGORITHMS_ENABLE_PEAK_BINDING_BENCH=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Full CTest intentionally returns failure for the three documented legacy tests.
The dedicated CI checks the exact failed-test set and uploads complete logs,
compile commands, compiler/CPU provenance, legacy differentials and paired timing.
It independently checks out immutable Stage A and builds both versions.

`experiments/peak_model_binding/compare.py OLD_PROBE NEW_PROBE OLD_BENCH NEW_BENCH`
compares exact output and runs nine alternating-order paired trials. Benchmarks
use identical standalone sources, 129 samples, identical 3969 multitrace windows,
20 passes per measurement. Compare only matching kernel/window counts. Single-
trace collecting fit and m=1/4/16 streaming fits are measured separately. Geometry
preparation, direct research hot fits, raw amplitude postfit, prepared binding,
full request (validation/cache), explicit collection, and selected properties
(including Bessel area) have separate labels/costs. Units are ns/window except
geometry (ns/cache preparation) and selected postfit (ns/trace).

Local legacy Stage-B/Stage-A median ratios: collecting single .9707; streaming
m=1 .9719, m=4 .9796, m=16 .9781. No local legacy slowdown detected. Release
uses -O3 -DNDEBUG plus existing supported CMake diagnostic/math flags. Final CI
numbers and exact tested HEAD are recorded in the Draft PR, with variability.
These are repeated shared-host measurements, not a zero-overhead proof.

## Deliberately outstanding

No production detection integration or quantitative certification. Stage C must
define statistical gates/noise/DOF, ambiguity handling and full corrected-property
covariance. Physical calibration/Jacobian, irregular or missing samples, weighted
fits, robust scaled-log Bessel evaluation and capability/status stabilization
remain explicit future work. Current numeric-range failures are reported rather
than approximated. Match cached step and raw-trace provenance at the call site.
