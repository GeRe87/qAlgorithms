# Stage A architecture and evidence

Starting source HEAD: `33cf85e472559bf4ee7cfd2c83d904a3371d1973` (briefing only).
Stacked target / immutable pre-briefing base: `8c6b51087205e11d9c09296e74b8a1f8ecf5ea59`.

The original qpeaks_find symbol remains unchanged. The options overload switches
once before input processing. Hyperbolic and invalid enum values throw
invalid_argument without touching result vectors, even for invalid input.
The incomplete Model primary template prevents accidental reuse of legacy math
for another model. Its sole specialization uses RegCoeffs, the original running
kernel, six-entry inverse cache, raw correctB0 and analytic property/gradient
functions. There are no new per-window allocations, inversions or virtual calls.
PeakProperties is an explicit read adapter; callers specify axis units. Intensity
units are inherited from the caller. Width uncertainty is unavailable in legacy.
No existing layout, serialization, CLI options or output format changed.

This intentionally uses the briefing's permitted interface/adapter vertical slice.
Math implementation bodies, cached geometry, covariance gradients and F/t gates
remain coupled inside qalgorithms_qpeaks.cpp. The interior-apex support gate is
named separately; it is not a universal tail-integrability predicate. The legacy
integral supports historical mixed tails: imposing a new two-concave-tail test
would change acceptance. Candidate support and legacy statistical validation stay
in their original order. Do not instantiate these statistical gates for Hyperbolic.

## Local validation (Windows 11, Intel i7-11850H, MSYS2 GCC 15.1.0)

Before production edits: Debug and Release full CMake builds succeed; CTest 10/13
pass in each. After: full builds succeed; 11/14 pass in each. Exactly the known
failures remain: test_qpeaks, test_qpeaks_pt03 (expected 0, found 2),
test_qpeaks_pt04 (expected 3, found 2). No expectations changed.

Separately linked baseline-source/library versus refactor probe: 40 deterministic
single/asymmetric/overlapping/noisy profiles, lengths 17..134, with/without DF and
collective descriptions. Every running coefficient/window order and every retained
RegressionGauss field/span/uncertainty plus group metadata serialized using hexfloat:
5,646,223 bytes identical in both Debug and Release, zero tolerance. Explicit and
default APIs checked field by field; unsupported selection checked before processing.
Existing independent dense-OLS multitrace and folding tests pass in both builds.
Synthetic mzML integration fixtures pass against both binaries; generated feature
and unresolved CSV files compared byte-for-byte separately. These fixtures exercise
centroid/EIC/feature orchestration; this is not exhaustive real-instrument validation.

Release paired running-fit benchmark: 512 identical samples, all scales/windows,
30 passes per measurement, 5 alternating process-order pairs x 11 measurements.
Baseline median 11.0127 ns/window (min 10.6182, max 12.9219, SD 0.4133);
adapter median 11.4038 (min 10.5678, max 13.0960, SD 0.4890), ratio 1.0355.
Timing includes vector collection and consumption, excludes peak validation and
postfit; single trace only, multitrace implementation unchanged. Library flags:
Release -O3 -DNDEBUG, CMake supported diagnostic/math flags (including
-fno-math-errno -fno-signed-zeros), native arch OFF; probe -O3 -DNDEBUG.
The 3.55% median increase lies within overlapping observed timing ranges; do not
claim proof of zero overhead from this noisy shared-host measurement. The running
fit body is unchanged and the adapter adds one call per complete batch, not window.
CI repeats the independent source comparison and timing on Linux Debug/Release;
its uploaded evidence must be consulted separately from these local measurements.

Reproduce: build baseline at the recorded source SHA after copying the current
`test/test_peak_model.cpp` into its test directory (without STAGE_A defined), build
current HEAD, then run `python test/stage_a/compare.py BASELINE_PROBE CURRENT_PROBE`.
The Stage-A CMake test target automatically defines STAGE_A.

## Stage B/C open decisions

- Bind research kernel with its own coefficient/geometry types only in Stage B;
  decide single/multitrace production entrypoint and physical axis/Jacobian first.
- Separate true shape integrability from observed support; preserve legacy semantics
  while defining Hyperbolic rejection reasons and model-specific statistical DOF/gates.
- Stage C must propagate full post-amplitude-refit covariance into area, including
  shared coefficients; current legacy uncertainty limitations remain intentionally.
- Decide optional property/capability and status taxonomy beyond the minimal adapter;
  no fabricated width uncertainty or shared statistical thresholds.
- Wider real-data golden centroid/EIC/feature corpus and isolated CPU performance
  measurements remain useful follow-up evidence, not certified quantification.

Hyperbolic production integration is absent. PRs #1–#4 and their branches are untouched.
