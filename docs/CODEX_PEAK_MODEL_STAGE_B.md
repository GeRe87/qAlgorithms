# Codex implementation briefing: modular PeakModel binding (Stage B)

Repository: `GeRe87/qAlgorithms`.
Status: specification only; this document does not authorize an implementation,
merge, release, or production quantification approval.
Stage A implementation baseline: `1039ad0eceb7b74822b252ecb354140e8a9eae0f`.
Stage A Draft PR: https://github.com/GeRe87/qAlgorithms/pull/5
Verified Stage A CI: https://github.com/GeRe87/qAlgorithms/actions/runs/42842039115

## Objective and branch boundary

When Stage B is separately commissioned, bind the existing hyperbolic running
kernel to the typed model interface and expose explicit, opt-in fit diagnostics.
Keep default piecewise-quadratic detection, centroiding, EIC/feature construction,
CLI, data layouts and persisted output byte-compatible with Stage A.

Create `agent/peak-model-binding-stage-b` from the exact reviewed Stage A HEAD.
If Stage A has advanced, record the new starting SHA and review its diff first.
Open a separate Draft PR targeting `agent/peak-model-contract-stage-a` while the
stack remains unmerged. Do not retarget or modify PRs #1–#5, their base branches,
or the research branch. No merge. Do not start Stage C implicitly.

## Decisions for the implementation

### 1. Bind fitting before statistical peak detection

Use `PeakModelKind::Hyperbolic` in a NEW opt-in fitting/diagnostics entrypoint.
It returns model-tagged fit results, not accepted `RegressionGauss` or
`FeaturePeak` objects. Support one aligned trace as m=1 of the same multitrace
kernel, and multiple traces with shared b and local a,C,D.

Keep existing `qpeaks_find` overloads and their default behavior unchanged.
In particular, `qpeaks_find(..., PeakFitOptions{Hyperbolic}, ...)` continues to
throw before processing or mutation during Stage B. This deliberately refines
the earlier shorthand 'Stage B binds Hyperbolic': the binding is usable for
explicit fits, but statistical detection cannot honestly accept Hyperbolic
before Stage C supplies its gates. Document this distinction in public headers.
Do not silently route Hyperbolic into the legacy detector.

A minimal API sketch (names may follow repository conventions):

- A request owns/borrows aligned finite log traces and matching raw traces,
  length, maximum half-width, and an externally fixed dimensionless step.
- The request has an explicit model enum, with no default change to old APIs.
- Dispatch once at entry to a templated implementation; no per-window enum switch.
- A streaming visitor receives a borrowed model-specific fit view with documented
  lifetime/order. An optional collecting adapter explicitly owns copied results.
- Use distinct coefficient/geometry types for each specialization; never encode
  hyperbolic a,b,C,D in RegCoeffs b0,b1,b2,b3 or a legacy six-value inverse.
- Reject unknown enum values and malformed requests before invoking the visitor.

Keep model-neutral orchestration small. Do not generalize legacy statistical
functions into a shared base class. No virtual calls, generic matrix inversion,
or per-window allocation in the running kernel. Prepare geometry/scratch once.

### 2. Fixed dimensionless coordinate, explicit units

Stage B uses only `u = fixedStep * (sampleIndex - windowCenter)` and aligned,
uniform grids. `fixedStep` must be finite, positive, supplied by the caller and
fixed independently of the fitted signal. No fitted lambda, grid search,
resampling, imputation, automatic instrument transform, or data-driven step tuning.

Return apex/width in dimensionless u and analytic full area as intensity*u.
Keep window-center sample index separately. Never label these values m/z, RT,
physical intensity*m/z, or intensity*time. Do not map nonlinear axes using one
constant Jacobian. Physical conversion and its integration/uncertainty rules are
a future explicit design step; do not infer calibration from research defaults.
For a future proven affine axis x=x0+s*u, area conversion is |s|*areaU; Stage B
need not expose that conversion.

### 3. Model capabilities and result validity

Add a Hyperbolic specialization using the existing research equations and cached
geometry. Avoid independent copies of the kernel: prefer a single internal math
implementation used by both the binding and the research wrapper, with a small
compatibility shim if relocation is needed. Keep research executables working.
Production headers must not require an experiments/ path. Any move must preserve
research API behavior and have a mechanical diff plus independent regression proof.

Represent capability availability separately from numerical values. A fit result
must distinguish at least:

- finite coefficients and identifiable geometry;
- integrable tails: C>|D| (equivalently both exponential tail coefficients >0);
- finite analytic apex, and separately whether it lies in observed support;
- valid raw-amplitude refit;
- available finite whole-line area and its numerical evaluation status;
- statistical acceptance: unavailable until Stage C.

Unidentifiable/unintegrable fits can remain diagnostic results but cannot become
accepted peaks. A full integral may exist even when the apex is outside observed
support; do not combine these statuses. A tiny/nonfinite Bessel result must be an
explicit numerical property failure, never zero area or a fabricated valid peak.
No new arbitrary shape thresholds beyond the documented kernel guard without
evidence and documentation. Preserve the old model's historical mixed-tail rules.

Stage B need not implement a new analytic FWHM formula. Mark hyperbolic width
unavailable unless an independently verified method is explicitly included in
scope. Missing widths or uncertainties must not appear as valid zeroes.

### 4. Raw-amplitude refit and covariance boundary

Apply the existing normalized scalar raw-intensity refit only to explicitly
selected windows, using matching raw samples. Never repeat it for every running
window. Preserve b,C,D and fit geometry; recompute height and area from corrected a.
Keep original log-OLS coefficients/residual diagnostics distinct from corrected a.
Return invalid status for negative/nonfinite raw data, zero signal or bad geometry
according to the existing kernel's documented behavior. Do not manufacture logs
for zero raw intensities. Validate finite supplied log samples before running fit.

Expose coefficient covariance factors with their parameter order, trace identity,
shared cross-trace terms, and homoskedastic log-noise assumptions. A factor is not
a standard error. The kernel's corrected-log-amplitude influence variance is not
full corrected-area variance. Mark postfit area uncertainty unavailable.

Stage C alone defines pooled residual variance, model-specific statistical gates,
full postfit property gradients/covariance and peak acceptance. Do not copy legacy
F/t thresholds, degrees of freedom, baseline/noise or ambiguity decisions.

## Mandatory acceptance evidence

1. Record original Stage A HEAD, working tree state, compiler/hardware and full
   baseline Debug/Release CMake/CTest results before production changes. Baseline
   is 12/15 passing; test_qpeaks, test_qpeaks_pt03, test_qpeaks_pt04 are existing
   failures (SIGTRAP on Linux). Do not hide them or modify their expectations.
2. Repeat Stage A exact differential checks for all running coefficients/order,
   retained fields, raw-corrected heights, areas, widths, uncertainties, spans,
   rejection/selection behavior and group metadata. Default/explicit legacy paths
   must agree. Compare synthetic mzML centroid/EIC/feature output files and CLI
   behavior against separately built Stage A binaries; add fixtures only as needed.
3. Compare the hyperbolic binding to independently assembled dense full OLS and
   inverse covariance for m=1,3,10 and multiple steps/scales. Check shared and
   cross-trace factors, parameter ordering, borrowed-view lifetime and output order.
   A comparison to the same kernel alone is insufficient independent evidence.
4. Check analytic apex against an independent derivative/root calculation and
   Bessel area against independent numerical quadrature. Include near-tail-boundary,
   unidentifiable geometry, nonfinite inputs, extreme amplitudes and apex-outside-
   support cases. Define justified per-property tolerances before assessing results;
   legacy differential tolerance remains zero on each matched toolchain.
5. Validate raw refit against direct one-parameter least squares: raw RSS cannot
   increase within justified rounding tolerance, shape/apex unchanged, area scaling
   matches exp(delta a), original log fit retained. Reuse independent finite-difference
   influence checks; do not claim full area uncertainty is validated.
6. Check that legacy qpeaks_find still rejects Hyperbolic before input processing
   or mutation; new fitting API explicitly dispatches Hyperbolic and rejects invalid
   models/requests. No routing to centroid/feature production pipelines or CLI.
7. Debug/Release CI must build production and isolated research targets, run all
   new independent checks and report the complete suite separately from the known
   failures. Upload exact SHA, compiler flags, logs and differential summaries.
8. Benchmark identical standalone sources against separate Stage A/Stage B builds
   on the same hardware/compiler. Cover single trace and m=1,4,16 multitrace,
   equal windows/scales; separate geometry preparation, hot fits, output collection
   and selected-window postfit costs. At least nine alternating-order paired trials;
   report median, spread and samples, not only best time. Do not mix API test code
   into timed executables. Repeat on isolated hardware if variability prevents a
   conclusion. A reproducible legacy slowdown >2% requires investigation and an
   explicit unresolved acceptance item; timing within that band is not proof of
   exact equality. Characterize binding overhead separately from inherent model cost.

## Delivery and stop conditions

Deliver the separate Draft PR with exact base/head SHAs, changed files, API/unit
examples, evidence, performance table and remaining questions. No merge/release.
Stop and report a concrete design blocker if binding requires fabricated legacy
results, reuse of legacy statistical gates, implicit physical calibration, or a
change to default legacy behavior. Do not quietly expand into Stage C.

Remaining future decisions: physical axis/calibration/Jacobian semantics; missing
or unequal-weight samples; statistical selection and interference; full corrected
area covariance; robust scaled-log Bessel evaluation across the required numeric
range; capability/status API stabilization. Stage B must identify unsupported
ranges explicitly rather than claiming certified physical quantification.
