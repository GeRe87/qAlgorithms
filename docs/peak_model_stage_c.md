# Stage C: selected-window inference and corrected-property covariance

## Boundary and code style

Base: `b00d894bde26c6260ffe06653932465d03b9e8ea` (Stage B, clean tree).
Branch: `agent/peak-model-inference-stage-c`; separate Draft PR against
`agent/peak-model-binding-stage-b`. No edits to PRs #1–#6, their branches,
legacy detector/CLI/output formats, or the shared mathematical kernel.

New C++ follows the existing four-space indentation, brace placement and
camelCase function naming. Changes are in separate inference header/source,
independent tests, optional benchmark, CMake opt-in target and dedicated CI.
Existing compact research code is not reformatted or mixed into this diff.

This implements a selected-window statistical/uncertainty layer, not scan-wide
peak detection or production quantification approval. `qpeaks_find(Hyperbolic)`
continues to fail before processing. Stage B APIs and their unavailable-statistics
semantics remain unchanged. No implicit adaptation into legacy RegCoeffs or outputs.

## API and assumptions

Include `qalgorithms_peak_inference.h`:

```cpp
using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;
HyperbolicInferenceOptions options;
options.selection = WindowSelection::prespecified;
options.familySize = numberOfWindowsDeclaredBeforeObservingData;
auto result = inferHyperbolicWindow(rawTraces, length, center, halfWidth, fixedStep, options);
```

The function derives logs from the SAME positive finite raw samples. Zero,
negative, clipped, interpolated or independently supplied log/raw pairs cannot
be substituted into this inference contract. Only the selected raw window is
read. Fixed step is external and dimensionless. The original running kernel is
reused without changes on this local window; only the central full-window fit
is retained. Extra local regression windows are discarded. This selected-fit
cost is outside the hot running API and is benchmarked explicitly.

Reference random model: independent homoskedastic Gaussian errors in log
intensity, fixed aligned grid, full-rank design, predetermined window and trace
pool. Full-window residuals come from the ORIGINAL unconstrained log-OLS fit,
before raw-amplitude refit. No model checking procedure can establish these
assumptions merely by taking logarithms. The single-precision fitting path also
has its normal rounding limitations; it is not an exact real-arithmetic solver.

With n=2h+1 and m traces, full model parameter count is 3m+1 and

    nu = mn-(3m+1) = m(n-3)-1
    sigmaHatSquared = RSS_full/nu.

The nested null is a_j + b*u (individual intercepts, common log-linear trend),
with m+1 parameters. Thus the added hyperbolic terms have 2m degrees of freedom:

    F = ((RSS_null-RSS_full)/(2m)) / (RSS_full/nu).

This tests additional curvature over that specific common-trend null. It is
NOT a test of model correctness, per-trace significance, independent chemical
identity, or absence of interference. In particular pooled evidence does not
establish each trace as an individually detected analyte. No trimmed residual
span or mechanically copied legacy degrees of freedom/F/t thresholds are used.
The model-specific counts and nested comparison above are derived for this API.
The general fixed-design least-squares reference is
[NIST](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd431.htm).

The F survival is I_x(nu/2,m), x=nu/(nu+2mF). Since m is integer, it is computed
as a finite positive series in log space. Both log p-values and ordinary p-values
are returned; the gate compares logarithms so p-value underflow is not treated
as evidence of exact zero probability. Bonferroni uses the caller's declared
complete family: adjustedP=min(1,familySize*p). It permits dependence between
fixed-window tests; it does not repair selecting the tested family from data.

The default is exploratory: nominal reference statistics and covariance can be
inspected, but groupDecision remains unavailable. Only an explicit prespecified
window declaration enables the conditional group gate. A passing group needs
finite/identifiable/integrable fits, valid refits/properties, sufficient observed
support on BOTH sides of every apex, and adjusted p <= declared family alpha.
The support policy is caller-visible (default two samples/side); it is not
another significance test. This is no automatic ambiguity/interference gate.
If those scientific conditions are unresolved, do not interpret the decision
as production peak acceptance. Failed numeric/shape inference has no covariance.

## Complete first-order postfit influence and cross-trace covariance

For a trace's fixed fitted shape s_i=b*u_i-C*cosh(u_i)-D*sinh(u_i), normalized
raw-amplitude weights and squared-shape weights are

    w_i = I_i exp(s_i) / sum(I exp(s)),
    v_i = exp(2s_i) / sum(exp(2s)).

The raw refit a*=log(sum(I exp(s))/sum(exp(2s))) has shape gradient

    g_b = E_w[u]-2E_v[u],
    g_C = -E_w[cosh(u)]+2E_v[cosh(u)],
    g_D = -E_w[sinh(u)]+2E_v[sinh(u)].

Its derivative also includes the direct observed-intensity term w_i. Treating
a* as independent of shape, or merely adding Var(a*) to an area gradient, would
omit essential cross terms. The implementation differentiates through BOTH.

Let r_i=u_i-kappa*sinh(u_i), Vc=sum((cosh(u)-meanCosh)^2),
Ss=sum(sinh(u)^2), R=sum(r_i^2). For a corrected log-height or log-area, add its
fixed-amplitude property gradient to g, yielding h. The full influence on each
log observation is a local term on this trace plus a shared term on EVERY trace:

    L_i = w_i-h_C*(cosh(u_i)-meanCosh)/Vc-h_D*sinh(u_i)/Ss,
    t = (h_b+kappa*h_D)/(mR),
    influence_ki = 1(k==j)*L_i + t*r_i.

For apex use its shape gradient and omit the raw-amplitude/direct w term.
Covariance is the complete dot product of these influences over all mn
observations, multiplied by sigmaHatSquared. This includes all cross-property
and cross-trace shared-b terms, not only diagonal area variance. Internally
local vectors/shared scalars avoid storing a dense (3m)x(mn) Jacobian. The
returned row-major covariance uses NATURAL properties in order
`[apexU_0,height_0,areaU_0,apexU_1,height_1,areaU_1,...]` after delta transformation.
The returned areaUncertainty is sqrt(the corresponding covariance diagonal).

Apex gradient follows implicit differentiation of b-C*sinh(x)-D*cosh(x)=0.
Corrected-height gradient uses stationarity at the apex. Area derivatives use
the closed Bessel-K integral; C,D derivatives use recurrence relations and the
b/order derivative uses central differences with Richardson extrapolation and
an explicit convergence check. Numerical derivative/Bessel failures leave
inference unavailable rather than fabricating covariance. See
[NIST DLMF integral representations](https://dlmf.nist.gov/10.32) and
[recurrence/derivative identities](https://dlmf.nist.gov/10.29).

These are first-order fixed-window reference uncertainties with estimated
log-noise variance. They exclude model mismatch, selection uncertainty, sample
correlation, calibration error and non-Gaussian/raw additive noise. Exploratory
results do not acquire calibrated coverage by returning this matrix. Width
and width uncertainty remain unavailable; no unsupported analytic FWHM is added.

## Physical coordinate/Jacobian interpretation

Only an externally supplied AFFINE map is supported:

```cpp
// origin is the physical coordinate of the selected window center.
auto physical = mapHyperbolicAffine(result, {origin, unitsPerU, AxisUnit::RetentionTime});
```

For x=origin+s*u, position=origin+s*apexU, height is unchanged and
areaX=abs(s)*areaU. Covariance transforms as D*Cov*D^T with per-trace diagonal
D=(s,1,abs(s)). Negative axes are valid; signed position cross-covariance is
preserved. Zero/nonfinite scale, unavailable inference and numeric overflow
are rejected. This map does not infer calibration, propagate calibration
uncertainty or approximate a nonlinear Jacobian with a constant. Nonlinear
instrument coordinates and irregular/missing/weighted samples need a separate
explicit model; no physical accuracy claim follows from the dimensionless area.

## Verification and performance protocol

Before edits: Debug/Release with research ON, 14/17 tests pass. After: 15/18,
with exactly the old test_qpeaks, test_qpeaks_pt03 and test_qpeaks_pt04 failures.
No expectations were changed. Build/test logs and immutable Stage B comparison
are collected by `.github/workflows/peak-model-stage-c.yml` in both configurations.

New independent checks use full dense OLS, direct raw scalar LS, apex derivative
bisection, and Simpson quadrature. Every covariance entry for m=1/3 is checked
against observational finite differences of this separate reference, including
cross-trace corrected-area covariance and all cross-property blocks. Tolerance
is .003 times the geometric mean of reference diagonal variances plus 1e-10;
this accounts for finite-difference step .002 and float log quantization. Property
values use 1e-7 absolute+relative tolerance; RSS 1e-9. Exact legacy differential
retains zero tolerance. F survival is checked against an independent incomplete-
beta continued fraction (2e-7) and a deterministic 1500-window Gaussian-null
simulation (local nominal 5% rejection rate: 5.1333%). This simulation tests the
fixed-window nominal comparison, not post-selection detection accuracy.

Tests also cover exploratory versus prespecified decisions, declared family alpha,
Bonferroni adjustment, observed support, invalid raw/step/options, unavailable
mapping, negative affine scale and signed covariance. Existing Stage B dense
OLS/shape/refit/numeric guard tests continue unchanged. Local legacy qPeaks
Release object hash is unchanged from Stage B; all 27 synthetic centroid/EIC/bin/
feature/group CSVs match, as does the complete serialized legacy probe output.

Common standalone Stage B/C benchmarks use nine alternating paired trials and
identical source/settings on the same hardware. Geometry, hot kernels, collection
and postfit are separate costs. The optional `bench_peak_inference` measures the
complete selected-window layer (49 points, m=1/4/16, nine trials x five calls),
including geometry, retained fit, Bessel/order gradient and covariance, in ns per
group. It is not a hot-window timing. CI uploads raw inference timings, paired
common summaries, full suites, compiler/CPU provenance and compile commands.
Exact final SHAs, CI links and final measured medians/spread are in the Draft PR.

No merge or pipeline/quantification approval. Remaining work before production
use includes ambiguity/interference and data-selected search calibration,
real-instrument correlated-noise validation, nonlinear coordinate/calibration
models, numerical Bessel range extension and performance of any future detector.
