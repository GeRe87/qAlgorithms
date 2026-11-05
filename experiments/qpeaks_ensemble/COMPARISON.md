# One versus two additive components

This is a reference comparison, not a physical peak-count decision. Production
sources, existing tests, default builds and feature CSVs remain unchanged.

## Matched model families

Both models use the same concave piecewise-quadratic log-profile family:

    u=(i-h)/(n-1)
    q=b0+b1*u-cL*u^2(left)-cR*u^2(right), cL,cR>=0
    M1: mu=q
    M2: mu=logaddexp(q1,q2)
    log(y)=mu+epsilon, iid epsilon~N(0,sigma^2)

Noise applies to the total intensity. Original observations enter once. The
components add in the raw domain, not by adding their log profiles. Exp(mu)
is the fitted median intensity under this observation model.

Concavity is a declared reference assumption, narrower than the preceding
unrestricted profile. It is not a theorem about all instrument peaks. Its
flat/log-linear boundary includes background-like functions. Components can
have apices outside the observed support. There are no minimum separation,
amplitude, width or fit-vote thresholds.

The first study fixes both hinges at the support midpoint, equally in M1 and
M2. Supports are prespecified; the synthetic overlap support uses the known
layout. This does not solve automatic peak boundaries. Hinge-search sensitivity
is evaluated separately below. Areas, interpolation, resampling and production
height corrections are unchanged.

## Optimization and nesting

Normalized geometry is cached. M1 is convex bound-constrained linear least
squares, solved by BVLS. Independent enumeration of all four active curvature
bound sets with SVD free solves verifies all 149 single-hinge references.
The largest RSS difference is about 3.78e-15; maximum checked KKT residual is
about 1.76e-13.

For M2 the analytic residual Jacobian is

    J=[w1*G1, w2*G2], wj=exp(qj)/(exp(q1)+exp(q2)).

The implementation uses stable logaddexp and bounded nonlinear least squares
with ordinary quadratic loss. Starts are a balanced duplicate of M1 plus
Gaussian fits on every contiguous split with at least three points per half.
Three points are the initializer's rank requirement, not an accepted peak
width. Every final objective uses the full support.

All starts, finite unsuccessful candidates and errors are retained. Best-found
does not mean a certified global nonlinear optimum. A tighter refinement is
recorded separately. The exact nested M1 representation is always a candidate;
constant data have an analytic zero-RSS lower bound and need no optimizer to
invent a decomposition. Component IDs retain parameter-block order across
files; curve component1/component2 correspond to IDs 0/1.

Numerical tolerances derive from float64 epsilon: sqrt(eps), then eps^(2/3).
Budgets are 800 evaluations per start and 3200 for refinement. These are
computational limits, not detection rules. The n>8 requirement defines this
unsaturated eight-parameter reference with unknown variance; it is not a
universal statement about other restricted models.

## Why ordinary chi-square/F release is not justified

Two identical half-amplitude components exactly reproduce M1. At this null
representation J=[G/2,G/2] has rank four, not eight: decomposition is not
identifiable. Vanishing components and active curvature bounds add further
nonregular cases. Full numerical rank at a fitted alternative does not repair
regularity of the null model.

For positive residual sums and unknown common variance the reported quantity is

    D_best=n*log(RSS1/RSS2_best_found).

Known-sigma controls also report (RSS1-RSS2)/sigma^2. Neither is converted to a
p-value or release cutoff. Exact zero RSS receives no invented finite ratio.
Calibration must include the actual support, hinge, branch and optimization
procedure under an explicit null/noise model.

The nonregular component-test issue has a primary-method precedent in
[Protassov et al.](https://arxiv.org/abs/astro-ph/0201547). The duplicate-Jacobian
argument above is derived for this model. Solver definitions and local/global
distinctions are documented in
[least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
and [lsq_linear](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.lsq_linear.html).
These references do not calibrate the present score automatically.

## Complete first protocol

The ledger accounts for all 6425 preceding cases: 105 selected, 3415 outside
the first-trials protocol, 2604 nonpositive log inputs, and 301 with n<=8.
Four generated controls bring this comparison to 109 cases. This is not a
claim that all 6425 cases or the real measurement were compared.

The main sets are the first 32 predefined trials each of one Gaussian, one
asymmetric peak, and an overlapping pair. Other inputs include unchanged
pt_03/pt_04, constants, amplitude controls and paired reflections. New controls
are an exact overlap, a close pair, reflection of the same sliced overlap and
one peak plus constant background. Truth is used only for evaluation, never
fit starts. Background is a deliberate misspecification example, not an extra
valid null replicate. Counts and seeds specify experiments, not release rules.

All 109 cases have results. Best candidates converge or are analytically
supplied. All 5156 candidate/start records remain, including 82 budget-limited
attempts. Global nonlinear optimality is not certified.

| Main group, 32 each | Median RSS1 | Median RSS2 | Median D_best |
| --- | ---: | ---: | ---: |
| One Gaussian | 0.73840 | 0.69452 | 7.91465 |
| One asymmetric peak | 0.39995 | 0.33137 | 6.48583 |
| Overlapping pair | 10.00827 | 0.05160 | 152.48033 |

RSS is compared within matched supports. Groups have different lengths/noise.
Overlap component-apex RMSE is about 0.07559 sample over 64 errors from 32
traces, not 64 independent trials. No coverage or real-file accuracy follows.
Two interior component apices occur in 31/32 Gaussian single-peak fits and all
32 asymmetric single-peak fits. These are not released false-positive counts;
there is no release rule. Geometry alone cannot establish two physical peaks.

One peak plus background gives RSS1 about 6.84297 and RSS2 about 0.41147,
D_best about 115.26. The second function absorbs background-like structure and
even obtains a broad interior apex under noise. A large gain can therefore
indicate missing background/shape terms. Constants have zero RSS in both
models, rank-four duplicate components and no unique apex.

## Hinge sensitivity

Overlap trial zero and pt_04 include all interior integer hinges plus the fixed
midpoint reference. Every unordered pair including duplicates is covered:
378 overlap pairs and 91 pt_04 pairs, 469 total. Warm projected components in
both label assignments and a geometric midpoint Gaussian split give finite
local starts, not a global certificate. There are 1836 pair-level candidate
records plus two refinements; 14 pair-level attempts are unsuccessful.

| Example | Fixed RSS1 / RSS2 | Searched RSS1 / RSS2 | D_best fixed / searched |
| --- | --- | --- | --- |
| Overlap | 9.85493 / 0.04565 | 7.24155 / 0.04213 | 155.869 / 149.259 |
| pt_04 | 27.65404 / 8.96670 | 27.28971 / 2.23412 | 15.768 / 35.037 |

Selection can materially change the score. A one-versus-two study does not
validate the original three-peak expectation; neither legacy test is changed.

## Verification and reproduction

Fourteen new tests cover raw addition, label symmetry, analytic Jacobians at
equal/different hinges, initializer translation, duplicate rank, convex KKT,
independent active-set solutions, exact two-Gaussian recovery, constants, units,
component-file identity and explicit failed/unsupported cases. The preceding
seven profile and five structure tests pass. No old expectation is relaxed.

The independent audit recomputes every finite objective (5156 conditional and
1838 grid records), all 469 pairs, 149 convex null profiles by another algorithm,
input selection, generated controls and parameter/curve identity. Maximum
conditional RSS discrepancy is about 1.14e-13. SVD rank uses machine precision
as a numerical diagnostic, not a physical signal threshold.

Optional dependencies used: Python 3.14, NumPy 2.5.4, SciPy 1.18.1, Matplotlib
3.11.2. No default C++ dependency changes. From the repository root:

    python experiments/qpeaks_ensemble/test_components.py
    python experiments/qpeaks_ensemble/compare_components.py COMMON COMPARISON
    python experiments/qpeaks_ensemble/hinge_sensitivity.py COMPARISON GRID
    python experiments/qpeaks_ensemble/audit_components.py COMMON COMPARISON GRID AUDIT

COMMON requires the preceding common-inputs.bin and supports.json, supplied in
the local package. Local installations can use --numpy-root and --scipy-root.
The main comparison supports --no-plots. --trials is the explicit computation
budget; --pilot marks a three-case debugging subset, not the full study.
The full run uses one BLAS/OpenMP thread for small dense solves.

Remaining method work includes the physical mean family/background, automatic
support/branch selection, optimizer reliability and calibration of the complete
comparison. An RSS ratio, minimum distance, vote count or nominal chi-square
cutoff is not substituted for that work.
