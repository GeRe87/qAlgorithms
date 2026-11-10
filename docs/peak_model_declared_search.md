# Declared-window hyperbolic search and unresolved overlap groups

## Scope and API

Base: `5807cc81661000d0d8e467ef50440cbf3c8ddc28` (Stage C, clean tree).
Branch: `agent/peak-model-declared-search`. Separate Draft PR targeting
`agent/peak-model-inference-stage-c`; no edits to PRs #1–#7 or their branches.

The new optional search layer calls the existing selected-window inference API
on a fixed, explicit family. It changes no legacy detector, CLI, output formats,
mathematical kernels, Stage B binding, or Stage C inference implementation.
New C++ follows the existing indentation, brace placement and camelCase naming.

```cpp
#include "qalgorithms_peak_search.h"
using namespace qAlgorithms::peak_model;
HyperbolicSearchOptions options;
options.selection = WindowSelection::prespecified;
options.familyAlpha = 0.05;
// Configure before observing intensities. This example declares all fitting
// scales 6..12 on a center grid anchored at 6 with stride 4.
const auto plan = makeUniformHyperbolicSearchPlan(length, fixedStep, 6, 12, 4, options);
const auto result = plan.evaluate(rawTraces);
```

Alternatively construct `HyperbolicSearchPlan(length,step,explicitWindows,options)`.
Windows are canonicalized by center, then scale; duplicates, invalid bounds,
invalid geometry/options and empty families are rejected before any raw data is
used. The plan exposes no mutation during evaluation. The constructor cannot
prove the caller did not choose the plan using observations: prespecification
remains an explicit scientific contract, not a runtime certificate.

Before fitting, evaluate checks every supplied sample in every full-length trace
for finite positive intensity and rejects empty/null traces. It returns a complete
value on success, with no output-file writes or visitor-side partial output. Raw
array lengths/alignment remain pointer contracts. There is no zero imputation,
background estimation, clipping correction or irregular-coordinate interpolation.

The result carries length, trace count, fixed step, options, actual family size,
EVERY window/inference result, candidate indices and unresolved groups. Raw
observations are not retained. All coordinate/area/covariance semantics are those
of Stage C: apex is relative dimensionless u, area intensity*u, natural-property
covariance is a fixed-window first-order reference. Selecting a window does not
turn that matrix into post-selection or simultaneous confidence intervals.

## Family accounting and interpretation

Each inference call uses the complete plan size as its family denominator,
including shape-invalid, nonintegrable, unsupported and numerically unavailable
fits. No prescreen, result filtering or best-window choice shrinks the denominator.
A nonzero totalFamilySize may declare additional tested trace pools, steps or
plans. It cannot be smaller than the plan. Repeated plans or data-driven step/pool
choices must be included in a larger genuinely declared family; separately
running many plans at alpha does not provide one experiment-wide alpha bound.

Default selection is exploratory, with no candidates/conditional decisions.
Explicit prespecified mode uses Stage C's pooled curvature comparison and shape/
support gates at the declared family alpha. Candidate means a window/pool that
passes those conditional criteria; it is not an individual trace significance,
chemical identity, true area reference, model-adequacy test or production approval.

The model remains independent homoskedastic Gaussian log errors, fixed grid and
fixed trace pool. Under valid marginal null tests, Bonferroni's union bound permits
dependence between overlapping window tests; it does not permit unmodelled serial
correlation in the original observations. See the general
[NIST Bonferroni reference](https://itl.nist.gov/div898/handbook/prc/section4/prc473.htm).
Our application to this declared pooled hyperbolic family uses Stage C's derived
numerator DF=2m and residual DF=m(n-3)-1; no legacy statistical gates are reused.

This resolves accounting over a genuinely declared family. It does NOT calibrate
an adaptively constructed family, validate a physical peak-detection search,
repair misspecified noise, or account for selection in estimated areas/covariance.

## Ambiguity is retained

Every conditionally passing window is retained. Candidate support intervals are
inclusive `[center-scale,center+scale]`. A sorted interval sweep forms connected
components, including touching endpoints and transitive overlap. Components of
size >=2 become unresolvedGroups with sample envelope and canonical window indices.
Each member points back to that group; single candidates retain a sentinel index.
No winning scale, peak merge, area sum or number of analytes is inferred.

The same physical peak may produce many overlapping window descriptions; different
physical peaks can also have connected support intervals. The groups therefore
indicate overlapping hypotheses only, NOT proof of interference or of a single
physical peak. The conservative grouping can join distinct peaks through a chain.
It cannot discover interference in rejected fits or certify unlinked candidates
as isolated. Failing windows remain visible with their inference statuses.

## Validation and performance

Before edits, full Debug and Release with research ON: 15/18 passing. After:
16/19. Exactly the known failures remain test_qpeaks, test_qpeaks_pt03 and
test_qpeaks_pt04; no existing expectations were modified. The new test includes:

- Complete family versus direct independent fixed-window API calls, including
  equality of statuses, decisions, adjusted p-values and full covariance.
- Independent all-pairs overlap graph/traversal versus interval-sweep components;
  separated groups, singleton, touching intervals and transitive chains. No
  production grouping helper is used as its own test oracle.
- Reordered requests, independently enumerated uniform grid, larger global family,
  duplicate/boundary/options/empty/null errors, and invalid data outside requested
  windows caught by the full-trace prescan.
- 800 simulated null traces with a declared family of 32 overlapping windows each.
  Reference family decisions use an independent incomplete-beta evaluator. Local
  any-window nominal false rejection: 68.75% uncorrected versus 5.0% Bonferroni;
  shaped candidate families: 0/800. Monte Carlo tolerance is .08 for nominal .05;
  these finite synthetic rates are not a guaranteed measured real-data FWER or
  a claim of zero false chemical peaks.

Local existing legacy serialized probe and all 27 centroid/EIC/bin/feature/group
CSV files match Stage C exactly. Legacy Release qPeaks object SHA256 is unchanged:
`26484C27C0BB1992EEC0D72D0E52C978C1F58432BEC1FC2C7C03793740A74FCA`.
The entire common local Release benchmark executable .text section is also
identical: `20E87E2B1158CCE577D879E7AE71C39789F4BCCF68062BA55F62A2707CB1F7E1`.
Timing variation between identical code remains possible on shared hardware.

Dedicated Debug/Release CI builds immutable Stage C separately, runs full suites,
exact numerical and 27-file CLI differentials, and the same nine paired common
benchmarks with recorded compiler flags and CPU. It also records executable code
section hashes, raw bounded-search timings and null-family test output. Both
checkouts use matched -ffile-prefix-map flags to normalize embedded source paths:
Debug assertions otherwise embed different baseline/current checkout filenames,
which can change executable layout and confound timing despite unchanged sources.
This changes CI provenance paths only; normal project compiler flags are untouched.
The common paired benchmark also copies each binary to the SAME temporary executable
path before invocation (outside timed regions), keeping argv[0] and the working
directory identical. This removes another process/stack-layout difference without
changing either benchmark program, category names or measurement counts.
Final
HEAD, CI links and performance distributions are in the Draft PR.

The optional bench_peak_search measures the whole new diagnostic plan evaluation:
49 samples, scales 6..12, center stride 4 (54 windows), m=1/4/16, nine trials. Plan
construction is excluded; sample validation, ALL window fits/properties/covariances,
collection and overlap grouping are included. Local Release medians were
1.313/2.539/7.795 ms per complete plan on Intel i7-11850H, MSYS2 GCC 15.1.0.
These are bounded diagnostic-batch costs, not running-fit ns/window.

This initial collected interface intentionally reuses Stage C without optimizing
or changing its kernel: it prepares selected-window geometry and retains complete
covariance per window. Work scales with the full family; retained memory includes
O(number of windows * traces^2) covariance storage. Large scan-wide production
searches need a separate performance/memory design, not a claim that these costs
are allocation-free. The old hot loops and their allocation/cache behavior remain
unchanged. No default CLI or automatic feature/centroid pipeline routing is added.

## Remaining before production use

Independent real-instrument validation of noise/correlation, background and missing
samples; calibrated adaptive selection/post-selection uncertainty; interference
resolution and identifiability of multiple physical components; nonlinear physical
coordinates/calibration errors; numerical Bessel range extension; scalable search
and detector-level output validation. Draft only, no merge or release.
