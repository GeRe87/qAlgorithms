# Paired cascade versus native-ensemble ablation

The first paired pilot finds useful common-fit localization improvements but
does **not** establish that the ensemble can replace the statistical cascade.
It also exposes a concrete limitation of the center-slope bundle construction.
The detector and legacy expectations remain unchanged.

## Frozen comparison and matched data domain

`prepare_paired.py` freezes the first 32 trials of each longer scenario, all
shorter scenarios, and complete median-/maximum-length EICs from each real
polarity. The selected workload has 280 inputs: 266 generated-reference cases,
four legacy regression expectations and ten real EICs with unknown truth.
Reflected cases remain paired controls, never independent extra trials.
All 5220 preceding inputs appear in a selection ledger; 4944 are explicitly
outside this pilot. The trial cap is a workload choice, not an algorithm rule.

`paired_native.cpp` executes the unchanged strict and collective qPeaks APIs
and exports the actual precomputed native running fits and float log values.
Failures retain status/error and discarded partial-result counts; failed
whole traces do not contribute a successful partial output list or a zero count.

Two ensemble input domains are compared without changing any fit coefficient:

* **Positive domain:** all positive blocks, as in the preceding ensemble method.
* **Matched domain:** the identical blocks above the cascade's inherited
  intensity floor exp(1). The original cascade is not retuned.

The native helper exports both domains explicitly in separate runs. An
independent full-cohort check proves every matched-domain regression is the
exact corresponding original-kernel coefficient tuple from the positive field.
All 280 cascade outcomes are unchanged between domain runs. This control is
necessary: better localization using many low-intensity tails is not evidence
for an ensemble improvement on the same information.

## What counts as an ensemble output here

Every atomic bundle is fitted, not only the earlier largest-bundle pilot.
The fit remains the same cached-geometry, prefix-moment piecewise log model.
For this explicitly **geometry-only ablation**, retain a proposal when the
unique numerical winning hinge has a strict interior local maximum whose
fitted height exceeds both fitted support edges. These are analytic shape
conditions; there is no fitted significance, distance, vote or persistence
cutoff. Merge/split proposals remain flagged as junction diagnostics, not
resolved physical components. All failed or ineligible fits remain recorded.

This is a completely specified uncalibrated ablation, not a completed ensemble
detector. A screening alarm means at least one proposal, versus at least one
native cascade description. A processing error is a third outcome, not no alarm.
Proposal multiplicity describes fragmentation/ambiguity; it is not a validated
physical peak count. Count agreement alone does not demonstrate reconstruction.

The radius-free one-to-one minimum squared-distance assignment reports
localization and unmatched counts. It cannot define precision/recall: a distant
output is still paired when it is in the smaller set. Primary localization
summaries therefore condition on the same generated inputs with the reference
number of nonjunction outputs. No oracle selects proposals or tunes the method.
Generative component centers are references, not necessarily raw-mixture modes.
The pt_03/pt_04 locations are separate legacy expectations, not new physical truth.

## Results on the fixed 32-trial base scenarios

Same floor / same native regression domain:

| Control | Cascade strict | Ensemble geometry, matched domain |
|---|---:|---:|
| Flat log noise | 31 alarms, 1 no alarm | 32 proposal alarms |
| Flat raw noise | 29 alarms, 1 no alarm, 2 errors | 32 proposal alarms |
| Zero-clipped raw noise, scale 24 | 1 alarm | 5 proposal alarms |
| Zero-clipped raw noise, scale 2400 | 1 alarm | 9 proposal alarms |
| Broad single component | one output in 32/32 | one proposal in 20/32; 2–4 otherwise |
| Narrow single component | one output in 32/32 | no proposal in 32/32 |
| Three separated components | three outputs in 32/32 | three proposals in 32/32 |
| Two overlapping components | two outputs in 32/32 | two proposals in 20/32; 1–6 otherwise, junction ambiguity retained |

These are descriptive pilot counts, not claims of calibrated error rates or
statistically established superiority. In particular, 31/32 versus 32/32 is
not evidence that either approach controls a per-trace false alarm level.

On the **same 20 broad-single-component inputs** where all three methods have
one nonjunction output/proposal, apex RMSE is 0.16484 samples for the cascade,
0.09044 for the matched ensemble and 0.01113 for the all-positive ensemble.
The large apparent unrestricted improvement includes additional tail data.
For the separated three-component case, both cascade and matched ensemble have
three nonjunction outputs on all 32 inputs: across 96 reference positions,
RMSE is 0.11165 versus 0.06117 samples. These are conditional localization
comparisons, not unconditional detector accuracy or calibrated uncertainty.

With all positive tails, the narrow case gives one proposal on 32/32 and RMSE
0.00484 versus cascade 0.04796. That is **not** a matched-domain result. At the
common floor the proposed bundling method loses every narrow case, so there
is no corresponding matched-domain localization RMSE; it is not reported as zero.

### Concrete structural counterexample

The common floor leaves five samples and exactly **one valid regression
center** in each narrow control. Its existing native piecewise model has a
strict interior maximum inside that five-sample data window. Yet its center
slope is only one value: no neighboring +,- center pair exists. The current
basin graph labels it a boundary exit and creates no maximum bundle.

Thus a boundary of the **valid-center grid** was implicitly treated like a
boundary of the **observed signal window**. These are different domains.
`overview.json` retains all 33 generated single-peak counterexamples including
the reflected control, native coefficients and their valid model apex.
This is a method limitation, not an arithmetic discrepancy or a reason to
change the reference expectations. It cannot be cured by calibrating a
significance threshold after the missing bundle has already been discarded.

The next association revision must use the stationary-point/support geometry
of the local piecewise models, including singleton/boundary-center cases,
and preserve merge ambiguity. It must be compared as a new version, with
additional trials/seeds, rather than overwriting this frozen benchmark result.

## Complete accounting, real traces and checks

The positive-domain run accounts for 1022614 native regressions, fits all
28322 bundles and records all 625565 hinge scores. The matched run accounts
for 878873 subset regressions, 28162 bundles and 611545 hinge scores. No fit
errors occur in either variant. Every best curve and every hinge RSS is
independently recomputed using an SVD-derived projection. Maximum RSS
discrepancy is about 9.44e-12; maximum best-curve discrepancy about 1.06e-10
log units. Fast-RSS subtraction differences are retained, not clamped.

Native Debug/Release field bytes and all cascade outputs agree, excluding
timing. All 28 earlier detailed native fields are bound exactly in both
domains. The ten whole real EICs (eight positive, two negative) are bound to
the fresh mzML EIC export. Real counts and runtimes are diagnostic: no true
compound/peak count, false-positive rate, recall or quantitative area accuracy
is inferred. This is a selected-EIC comparison, not a regenerated full-file
feature list. RT interpolation and area correction remain deferred.

Seven new benchmark contract tests and the eleven preceding bundle tests
pass. They verify implementation/accounting, not detector validity. The
benchmark's missed narrow controls remain explicit functional limitations.
The unchanged production pt_03/pt_04 expectations remain unresolved.

## Runtime interpretation

A dedicated final pass runs the Release native helper, then the two Python
variants sequentially with one BLAS/OpenMP thread. Native cascade timings
exclude CSV output; the Python stage includes full graph construction,
cached geometry, all fits and score logging. Cache warmness follows the
fixed case order. The native coefficient stage is timed separately.

The current Python reference is substantially slower than the C++ cascade
on this workload. These timings describe the current prototypes; they do not
compare two optimized native algorithms or establish an asymptotic speed result.
Build/debug conformance and independent-audit runtimes are not inference timing.

## Reproduction

Use the Python 3.14/NumPy 2.5.4/Matplotlib 3.11.2 environment recorded in the
local package. On Windows select Python explicitly if compiler PATH changes
would select MSYS Python instead of the matching interpreter.

```text
python experiments/qpeaks_ensemble/prepare_paired.py STUDY FRESH_EICS.bin INPUT --trials 32
g++ -std=c++20 -O3 -DNDEBUG -Iinclude -Iexternal/libcerf experiments/qpeaks_ensemble/paired_native.cpp build-release/libqalgo_lib.a -lz -o paired-native.exe
paired-native.exe INPUT/inputs.bin NATIVE_POSITIVE
paired-native.exe INPUT/inputs.bin NATIVE_MATCHED --cascade-domain
python experiments/qpeaks_ensemble/compare_bundle_cascade.py INPUT NATIVE_POSITIVE POSITIVE_RESULT --numpy-root NUMPY
python experiments/qpeaks_ensemble/compare_bundle_cascade.py INPUT NATIVE_MATCHED MATCHED_RESULT --numpy-root NUMPY
python experiments/qpeaks_ensemble/paired_overview.py INPUT NATIVE_POSITIVE NATIVE_MATCHED POSITIVE_RESULT MATCHED_RESULT OVERVIEW --numpy-root NUMPY
```

Repeat the native build with -O0 -g and the Debug library for conformance.
`audit_paired.py INPUT RELEASE DEBUG RESULT OLD_NATIVE_COLLECTION FRESH_EICS.bin AUDIT`
verifies each variant independently; add --numpy-root as needed. Complete
inputs, native fields, scores, all fits, errors and the selection ledger are
supplied locally. No measured data are committed to the repository.
