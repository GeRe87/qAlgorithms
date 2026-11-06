# A conservative singleton base case for native regression bundles

This checkpoint repairs a specific representational gap in the experimental
ensemble method. It does not replace the statistical detector, change the
native running regressions or optimize/port the method to C++.

## The base case follows from the model, not a fitted numeric cutoff

Let F_B contain every available native running regression in one contiguous
input block B, before local statistical acceptance. The center-slope field
cannot contain a neighboring +,- crossing if F_B contains only one model.
Nevertheless that existing C1 piecewise model may have an observed interior
maximum. A singleton ensemble must not be ruled out merely for lacking another
regression center.

For its existing coefficients and half-window s, compute

    q'(-s) = b1 - 2*s*b2,
    q'(+s) = b1 + 2*s*b3.

If the first is positive and the second negative, the continuous piecewise-
linear derivative encloses exactly one strict interior maximum. Its location
is -b1/(2*b2) on the left, -b1/(2*b3) on the right, or zero at the hinge when
b1=0 and both curvatures are negative. The root must lie in the observed native
window; no root is extrapolated beyond measured data.

Exactly constant **native float log observations** cannot support a peak shape.
They are detected by exact equality, not by an amplitude/curvature epsilon.
This protects the base case from the known tiny curvature artifacts of the
printed inverse table/native arithmetic. Nonfinite or missing observations
remain explicit errors. Near-zero nonconstant signs are not numerically
certified; exact equality is not a substitute for a future error-bound analysis.

`model_singleton_v3` adds this model-supported node only when |F_B|=1 for the
**entire input block**, not when one scale row has one center and not when just
one regression passes the old cascade. Every multi-model block follows frozen
v1 unchanged. The shared fit is still the cached four-parameter linear model
on original native log observations, counted once.

Under the current native minscale this base case has five observations and
four coefficients, but the implementation tests the available model-family
cardinality, not a hand-picked width or signal label. There is no minimum vote
count, tuned distance, significance level or intensity cutoff in this addition.
The inherited cascade intensity floor is used only by the matched benchmark.

For one window the common fit has no extra information beyond that window.
Its apex precision should agree with the original regression, apart from native
versus double-reference rounding. Multiple overlapping regressions are not
treated as independent confirmations.

## Why the broader boundary extension was retained only as an experiment

`model_boundary_v2` also closes general center-grid exits using an endpoint
model whose maximum is inside its observed window but outside the valid center
domain. Its geometry is analytic, but adding these nodes to the existing
basin-transport graph does not establish physical identity of neighboring
branches. On development trials 32–63, the two-component overlap control had
two proposals in 23/32 v1 cases but only 3/32 v2 cases; extra competing proposals
were introduced. Count equality itself still does not establish deconvolution.

Accordingly this broader rule was not selected as the preferred correction.
Its complete case results and frozen source snapshot are retained. V3 fixes
only the singleton base case while preserving all existing multi-model bundles,
supports, fit coefficients and hinge scores. This is a scoped correction,
not a claim to solve general boundary association.

A remaining explicit counterexample has six observed points, two available
local fits and a peak near the observed left edge. V1/v3 can still have no
center-slope crossing, while v2 produces a model-supported boundary node.
Resolving such cases without the observed overlap fragmentation requires a
better multi-model association rule. It is not hidden by a new numeric filter.

## Development and fresh validation are separated

* Trials 0–31: earlier frozen paired benchmark, 256 randomized inputs plus
  24 fixed diagnostics, including ten real EICs.
* Trials 32–63: 256 additional randomized inputs used to evaluate v2. Those
  findings informed the decision to prefer the conservative singleton scope.
* Trials 64–95: 256 further randomized inputs used to validate the frozen v3
  rule. The same 24 fixed diagnostics are reused and explicitly identified;
  they are not 24 new independent observations.

Source hashes were recorded before each additional evaluation. V2 source is
archived unchanged; v3's core hashes match the final validated implementation.
No numerical threshold was optimized. Method selection used development
results; this is not a claim that no earlier data informed the design.

The benchmark continues to fit **every** bundle and retains every ineligible
fit/error. Ensemble results are geometrically complete proposals, not
statistically confirmed physical peaks. The old cascade is unchanged in both
strict and explicitly collective modes. Primary comparisons use exactly the
same native windows above the inherited cascade floor.

## Main v3 results

On both reference trials 0–31 and validation trials 64–95, the narrow single-
component control changes from 0/32 proposals in v1 to one proposal in 32/32
v3 cases. On the fresh validation cases, cascade versus v3 apex RMSE is
0.04308225 versus 0.04308227 samples, as expected for identical singleton data.

All 28162 reference v1 bundle fits and all 28169 validation v1 bundle fits are
preserved exactly. Their 611545 and 610761 hinge scores, respectively, are
exact subsets of the independently audited v3 scores after removing bundle-ID
renumbering. V3 adds 43 reference and 38 validation singleton fits.

On validation trials 64–95, broad single peaks have conditional apex RMSE
0.24947 for the cascade and 0.05850 for the ensemble on the same 16 inputs
where both produce one nonjunction output/proposal. This does not describe
unconditional recognition: v1/v3 produce multiple proposals on the other
16 inputs. For three separated components, both have the reference count on
all 32 cases; across 96 positions RMSE is 0.12352 versus 0.05315 samples.
These multi-model estimators are unchanged from v1, not new gains caused by
the singleton patch.

Noise remains unresolved. On the fresh log-noise control, strict cascade has
30 alarms, one no-alarm and one processing error; v1/v3 have proposal alarms
on all 32 cases. On zero-clipped noise, v3 increases proposal-alarm cases from
1 to 3 (scale 24) and from 5 to 7 (scale 2400), while the cascade returns no
outputs on these 32-case groups. These additional hypotheses are visible, not
accepted peaks or evidence of calibrated error control.

Overlap remains unresolved: v1 and v3 have identical proposal multiplicities
on the fresh overlap cases, with two proposals in 17/32 cases versus two
cascade descriptions in 32/32. Junction flags remain; matching the reference
count does not prove component separation. The inherited pt_03/pt_04 tests
and their failures are not redefined or weakened.

## Verification and remaining limits

V3 refits all 28205 reference and 28207 validation bundles without fit errors.
All 1222387 hinge RSS values and 56412 best curves are independently checked
by SVD-derived projections. Maximum RSS discrepancy is about 7.96e-13 and
best-curve discrepancy about 1.71e-13 log units.

The structural audit verifies a unique owner for every native fit, retention
of all v1 maximum nodes/edges and exact preservation of all multi-model bundle
members. It checks reflection of all 560 supplied coefficient fields and the
actual native reflected controls separately. Native Debug/Release field banks
and cascade outputs agree except timing; unchanged exporters were rerun on
the additional data, not reimplemented or optimized.

All ten real EICs remain fixed diagnostics from the same fresh mzML export.
No real recall, false-positive rate or quantitative area accuracy is inferred.
There is no full-file ensemble feature export, area correction or RT resampling.
The existing complete production CSV lists remain unchanged reference outputs.

Twelve analytic boundary tests, five singleton-scope tests, eleven unchanged
v1 tests and seven benchmark contract tests pass. Those checks establish
implementation properties; they do not make the geometry a validated detector.
The original v1 case results, fitted bundles and hinge-score file reproduce
byte-for-byte. All known method limitations remain explicit.

Next work is multi-model association (including genuine boundary cases and
alternative branches), shared-model adequacy/overlap and a statistically
calibrated **ensemble-wide** decision. Neither member count nor the number
of verified fits is an independent sample size. Formal numerical sign
certificates, selected-model uncertainty and the later native C++ optimization
remain open.

## Run the reference

```text
python experiments/qpeaks_ensemble/test_singleton_method.py
python experiments/qpeaks_ensemble/prepare_paired.py STUDY FRESH_EICS.bin NEW_INPUT --trials 32 --trial-start 64
paired-native.exe NEW_INPUT/inputs.bin NEW_NATIVE --cascade-domain
python experiments/qpeaks_ensemble/compare_bundle_cascade.py NEW_INPUT NEW_NATIVE RESULT --topology model_singleton_v3 --numpy-root NUMPY
python experiments/qpeaks_ensemble/audit_model_boundaries.py NEW_NATIVE STRUCTURE --scope singleton_blocks
```

Native exporter code and build instructions are in [PAIRED_COMPARISON.md](PAIRED_COMPARISON.md).
The default topology remains v1 for reproducible comparisons. V2 and v3 are
explicit opt-in research variants; none changes the production detector.

## Reuse a native regression

The optional `--estimator native-representative` keeps the v3 grouping
and chooses one existing native regression per bundle. It ranks candidates by
log RSS on the same observed points, the union of their native windows. This
can extrapolate during ranking. The output records how many points were
outside the winner's native window, and the apex and edge checks use only
that window. An exact score tie remains an explicit ambiguous outcome.
The batched scorer evaluates polynomials directly. It creates no inverse
matrix and solves no additional regression.

On the 280 validation cases this evaluates 28,207 bundles, with 8,340
geometric proposals and no scoring errors. Of the winning models, 28,169
were scored outside their native window. An independent scalar sum of the
squared residuals agreed on all selected candidate IDs and outcomes.
Debug and Release banks produce the same results apart from timings.
The original refit outputs reproduce byte-for-byte.

On the same 16 wide-peak cases with one nonjunction proposal from both
ensemble variants, apex RMSE in scans is 0.06977 for native reuse, 0.05850
for the refit and 0.24947 for the strict cascade. All three methods agree
at 0.04308 scans on the 32 narrow cases. Exactly three proposals occur
in 26/32 triple cases with native reuse and in 32/32 with refit and
cascade. Exactly two occur in 14/32 overlap cases with native reuse,
17/32 with refit and 32/32 with cascade. These counts do not establish
correct component separation. The ten real EICs have no reference peak
labels. This remains an experimental proposal generator without statistical
acceptance or quantitative areas.
