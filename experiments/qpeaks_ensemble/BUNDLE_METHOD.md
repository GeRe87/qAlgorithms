# Two-stage method proposal based on the native running regressions

This checkpoint addresses the intended method: native running fits form an
ensemble **before local statistical validation**, and the ensemble determines
a common piecewise fit. It introduces an executable geometry reference, not
a calibrated replacement detector. Production and old tests are unchanged.

## Model and retained information

On the equidistant measured-point index, every existing local fit is

    q(c,s;i) = b0 + b1 (i-c) + b2 (i-c)^2 1(i<c)
                                  + b3 (i-c)^2 1(i>c),
    i = c-s,...,c+s; z_i = float(log(y_i)).

The existing inverse tables, scales and native coefficients are retained.
The hinge c is not generally the apex. In particular, q'(c)=b1 exactly for
this C1 model. Each original coefficient, support and source observation stays
available. Nonpositive points, short blocks and endpoints retain the native
collector's explicit accounting. No epsilon log, interpolation or RT regularization
is introduced. The code never reads current_accepted to form an ensemble.

## Stage 1: basins in the center-slope field and a scale DAG

At each scale, define d_s(c)=b1(c,s) over its consecutive valid centers.
Follow the sign of this estimated slope toward increasing signal: + goes
one enumerated center right, - one left. A neighboring +,- pair is a maximum
**corridor [c,c+1]**, not a newly interpolated measurement or a precise apex.
Every center flowing to that corridor belongs to its attraction basin.

Exact zero runs have explicit states. +,0...0,- is a stationary maximum
plateau; equal flanking signs continue in that direction; -,0...0,+ is a
stationary minimum; zeros touching a missing boundary stay unresolved.
All-zero data are flat/boundary-stationary, not a maximum. Slopes leaving the
valid center domain are recorded as boundary exits. Every native regression
belongs to exactly one recorded basin, including nonmaximum basins.

Between consecutive native scales, link every pair of maximum basins sharing
at least one valid center. This transports the same source-center coordinates
through scale and tolerates diagonal movement of the root corridors. Comparing
centers avoids joining models merely because their broad data windows overlap.
The shared-center count is bookkeeping, not an independent vote or probability.

The resulting graph is directed and acyclic because scale increases along
every edge. Contract an edge only when its source has one outgoing neighbor
and its target one incoming neighbor. A contracted path is an **atomic bundle**.
Retain every merge, split, birth and termination. In particular, do not turn
the entire connected graph into a single physical peak. A coarser branch
after a merge may describe an unresolved mixture and is labelled accordingly.

This is an explicit geometric proposal, not a theorem that basin overlap means
the same physical peak. It has no amplitude, apex-distance, persistence or
minimum-member cutoff. These properties alone do not prove it optimal.
The native center-slope field is not claimed to be the derivative of one
globally fitted smooth curve. Its local-polynomial scales are not a Gaussian
semigroup; births and splits cannot be ruled out by Gaussian scale-space theory.
Topology and physical interpretation must therefore remain separate.

Structural precedents are the derivative/location/scale view in
[Chaudhuri and Marron's SiZer](https://doi.org/10.1080/01621459.1999.10474186)
and [Lindeberg's scale-space axioms](https://www.csc.kth.se/~tony/abstracts/CVAP159.html).
The reference adopts neither SiZer's significance procedure nor Gaussian
scale-space guarantees for this different kernel.

## Stage 2: a shared piecewise fit determined by the bundle

For one bundle E, use the union S_E of its native fit supports. Retain the
bundle's distinct native centers as the discrete hinge proposal set H_E.
Both come from the unfiltered running fits, not old accepted feature windows,
known synthetic apex locations or manually chosen boundaries. A disconnected
support is an explicit error; it must not be silently bridged.

For each h in H_E, fit the same four-column model on the **same** S_E:

    X_h(i,:) = [1, i-h, (i-h)^2 1(i<h), (i-h)^2 1(i>h)],
    b_hat(h) = (X_h' X_h)^(-1) X_h' z,
    h_hat = argmin_h sum_{i in S_E} (z_i - X_h(i,:) b_hat(h))^2.

Do not compare unnormalized RSS values from different supports to choose a
winning physical bundle. The h comparison above is valid as a point-estimation
objective because its support is held fixed. It is a conditional least-squares
estimate; under fixed support and iid Gaussian log errors it is also the
discrete-hinge likelihood estimate. Real noise and selection are not thereby
calibrated. Coefficients remain unrestricted, matching the existing model
family; a fitted curve without an interior maximum remains such a curve.
Do not substitute a shape-conforming runner-up and call it an accepted peak.

For a right-side maximum, apex=h-b1/(2*b3); on the left use b2. The root must
lie on the corresponding side and inside the support. Hinge, apex, fitted
log height, boundary dominance and possible equal numerical minimizers are
separate outputs. A merge-containing bundle is a collective diagnostic fit,
not evidence that one component explains its finer branches. No area is fitted.

### Why the coefficient ensembles are not independent data

Write a local linear fit as beta_j=A_j z on the original observations and
stack them as beta=A z. Their joint covariance is sigma^2 A A', not block
diagonal: overlapping fits reuse the same noise. Their common-model means are
A_j X_h b, not coefficients translated and averaged as if every hinge were
the same. With a known general noise covariance use A Sigma A' instead.

For iid errors and the ideal linear operators, the joint generalized
least-squares objective is proportional to

    (z-X_h b)' P_A (z-X_h b),
    P_A = A' (A A')^+ A.

If A spans the original-data space, P_A=I and this is exactly the original-data
least-squares fit. If it does not, the coefficient stack has discarded some
information. The reference deliberately retains the original observations;
it does not invent missing information or claim every bundle has full rank.
Finite native/table rounding must also be distinguished from this exact-linear
identity. The ensemble contributes structural/support information; repeating
the same measurements does not create additional statistical precision.

### Preserving the precomputation strategy

There is no nonlinear solver in this shared single-profile step. For fixed
(support length, hinge offset), X_h' X_h is data-independent and cacheable.
Using prefixes of z, i*z and i^2*z gives its right-hand side in constant work:

    p0 = sum z,
    p1 = sum i*z - h sum z,
    pL = sum_{i<h} i^2*z - 2h sum_{i<h} i*z + h^2 sum_{i<h} z,
    pR = sum (i-h)^2*z - pL.

Hence b=M_(N,h) [p0,p1,pL,pR]'. The existing symmetric local tables are
unchanged. An asymmetric union needs an **additional** geometry cache; it
cannot generally be evaluated with an unchanged symmetric table alone.
This extends the same equidistant, data-independent precomputation principle.

`MomentGeometry` builds a column-scaled 4x4 inverse from QR once per geometry
and uses centered prefix moments for its coefficients. It explicitly scores
residual curves to audit the fast identity RSS=z'z-b'p. Subtractive fast-RSS
cancellation is reported, not clamped. This Python reference is not a native
throughput implementation or a production memory/speed guarantee.

## Checked pilot and what it establishes

`study_bundle_method.py` reads **all 28 detailed trial-zero native inputs**:
524349 local regressions, including controls, pt_03/pt_04, paired reflections
and six whole real positive EICs. Other random trials in the older collection
are outside this declared topology pilot. Neither old local acceptance nor
old feature windows enter grouping or fitting.

The geometry is run on both the original native kernel field and the preceding
centered-arithmetic native reference (same model and inverse table). All basin,
edge, bundle, hinge-fit and fitted-curve outputs are byte-identical between
these two fields on this dataset. This does not prove universal sign stability.
For the earlier original collection without serialized float logs, native log
observations from the centered collection are supplied separately; all original
x, y and measured-point counts are checked one-to-one before fitting.

Every detailed regression is copied into native-membership.csv. Every basin,
edge and bundle is retained. The fit workload is explicitly **all bundles tied
for the largest member count per input**, not acceptance of those bundles.
All other bundles are marked not_fitted_in_this_pilot. There are 27 conditional
shared fits and no fit errors; the exact constant has no maximum bundle.
No physical peak count is released. This workload sample is not a final
scale/bundle selection policy and its fits do not stand for every bundle.

The overlapping-pair pilot explicitly illustrates why a merge is collective:
forcing one log-piecewise curve across its merged support overshoots the raw
height severely. It is labelled collective_diagnostic_fit, not a fitted physical
component. Likewise, a local interior maximum need not dominate both support
edges (the pt_04 pilot does not). Both failures remain visible. A symmetric
quadratic log mean has an unidentifiable hinge even when arithmetic chooses
one numerical winner; its fitted apex can still be well determined.

Every searched hinge is independently checked against direct normalized SVD.
Maximum prediction discrepancy is about 6.29e-11 log units; maximum fast-RSS
versus explicit-RSS discrepancy is about 1.37e-9. The independent audit checks
all copied coefficients, all ten reflected basin/edge graphs, exact arithmetic-
reference output identity and all six real vectors against the fresh mzML EIC
export. Eleven structure/algebra tests cover corridors, plateaus, flat input,
boundary accounting, reflection, diagonal motion, a merge, support construction,
prefix/SVD agreement, log-unit offsets, apex versus hinge and invalid inputs.

In this pilot a broad single Gaussian has one atomic bundle, the exact constant
none, pt_04 three, but the overlapping-pair control has nine and pt_03 168.
Those are graph objects, not detections. The overlap count already shows that
atomic paths cannot simply be counted as physical peaks. Substantial mode
structure persists in noise. The old pt_03/pt_04 expectations remain unchanged
and unresolved as statistical detection controls.

## Remaining decisions before a detector can be replaced

1. Check the basin-transport proposal against alternative associations and
   controlled overlapping/asymmetric peaks, especially rapid root motion,
   disappearing valid-center domains and mixed windows near junctions.
2. Establish a bundle/support selection criterion from the full data-dependent
   procedure. Member count is not evidence and pilot selection is not a rule.
3. Quantify model adequacy/background and unresolved mixtures. A shared
   single-profile diagnostic cannot deconvolve arbitrary raw intensity sums.
4. Certify numerically ambiguous slope signs with floating-error bounds.
   The current exact floating signs introduce no tuned epsilon but are not a
   general certificate. Near-zero sign sensitivity must remain visible.
5. Calibrate any eventual existence/count/uncertainty decision over the **whole**
   grouping/support/hinge selection under a stated noise model. Conditional
   inverse-matrix standard errors do not have post-selection coverage.

These are method-development tasks, not permission to add arbitrary distance,
RSS, persistence or vote cutoffs. Statistical calibration cannot be removed
merely by replacing a per-fit cascade with an ensemble. No production detector,
area correction, reference expectation or interpolation policy changes here.

## Reproduction

```text
python experiments/qpeaks_ensemble/test_bundle_method.py --numpy-root NUMPY -v
python experiments/qpeaks_ensemble/study_bundle_method.py NATIVE_COLLECTION PRODUCTION_RESULT --observations CENTERED_COLLECTION/observations.csv --numpy-root NUMPY
python experiments/qpeaks_ensemble/study_bundle_method.py CENTERED_COLLECTION CENTERED_RESULT --numpy-root NUMPY
python experiments/qpeaks_ensemble/audit_bundle_method.py NATIVE_COLLECTION CENTERED_COLLECTION PRODUCTION_RESULT CENTERED_RESULT SCENARIOS.json FRESH_EICS.bin AUDIT
```

Use one BLAS/OpenMP thread for the checked small solves. The native collections
were previously verified in Debug/Release; this checkpoint executes Python
geometry and fitting, not another native build or full mzML feature run.
