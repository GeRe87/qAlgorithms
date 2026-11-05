# Centered arithmetic and a common-support point estimator

This is the second experimental checkpoint. It adds an optional arithmetic
reference and an actual shared-curve point estimator on **prespecified**
supports. It does not select supports automatically, decide peak existence or
count multiple peaks. Production code, tests, default CMake targets and CSV
feature output remain unchanged.

## Centering is an arithmetic identity, not a new noise threshold

For a fixed design X whose first column is one, with M=(X'X)^(-1),

    M X' (z-a 1) + a e0 = M X' z.

Thus subtracting a common log offset and restoring the intercept leaves exact
least squares unchanged. The offset is not an estimated physical background;
it introduces no statistical parameter or additional degrees of freedom.
`centered_kernel.hpp` uses the middle log value, or the mean of the middle two
for even blocks, to preserve reflection symmetry of the arithmetic reference.
It retains the existing inverse matrices and moving moment updates, with
centered double running sums. The existing float log observations are retained
for this comparison. No near-zero curvature tolerance is used.

The original kernel and its current-cascade annotations remain available as
the default collector mode. `--centered` changes only the exploratory fit
field. **`current_accepted` still annotates acceptance of that window by the
unchanged production kernel**, not a validation test on centered coefficients.
The strict API comparison likewise always runs the old production path.
Its unchanged outputs therefore establish isolation, not drop-in equivalence
of a future production kernel replacement.

On the same 5220 inputs and 6502840 windows:

* The exact constant control changes from 2301 formal interior maxima to zero,
  with exactly zero slope, curvature and RSS. This follows from zero centered
  observations, not a tuned tolerance.
* One additional weak noise window in the existing protocol changes its
  geometric maximum classification. This is recorded, not suppressed.
* All ten reflected controls retain their classifications. Their maximum
  raw apex reflection discrepancy decreases from about 0.001086 sample to
  7.11e-13 sample. This is a numerical observation on these controls, not a
  universal bound for nearly zero curvature.
* Release and assertion-enabled -Og Debug produce five byte-identical CSVs.
* Independent arithmetic checks cover 28295 windows: all detailed synthetic
  windows and the first/last center at each scale in each detailed real block.
  The native-versus-direct-table errors satisfy a componentwise bound derived
  from gamma_k=k*eps/(1-k*eps), weighted absolute moment sums and the table.
  The maximum used fraction of this bound is about 0.271. These numerical
  checks are not peak acceptance thresholds.
* Comparison to independently constructed QR fits has maximum coefficient
  discrepancy about 4.05e-10 on those checked windows. The existing printed
  inverse-table rounding remains; it is not hidden or recalibrated here.

`observations.csv` additionally exports `log_y_float`, the exact native float
log observation used by the arithmetic audit. Nonpositive observations have a
blank log and are explicitly counted; they are never assigned an invented log.

## A shared curve on original observations

`common_profile.py` implements the log-mean family

    mu_i = b0 + b1 (i-h) + b2 (i-h)^2 1(i<h)
                           + b3 (i-h)^2 1(i>h).

The common support is fixed before fitting. Every original measurement appears
once in the objective RSS=sum_i(log(y_i)-mu_i)^2. There is no independent vote
or extra likelihood contribution for an overlapping local regression.

For each interior integer hinge h with full-rank design, the four coefficients
are solved by column-normalized QR. Operators are cached by support length and
hinge; no per-trace matrix inversion is performed. All h=1,...,n-2 are retained,
including weak one-point tails, rather than adding an arbitrary width rule.
At least five observations are required to leave positive residual degrees
of freedom for the four-coefficient reference. This is a broader final-fit
hinge family than the symmetric local windows, not a relabeling of an existing
single candidate. The cache and Python QR work are experimental reference
calculations, not a production speed or memory guarantee.

The global RSS minimum over this **discrete** hinge set is selected. Under
independent, common-variance Gaussian log errors on that fixed support this
is a maximum-likelihood point estimate. It is not an optimization over all
continuous hinges or multiple-component intensity models. Coefficients and
apex positions are reported in sample-index units, as in the local kernel.

Coefficients are unrestricted during fitting. If the winning curve has no
strict interior maximum, the result says so; a runner-up is not forced into
an accepted peak. Exact constant data have no chosen hinge or peak. Multiple
exact numeric minimizers are exposed. The complete hinge RSS profiles for
trial-zero controls are retained. Exact numeric ties do not diagnose general
statistical nonidentifiability; in particular a symmetric quadratic log curve
has no identifiable physical hinge even if floating arithmetic picks one.

The input reference uses double log observations, whereas the centered native
collector preserves the existing float log observations. This distinction is
explicit; no bitwise equality between these two separate references is claimed.
`exp(mu)` is the fitted **median** intensity under the specified zero-mean
log-Gaussian model, not the arithmetic intensity mean. No raw-height correction,
area integration, RT interpolation or signal resampling is introduced.

A stack of local estimates has covariance A Sigma A'. When the stack spans the
original-data space, the earlier correlated-ensemble GLS reference reduces to
an original-data fit. This motivates counting original observations once.
It is not asserted that every selected bundle spans that space: a single
five-point local fit retains only four coefficient combinations. The new
reference deliberately keeps the original measurements and does not infer
missing information from a covariance inverse or coefficient average.

## Complete protocol and limitations

6425 fitting cases are serialized and accounted for:

* 5220 complete input traces from the first experiment;
* 301 fixed narrow-peak supports [18,22] and 301 fixed overlap supports [24,52],
  including their reflected trial-zero examples;
* 300 asymmetric one-peak simulations and their paired reflections;
* three noiseless asymmetric amplitude-unit controls.

The extra supports use the **known synthetic layout**. They are not an
algorithm for choosing real peak boundaries. The new simulated log mean is
log(1000)+0.1*(i-20)-0.025*(i-20)^2 on the left and curvature -0.08 on the
right, n=41, with independent log-error sd 0.1 and recorded seed 20261105.
Its true apex is 20.625. Unit factors 1,0.01,10000 specify controls, not
acceptance thresholds. Reflections are paired, not another 300 independent
noise realizations.

Every case has a result/status row. 2604 cases contain nonpositive observations
on the requested full support and are explicitly unsupported by this log
reference. They are not silently shortened, interpolated, changed to epsilon,
or counted as successful zero detections. `scope_counts` includes every
scenario/support, including unsupported cases. All actual fitting observations
are saved in `common-inputs.bin`, with `supports.json` binding them to source
indices and the generated controls.

For the 300 asymmetric cases, the selected-hinge apex RMSE is about **0.04174
sample**; its bias is about +0.00106 sample. Fixing the known simulation hinge
at 20 gives RMSE about 0.02917. Searching the hinge therefore has a real
variance cost in this control. A single global quadratic, which is misspecified
for this asymmetric log mean, has bias about -3.696 samples.

For the symmetric Gaussian control, the three-coefficient global quadratic
is the correct smaller model and has RMSE about 0.00430, compared with about
0.01012 for the selected-hinge four-coefficient family. More flexible shared
models are not automatically better. These are empirical model-reference
comparisons on prespecified supports, not real-file performance guarantees.

On overlapping peaks and the full pt_04 trace, a single shared model is
visibly inadequate. Its geometric apex does not establish peak existence,
correct shape or peak count. The retained legacy expectations are not rewritten
into one-model outputs; pt_03/pt_04 remain unresolved production controls.

No covariance-derived selected-hinge confidence interval is presented as
calibrated. Both support selection and hinge search must enter such a
calibration. A future multi-peak model must add component intensities in the
raw domain; taking a sum of individual log profiles would be a different model.

## Checks and reproduction

Nine collector tests pass in Release, standard unoptimized Debug and
assertion-enabled -Og Debug. They retain the preceding six tests and add exact
constant arithmetic, field enumeration/isolation, centered reflection and axis
units. Seven profile tests check an analytic asymmetric mean, amplitude-unit
invariance, independent SVD predictions for every hinge, noisy reflection,
exact constants, geometry reuse and explicit unsupported inputs. The preceding
five structure/covariance tests still pass. No original test was weakened.

From a repository checkout with the normal qalgo_lib already built:

```text
g++ -std=c++20 -O3 -DNDEBUG -Iinclude -Iexternal/libcerf experiments/qpeaks_ensemble/collect.cpp build-release/libqalgo_lib.a -lz -o ensemble-collect.exe
python experiments/qpeaks_ensemble/test_collect.py ./ensemble-collect.exe
./ensemble-collect.exe STUDY/inputs.bin REFINEMENT/centered-release --centered
python experiments/qpeaks_ensemble/common_profile.py STUDY REFINEMENT/common-profile
python experiments/qpeaks_ensemble/test_profile.py
python experiments/qpeaks_ensemble/test_analysis.py
python experiments/qpeaks_ensemble/analyse.py REFINEMENT REFINEMENT/centered-analysis --collection-name centered-release
```

For the last command, copy the unchanged `inputs.bin` and `scenarios.json` from
STUDY into REFINEMENT. NumPy/Matplotlib requirements and `--numpy-root` are as
in the first experiment. `common_profile.py --no-plots` needs only NumPy.
Use the original Debug build flags/libraries to build the native collector
with -g (standard Debug) or -Og -g, _GLIBCXX_ASSERTIONS and -ftrapv.

`audit_refinement.py STUDY REFINEMENT` checks the complete production-comparison
accounting, Release/-Og Debug CSV identity, independent native arithmetic and
all serialized common supports/results. It expects the fresh Debug output
under REFINEMENT/centered-debug-checked. The evidence package stores one copy
of the byte-identical native CSVs plus the verified Debug hash record; rerun
the Debug collector to repeat that comparison from files. `--matrix-header`
can point to the bundled exact inverse table when auditing outside a checkout.
No new whole-file mzML run or replacement of the default production pipeline
is claimed in this research-only continuation.
