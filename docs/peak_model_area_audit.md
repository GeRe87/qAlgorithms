# Area audit for matched TripleTOF windows

The input is the checksum-verified segmentation evidence at
69f0bc653252b09df9c254b33fd87884eb5ad583. A record enters this audit if its
original scan, center and half-width appear in both the positive-window
audit and the segmented search. Cohorts are retained for comparison.
The 2,027 records represent 1,985 distinct windows across both injections.
[Source CI and artifact](https://github.com/GeRe87/qAlgorithms/actions/runs/42999805753)

Each window uses the existing fixed step of 0.35 in dimensionless u. The
analysis reconstructs the scalar raw least-squares amplitude from the emitted
shape coefficients and original intensities. The observed area is the
trapezoidal sum from the first to the last sample. The fit is integrated over
the same bounds using adaptive quadrature. Both areas use raw intensity and
the same coordinate units. A second trapezoidal sum uses the fit evaluated at
the observation nodes to expose the effect of using a continuous integral.

For integrable fits the quadrature also evaluates both semi-infinite tails,
splitting at the analytic apex where applicable and normalizing by the fitted
height. The sum of the window and both tails covers the whole real axis. It
is compared with the Bessel area emitted by the native audit. The relative
comparison tolerance is 2e-8. Tests include symmetric and asymmetric
profiles, an apex outside the observed window and a constant profile with
a finite window area. [SciPy quadrature documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)

## Measured results

All 885 finite whole-line comparisons passed. The maximum relative
differences were 6.94e-14 and 4.83e-13 for injections 1 and 19. The other
1,142 records contribute observed-window comparisons. The following values
are percentages. The window difference is 100 times fit-window area divided
by observed trapezoid area minus 100. Tail fractions use the subset with
finite whole-line area.

| Injection | Cohort | Samples | Records | Window difference p10 / p50 / p90 | Finite whole areas | Tail p50 / p90 |
|---|---|---:|---:|---:|---:|---:|
| 1 | Fixed index | 9 | 235 | -4.43 / -0.12 / 2.84 | 72 | 7.59 / 28.12 |
| 1 | Fixed index | 17 | 42 | -14.96 / 2.93 / 14.97 | 16 | 0.09 / 0.29 |
| 1 | Observed maximum | 9 | 563 | -5.07 / 0.08 / 3.35 | 263 | 6.05 / 15.56 |
| 1 | Observed maximum | 17 | 93 | -18.29 / 2.70 / 9.82 | 56 | 0.33 / 7.34 |
| 19 | Fixed index | 9 | 223 | -6.50 / -0.19 / 2.56 | 55 | 6.25 / 57.15 |
| 19 | Fixed index | 17 | 48 | -15.73 / 4.13 / 12.63 | 26 | 0.50 / 4.65 |
| 19 | Observed maximum | 9 | 680 | -6.16 / 0.21 / 3.28 | 312 | 5.64 / 19.84 |
| 19 | Observed maximum | 17 | 143 | -15.34 / 2.36 / 10.23 | 85 | 0.27 / 6.23 |

The nine-point median window differences are close to zero across the four
cohorts. Some short fits extrapolate a substantial fraction of their area.
The seventeen-point groups represent different eligible windows and have
larger observed-window differences. This table does not estimate the causal
effect of window length. A paired width experiment requires holding the
anchor and observations fixed. Absolute error on the real profiles remains
unknown. The fixed-step u area requires a physical coordinate map and its
Jacobian for interpretation in mass or time units.

## Development direction

The model contract, typed hyperbolic fitting, full first-order postfit
covariance and declared-family search are implemented in stacked drafts.
Segmented evaluation now accepts all 160 audited complete scans and carries
the original indices, threshold-based exclusions and exploratory status.
The analytical and numerical area agreement supports the implementation of
the fitted-model integral. The real-data evidence measures observed-window
agreement and extrapolation sensitivity.

The proposed target is one peak-analysis framework with selectable models,
explicit input and uncertainty contracts, comparable diagnostics, and a
reviewable path from experimental fits to validated detector outputs.
The piecewise-quadratic baseline supplies the compatibility and performance
checks. The hyperbolic model can be compared in the same workflow under
declared experimental conditions.

The next experiment uses simulated profiles with known integrals. It
varies sampling, window width, background, overlap, noise and deviations
from the fitted shape. It reports bias, spread and uncertainty coverage by
condition, alongside tail fractions. A matched anchor experiment can
separately assess window-width sensitivity on real profiles. These results
can guide a decision about which diagnostic flags to expose. Statistical
calibration under selection and physical axis interpretation are separate
requirements for deployment decisions.

## Reproduction

Use the real-profile-segmentation-evidence artifact from run 42999805753.
Install requirements.txt in an isolated Python environment and run
test_area_audit.py and area_audit.py with the extracted artifact directory
and an output directory. The workflow uploads per-window CSVs and a summary
with the source commit, input hashes and grouped metrics.
