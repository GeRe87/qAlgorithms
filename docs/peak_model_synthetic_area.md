# Synthetic area recovery with paired native fits

The experiment measures area error against a known synthetic target, using the
same positive intensities and centered observed window for each model. It uses
the repository's running kernels and raw-intensity amplitude corrections. The
central window is given to each fit, so detection and selection effects are
outside this experiment.

The factorial design covers three target shapes and two sampling steps
(.175 and .35 in u). Windows have 9 or 17 points. The target is a hyperbolic
profile, a Gaussian or an asymmetric two-Gaussian mixture with total area
1,000. The Gaussian components have explicit normalized integrals. The
hyperbolic shape is normalized by independent numerical quadrature. The
background is zero or a constant 2% of the target intensity at the anchor.
The neighbor has zero area or 15% of the target area with center +.9 and
width .4. The target area excludes background and the separate neighbor.
Noise multiplies the composite signal with mean-one lognormal factors.
Scenarios use no noise, independent log noise with SD .04 or .12, or stationary
AR(1) log noise with SD .12 and lag coefficient .65. There are 40 replicates
per noisy condition and one deterministic case per noiseless condition.
The RNG seed is 74821. The paired comparison uses the same float input for
each model.

Whole-line area is reported for the piecewise model with two negative
quadratic coefficients and for the hyperbolic fit with a finite native
integral. Availability is counted before comparing errors. These fitted
areas use the full dimensionless u axis. The native models receive the
same fixed anchor and fit every condition once per replicate. Covariance
coverage is measured on each method's native scale. The legacy code
reports a log-area standard error, so its interval is symmetric in log
area. The hyperbolic postfit covariance gives a standard error in area
units. Both use a nominal 1.96 multiplier. Noiseless conditions are
excluded from coverage. Coverage is descriptive and conditional on the
availability of a positive finite standard error. The intervals use
different covariance approximations.

## Initial results

Windows attempted 5,808, with 5,446 finite whole-line estimates from both
models. The hyperbolic fit had smaller absolute error on 2,384 of those
paired windows. These counts reflect the chosen mixture of conditions.
They do not estimate a prevalence-weighted ranking on measured spectra.
All noiseless matched-shape recovery checks without background or neighbor
passed the relative tolerance 1e-5 for both sampling steps and widths.

Median absolute area errors in selected 17-point conditions with step .35,
no background or neighbor, and 40 independent log-noise draws at SD .12:

| Target shape | Piecewise quadratic | Hyperbolic |
|---|---:|---:|
| Hyperbolic | 14.81% | 4.39% |
| Gaussian | 3.44% | 9.25% |
| Asymmetric mixture | 11.01% | 6.40% |

For the matched hyperbolic target at step .35 and 17 samples, the nominal
95% interval covers 37/40 cases under independent log noise and 22/40
under AR(1) noise at the same marginal SD. This is a finite Monte Carlo
comparison of two declared noise conditions. It motivates testing covariance
against a measured error model before relying on these intervals.

The next table groups noisy fits with finite area by their fitted share
outside the observed window. The outcome is absolute error above 10% of
the target area. The bins are descriptive, with all nuisance conditions
pooled.

| Fitted tail fraction | Piecewise, count | Piecewise, error >10% | Hyperbolic, count | Hyperbolic, error >10% |
|---|---:|---:|---:|---:|
| Below 5% | 4,264 | 61.0% | 4,307 | 63.9% |
| 5% to below 20% | 909 | 24.3% | 836 | 24.2% |
| 20% or more | 581 | 74.9% | 255 | 72.2% |

A small tail fraction occurs with large area errors under these conditions.
It measures extrapolation, while contamination and shape mismatch can
also shift the integral within the measured support. The pooled rates
mix shapes, widths, noise and contaminants. For a diagnostic rule, the
next comparison should condition on those factors and evaluate held-out
scenarios. The finite-fit requirement also creates different subsets for
the two models. Neither the pooled rates nor this experiment establish
error prevalence or area coverage in TripleTOF measurements.

## Reproduction

Configure CMake with QALGORITHMS_ENABLE_SYNTHETIC_AREA_AUDIT=ON and build
synthetic_peak_area. Install the pinned dependencies from the matched
area audit, then run simulate.py with the native executable and output
directory. The script checks matched-shape recovery, window/whole area
consistency and paired row counts. The workflow publishes per-model trial
rows, all condition summaries, diagnostic counts and the tested commit.
