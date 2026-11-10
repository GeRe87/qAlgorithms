# Real TripleTOF profile-mode evidence: paired log-hyperbolic vs qPeaks

**Research status: no production acceptance.** This report records two real profile-mode mzML runs from [Zenodo 20729183](https://doi.org/10.5281/zenodo.20729183), analyzed on branch `research/hyperbolic-running-regression`, [GitHub Actions run 42831239113](https://github.com/GeRe87/qAlgorithms/actions/runs/42831239113).

## Provenance

- Source: *Profile-mode LC-MS data files from human serum sample pools*, public TripleTOF 5600+ QTOF profile-mode MS1 spectra, positive polarity.
- Injection 1: `20171016_POOL_POS_1_105-134.mzML`, MD5 `4130309a9283c8b492660dbad1fd3670` (validated during CI).
- Injection 19: `20171016_POOL_POS_3_105-134.mzML`, MD5 `8bce27a1fdbe2ea0f0659d398b0bfaaf` (validated during CI).
- Both are distinct injections from the **same analytical LC-MS run**, not independent calibration standards or different instruments.
- Input SHA-256 checksums, per-profile fit CSV, run logs and parsed JSON summaries are stored in CI artifact `hyperbolic-real-tripletof-evidence` from the run above.

## Predeclared protocol

- Native `qAlgorithms::XML_File` and `get_spectrum` decode raw profile mzML arrays. No centroiding, no intensity interpolation, and no synthetic missing-point filling.
- Sample 80 MS1 spectra evenly over the full run for each injection; no scan selection based on which model wins.
- Model-blind candidate selection requires an observed local maximum `>=100`, clean shoulders, positive observed samples, uninterrupted reasonably regular sampling in `sqrt(m/z)`; local native-step CV <=20%.
- Evaluate observed five-, seven- and nine-point windows: half-widths `h=2,3,4`; no larger window fabricated.
- Group **consecutive m/z-sorted candidate profiles within the same scan** in nonoverlapping batches of `m=1,2,4` traces. These are descriptive statistical groups only; they are **not identified as the same chemical species**.
- Compare the actual `qAlgorithms::forEachCoefficientsMultiTrace` and the isolated `hyperbolic_research::forEachRunning` on identical arrays, and perform separate *closed-form raw-intensity amplitude-only corrections* for both methods.
- Fixed dimensionless regression-axis step `u=0.35*(index-center)`. This is an *externally fixed exploration setting*, not a calibrated native m/z or TOF coordinate transform. The observed native sampling is only screened for regularity; we do **not** claim exact physical coordinate equivalence.
- Evaluate log RSS, AICc (with parameter-count penalty when defined), raw RSS before and after amplitude refit, integrability flags, apex offsets, and **full model-based area in u**. No independent reference area or exact mass is supplied.
- Timed kernels exclude mzML decoding, profile selection, and postfit special functions. Times are hardware/compiler/run dependent.

## Source data counts

| Metric | Injection 1 | Injection 19 |
|---|---:|---:|
| MS1 spectra decoded | 80 | 80 |
| Screened observed peak windows across all widths | 1,225 | 1,170 |
| Paired model groups across widths and trace counts | 1,985 | 1,904 |
| Evaluated per-trace fits (repeated across widths/groups) | 3,161 | 3,046 |
| Valid raw-intensity amplitude refits | 3,161 | 3,046 |

The counts are **not independent physical peaks**: one observed peak may appear in more than one width or pooling configuration.

## Matched results: nine observed samples, four pooled traces

| Metric | Injection 1 | Injection 19 |
|---|---:|---:|
| Evaluated pooled groups | 50 | 52 |
| AICc prefers log-hyperbolic | 43 groups (86.0%) | 45 groups (86.5%) |
| AICc prefers original piecewise quadratic | 7 groups (14.0%) | 7 groups (13.5%) |
| New model mathematically integrable fits | 96.0% | 97.1% |
| Old model mathematically integrable fits | 100% | 100% |
| Median postfit full-area ratio (new/old) among finite pairs | 1.0104 | 1.0112 |

AICc applies only to the assumed independent, identically distributed log-residual likelihood. Adjacent mass-profile samples can be correlated. Furthermore, the new model has `3m+1` free coefficients compared with `m+3` for the old shared-shape parameterization. Although AICc applies a nominal complexity correction, superiority on these descriptive metrics is not proof of a more accurate physical peak shape.

## Short profiles: an important failure mode

At only five observed samples per trace, the hyperbolic whole-line model often cannot be declared physically integrable (`C_j>|D_j|`):

| Five-sample pooling | Injection 1 integrable fraction | Injection 19 integrable fraction |
|---|---:|---:|
| One trace | 38.3% | 36.9% |
| Two traces | 39.7% | 45.3% |
| Four traces | 53.9% | 47.5% |

This is a meaningful limitation. Such fits must not be passed to a Bessel-based whole-line area computation. Full-area estimation at five points is model-sensitive regardless of speed.

## Kernel timing on the real candidate sets

From the two-injection Action run, Release on the hosted Ubuntu runner (all evaluated windows, repeated kernel calls):

| Timing | Injection 1 | Injection 19 |
|---|---:|---:|
| Existing qPeaks multitrace OLS | 13.71 ns/window | 13.78 ns/window |
| New hyperbolic OLS | 18.01 ns/window | 18.41 ns/window |
| New / existing | 1.31× | 1.34× |

These are microbenchmarks of short per-group window arrays, not full mzML throughput. Alternate CI runs produce different absolute times.

## Scientific interpretation

1. The smooth hyperbolic family improves the penalized in-sample log fit for **many sufficiently sampled pooled groups** under these assumptions. This is descriptive evidence of modelling flexibility, not validated accuracy.
2. The original qPeaks model is often preferable for **single-trace groups** and very short windows, and it yields many more integrable fits at five points.
3. The raw-amplitude step is an exact scalar intensity-domain OLS correction; it does not alter shape or apex. It never increases raw RSS on the retained fitting window, but it need not reduce full-area bias.
4. Similar model areas do not establish true peak areas. Both evaluate tails beyond observed samples and use the research dimensionless coordinate. Native m/z-domain integration requires a justified Jacobian.
5. A fixed `u` coordinate scale is required for completely linear regression. Its instrument interpretation is unresolved; do **not** optimize it against the same candidate peaks and then quote uncorrected OLS uncertainties.
6. Multi-trace pools consist of neighboring ions, not verified translated copies of one shared instrument response. Grouping by mass does not establish common shape.
7. There are no independent peak-mass, absolute recovery or reference-area annotations in this public dataset. Neither relative 1% model-area agreement nor AICc establishes a validated quantification improvement.

## Release gates (not yet satisfied)

- Independently justified instrument-coordinate scale and checks against true native sampling.
- Validation on certified mass/area standards across intensity, m/z, resolution and instrument classes.
- Repeatability with actual replicate standards rather than only serum injections.
- Study of correlated log errors, baseline/zero handling, overlap, finite-window area/tail instability and selection bias.
- End-to-end propagated uncertainty for the amplitude-corrected Bessel area.
- Thresholds for rejection of unidentifiable or nonintegrable hyperbolic fits.
- Only then decide whether any component belongs in the production qPeaks pipeline.

This is intentionally a research note; the original production code remains unchanged.
