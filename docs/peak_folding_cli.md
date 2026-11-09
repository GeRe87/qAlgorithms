# Experimental mzML peak-folding CLI

This **research-only** route uses the real mzML reader and the linear
folding solver already present in qPeaks. It is opt-in. The default
qCentroids -> qBinning -> qPeaks feature pipeline remains unchanged.

## Command-line usage

TOF/QTOF profile data:

~~~powershell
.\qAlgorithms.exe -i "C:\data\20171016_POOL_POS_1_105-134.mzML" -o "C:\data\folding" --peak-folding --folding-axis tof --folding-model asymmetric --folding-max-scans 10
~~~

Orbitrap profile data:

~~~powershell
.\qAlgorithms.exe -i "C:\data\my_orbitrap_profile.mzML" -o "C:\data\folding" --peak-folding --folding-axis orbitrap --folding-model symmetric --folding-max-scans 10
~~~

A single input or directory of mzML files is accepted using the existing
-i/-o interface. The output is input-stem_folding.csv, or separated into
input-stem_positive_folding.csv and input-stem_negative_folding.csv for
mixed-polarity files. **Without** --peak-folding, the default output
continues to be the standard feature CSV.

| Option | Meaning |
|---|---|
| --peak-folding | Select research folding instead of normal feature processing |
| --folding-axis tof | Native coordinate sqrt(m/z), required for TOF |
| --folding-axis orbitrap | Native coordinate -1/sqrt(m/z), required for Orbitrap |
| --folding-model asymmetric | Two shared half-curvatures (default) |
| --folding-model symmetric | One shared curvature |
| --folding-min-intensity 100 | Minimum absolute intensity of local candidate maximum |
| --folding-max-scans 10 | Process N MS1 spectra per polarity (default: all) |
| --skip-existing | Preserve previously written folding CSV |
| --skip-error | Continue on file-level error; overall exit remains nonzero |

The instrument sampling axis is **explicit**, never inferred from a name.
Folding cannot be combined with the usual centroid/bin/feature print flags.

## Processing semantics

For each selected MS1 profile spectrum and polarity:

1. Read m/z and intensity arrays with the existing mzML reader, without
   centroiding, interpolation, scan summation or artificial sample filling.
2. Transform sampled masses into monotone instrument coordinates. Estimate
   the local sampling-step scale using a lower quantile of neighbour gaps.
3. Select peaks with an observed maximum above threshold, no nearby
   competing maximum above the shoulder criterion, and at least two
   successive, positive-intensity native samples on each side. Use up to
   four original samples per side. Missing bins are never imputed as zero.
4. Perform one **linear pooled OLS fit per spectrum** for two or more
   selected peak profiles, with individual b0 and b1 but shared b2,b3
   (or b2=b3 in symmetric mode). The Schur system is only 1x1 or 2x2.
5. Emit one CSV row per selected profile, including its raw maximum,
   fitted peak coordinate when valid, its conditional model-based
   uncertainty, coefficients and pooled fit statistics.

If a scan contains too few eligible profiles, it is skipped. An input
with no valid grouped scan produces a header-only CSV and a warning.

## CSV contents

Includes spectrumIndex, retentionTimeSeconds, polarity, axis, model,
samplingStep, profilesInFit, anchorIndex, anchorMz, anchorIntensity,
samplesUsed, validApex, apexMz, apexUncertaintyMz, b0, b1, b2, b3,
b0StandardError, b1StandardError, b2StandardError, b3StandardError,
logRSS, pooledLogMSE and degreesOfFreedom.

Apex coordinates and uncertainties are blank for profiles that do not
meet the unimodal-concavity and in-window apex criteria. These results
are **not** a validated feature list or an externally calibrated mass
accuracy estimate.

## Research and numerical limits

- Sampled masses must be positive and strictly increasing; zeros and
  nonfinite intensities cannot enter the logarithmic model.
- Lower-quantile step estimation and neighbour/shoulder selection are
  heuristic instrument-profile screening, not validated centroiding.
- TOF uses sqrt(m/z); Orbitrap uses -1/sqrt(m/z). Check these transforms
  against actual native sampling before physical interpretation.
- The asymmetric curvature breakpoint is fixed at the observed local
  maximum; a profile-specific slope does **not** move it. Symmetric
  quadratic peaks are exactly translatable using local slopes.
- Classical OLS covariance assumes fixed coordinates and candidate
  selection with uncorrelated, homoscedastic log-domain errors. It does
  not account for calibration, sampling correlations, baseline,
  overlapping ions or selection bias.
- No peak linking across successive MS1 scans, chromatographic fitting,
  mass assignments or automatic IMS spectrum aggregation is attempted.
  mzML inputs with non-increasing selected scan retention times are
  rejected by the existing reader.
- All new outputs are explicitly research-grade. The qPeaks feature
  pipeline remains unchanged and should be run separately for comparison.

## Automated tests

The new CTest entry test_folding_cli generates deterministic synthetic
mzML profiles in test/integration/test_folding_cli.py, then verifies the
real executable's output, CLI handling, polarity splitting,
symmetric/asymmetric options, no-overwrite behavior, blank profiles and
centroid-only rejection.

Numerical regression and covariance correctness are covered separately
by test_qpeaks_folding and test_qpeaks_multitrace. Passing these tests
does not establish mass accuracy on independent real samples.
