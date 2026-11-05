# Real mzML extension of the component comparison

The real profile mzML is usable for model diagnosis without a known physical
peak count. It supplements the controlled signals in [COMPARISON.md](COMPARISON.md).
It does not provide a false-positive/recall calibration by itself.

`export_eics.cpp` uses the unchanged production reader, centroiding and qBinning
path. It exports **every** EIC observation, with original RT, centroid area,
mass, measured-point count and scan number. There is no qPeaks filtering in
this exporter. Missing polarities and blank scans have explicit count rows;
unsupported inputs and invalid/truncated observations fail explicitly.
The binary is little endian: `QEIC0001`, followed by records of three uint32
values (polarity 1/2, zero-based EIC ID, length), then length entries of float32
RT, float32 area, float32 mass, uint16 measured-point count and uint32 scan.
There is no padding between vectors. `inventory.csv` and `polarities.csv`
account for all records and observations.

The local real run selected 715 MS1 profile scans per polarity and exported
7017 positive / 9490 negative EICs (368267 / 370781 observations). Native
Debug and Release exports, including inventories, are byte-identical.
All 38524 existing production feature windows were bound to these freshly
extracted EICs by polarity, EIC ID, point range and the existing nine-significant-
digit RT serialization. Duplicate feature/audit identities are errors.
This is an input/window consistency check, not validation of old peak labels.

Before fitting, `real_comparison.py` selects the median-length and longest
eligible prior window per polarity and old resolved/collective class, plus
the first eligible collective feature in each historical positive conflict
EIC explicitly supplied by CLI. Equal windows are deduplicated. Eligibility
`n>8` permits an unsaturated eight-parameter reference, not a peak threshold.
The local run has 14 windows, including six historical conflict windows.
The complete feature-selection ledger accounts for every other prior window.
Sampling is deliberately diagnostic and is not representative prevalence.

The same fixed-support, fixed-midpoint one/two-component model, all-split
initialization and budgets as the controlled comparison are used. Each
original measurement is counted once. No resampling, RT-domain fit, area
correction or automatic peak-count decision is introduced. Real intensities
here are centroid areas within an EIC, not integrated chromatographic areas.

All 14 windows have result and curve records. Best-found two-component fits
converged in 13; one longest positive collective window remains unconverged.
There are 2422 candidate/start records and 216 unsuccessful starts, retained
explicitly. In the unconverged case RSS decreases from 323.8033 to 143.4038,
but adjacent-observation log-residual correlation remains about 0.977.
Another longest positive window converges, yet retains correlation about
0.853. These are descriptive indications of remaining structure; they do not
separate a deficient mean model, optimization failure and correlated noise.
A larger RSS gain alone cannot establish additional physical peaks.

`audit_real.py` independently recomputes every finite start score, the convex
single-component solution by all active sets/SVD, fitted curves and component
parameters, and verifies all original selected vectors and window identities.
Maximum RSS discrepancy is 7.28e-12; the independent convex reference differs
by at most 3.56e-15. Seven accounting/export tests pass with each native build;
the 14 component-model tests also pass. Existing production expectations,
including the unresolved pt_03/pt_04 failures, remain unchanged.

## Reproduction

With the same libraries/optional Python dependencies as the preceding study:

```text
g++ -std=c++20 -O3 -DNDEBUG -Iinclude -Iexternal/libcerf experiments/qpeaks_ensemble/export_eics.cpp build-release/libqalgo_lib.a -lz -o export-eics.exe
export-eics.exe INPUT.mzML EXPORT
python experiments/qpeaks_ensemble/test_real.py --exporter export-eics.exe -v
python experiments/qpeaks_ensemble/real_comparison.py EXPORT FEATURES INPUT.mzML COMPARISON --historical-positive-bin ID
python experiments/qpeaks_ensemble/audit_real.py EXPORT DEBUG_EXPORT FEATURES COMPARISON AUDIT
```

Repeat `--historical-positive-bin` for the historical EIC IDs recorded in the
local protocol. Without that option only the class/length strata are selected.
FEATURES contains complete positive/negative `*_features.csv` and
`*_area_timing_audit.csv` files. The local deliverable contains those original
lists, hashes, full extracted EICs, the exact selected vectors, all fits,
all failed starts, diagnostics and the verified experimental source.

The next method question is whether background and/or a more flexible shared
profile explains the long-window residuals, together with optimizer and
support sensitivity. Native candidate grouping, calibrated component-count
selection, pt_03/pt_04, and chromatographic area validation are not solved by
this conditional real-data comparison. The complete production qPeaks feature
pipeline was not rerun here because its source did not change; the native
reader/centroiding/binning path was rerun on the whole real mzML in both builds.
