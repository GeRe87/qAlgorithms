# Real-profile audit: zero handling and residual-model limits

## Scope and provenance

Base: `03649958394fecc759de802ec429d4ecd5f27d95`, clean tree.
Branch: `agent/peak-model-real-profile-audit`. Separate Draft PR against
`agent/peak-model-declared-search`. No production/model/detector/search/statistics
implementation changes, and no changes to PRs #1–#8 or their branches.

Inputs are the two checksum-verified public files from
[Zenodo 20729183](https://zenodo.org/records/20729183): TripleTOF 5600+ serum-pool
profile spectra, injection 1 and injection 19 from the same analytical run.
These are not certified standards, independent instruments or known-area samples.
Published MD5 values: 4130309a9283c8b492660dbad1fd3670 and
8bce27a1fdbe2ea0f0659d398b0bfaaf. Both full SHA256 values are in the committed
summary files and CI input manifest. Raw source files are not committed.

`experiments/peak_model_audit/PROTOCOL.md` fixes the actual protocol. There are
80 evenly scan-order-stratified positive profile MS1 scans per injection, up to
64 fixed index centers/scan and separate data-selected observed maxima >=100,
each with nine- and seventeen-sample windows. No interpolation, zero filling,
step tuning or physical-coordinate conversion. Dimensionless step is .35,
matching previous research. All attempted windows and rejection reasons are
recorded. Eligibility requires finite positive samples and sufficiently regular,
uninterrupted native sqrt(m/z) spacing. Eligible subsets are censored/selected;
they do not become prespecified null-testing families.

## Measured findings (local Release)

MSYS2 GCC 15.1.0, Intel i7-11850H, NumPy 2.5.4. Summaries in
`experiments/peak_model_audit/results/`; CI separately reproduces the audit and
uploads raw JSONL observations, residual CSVs, summaries and exact source HEAD.

| Input | Sampled scans | Profile samples | Exact zero intensities | Full-scan exploratory search |
|---|---:|---:|---:|---|
| Injection 1 | 80 | 142917 | 72279 | 80/80 rejected |
| Injection 19 | 80 | 146059 | 73945 | 80/80 rejected |

No negative or nonfinite intensities were observed in those sampled scans.
The new search requires complete positive uncensored traces, so its rejection
is expected contract enforcement. This is direct evidence that the current
interface cannot simply be connected to these complete raw scans. Replacing
zeros by epsilon, taking log(I+1), or silently deleting samples would introduce
a different observation/selection model and is deliberately not implemented.

| Cohort | Injection 1 eligible/attempted | Injection 19 eligible/attempted | Integrable I1 / I19 |
|---|---:|---:|---:|
| Fixed-index, 9 samples | 235/5120 | 223/5120 | 72/55 |
| Fixed-index, 17 samples | 42/5120 | 48/5120 | 16/26 |
| Observed maximum, 9 samples | 563/5049 | 680/5070 | 263/312 |
| Observed maximum, 17 samples | 93/5049 | 143/5070 | 56/85 |

For these files every rejected local window contained a nonpositive observation;
the geometry screen did not add rejections among the remaining positive windows.
All eligible raw-amplitude refits were numerically valid, but many fitted tails
were not integrable. A computable standard error is only a fixed-window model
reference, not validated coverage on these selected real profiles.

## Residual comparison: projection matters

Residual lag-1 is sum(e_i*e_(i+1))/sum(e_i^2), not a fitted AR coefficient.
For independent errors, OLS residuals still have covariance sigma^2 R with
R=I-QQ^T. The audit therefore compares residual shape to 10000 independently
seeded projected Gaussian draws R*epsilon per geometry. It does not compare
residuals to an incorrectly assumed zero-correlation baseline.

| Window | Projected-iid median lag-1 | Central 95% reference range |
|---|---:|---:|
| 9 samples (5 residual DF) | -.4460 | [-.8238, -.0013] |
| 17 samples (13 residual DF) | -.2348 | [-.6223, .1988] |

Observed fractions outside this reference range:

| Cohort | Injection 1 | Injection 19 |
|---|---:|---:|
| Fixed-index, 9 samples | 20/235 (8.51%) | 11/223 (4.93%) |
| Fixed-index, 17 samples | 29/42 (69.05%) | 34/48 (70.83%) |
| Observed maximum, 9 samples | 48/563 (8.53%) | 62/680 (9.12%) |
| Observed maximum, 17 samples | 72/93 (77.42%) | 103/143 (72.03%) |

These are DESCRIPTIVE reference-envelope comparisons, not adjusted p-values or
estimated raw-noise correlation rates. Windows overlap, scans are not asserted
independent replicates, positivity eligibility censors both cohorts, and observed
maxima are explicitly data-selected. Deviations can reflect shape misspecification,
interference/background, selection, heteroskedasticity, sampling or correlated
noise. This audit cannot identify one cause or certify a noise model. In
particular the scarce 17-sample positive subset is not representative of all peaks.
General residual-model assessment guidance: [NIST](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd441.htm).

Both kernel and reference fit exactly the emitted float log observations. An
independent Householder [QR fit](https://numpy.org/doc/stable/reference/generated/numpy.linalg.qr.html)
checks every eligible kernel coefficient (tolerance 2e-7 relative-plus-absolute).
Maximum normalized differences were 5.90e-14 and 3.52e-14. Thus large residual
patterns in these selected wider windows are not explained by a disagreement
between the running kernel and the independent least-squares solution.

Width sensitivity is also descriptive: among finite whole-line-area pairs for
observed-maxima windows, 17/9-sample area medians were 1.0666 (33 pairs) and
1.0835 (66 pairs); 10th–90th quantiles [.9076,1.1910] and [.9200,1.2898]. These
are different model extrapolations in u, not known-area accuracy or mass-domain
integration. No per-peak false-discovery, mass-accuracy or area-coverage target
is established by these files.

## Implementation and independent checks

Only an opt-in research executable, Python analysis/reference tests, fixed input
fetch/checksum script, dependency pin, protocol/results and dedicated workflow
are added. CMake option QALGORITHMS_ENABLE_REAL_PROFILE_AUDIT defaults OFF.
NumPy is an optional audit dependency only; default production builds are unchanged.
New C++ uses existing indentation/braces/naming and the native mzML reader and
existing typed APIs. No library mathematical equations are copied or rewritten.

Before edits: Debug/Release 16/19 passing with research ON. With audit enabled:
17/20, including the independent projection tests; exactly the existing
`test_qpeaks`, `test_qpeaks_pt03`, `test_qpeaks_pt04` failures remain. No expectations
are modified. Projection tests check orthogonality, idempotence, residual rank,
independent per-draw least-squares simulation and appropriate reference-envelope
frequency. A second simulation verifies that residual correlation differs from
raw iid noise, preventing the zero-correlation interpretation error.

Dedicated CI builds the immutable base/current sources in Debug and Release,
records full suites and exact legacy/27-file CLI differentials, runs nine paired
common benchmarks using normalized source paths/shared invocation path, and
records common executable .text hashes. A required subsequent real-data job
fetches/checksums both sources and uploads the complete evidence. It is not marked
continue-on-error; a data/network failure cannot be presented as completed evidence.
The exact final tested HEAD, timings, CI links and any local/CI numeric differences
are recorded in the Draft PR.

Reproduce with an isolated Python environment containing requirements.txt, CMake
option ON and Python3_EXECUTABLE pointing to that interpreter. Run real_peak_audit
INPUT.mzML OUTPUT.jsonl, then analyze.py OUTPUT.jsonl INPUT.mzML OUTPUT_DIRECTORY.
fetch.py downloads only the two declared files and validates published MD5 before
writing them. The result summaries retain NumPy version, seeds, input hashes,
all rejection counts and per-scan reference-envelope fractions.

## Decision

Keep production routing blocked. The next scientific design task is an explicit
observation model for zeros/detection limits/background and its selection rules,
followed by suitable real standards/replicate validation. This audit does not
justify silently changing the log model, whitening data, adding a baseline,
relaxing validity gates or interpreting wider-window residual structure as a
measured AR parameter. Nonlinear physical coordinates/Jacobians and calibrated
area uncertainty remain independent unresolved requirements. No merge/release.
