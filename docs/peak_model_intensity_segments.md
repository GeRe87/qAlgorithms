# Intensity segmentation for hyperbolic searches

The branch starts at 6cfd2847b9067f975b222997f4e13493d1b64de3, the head of Draft PR #9.
PRs #4 through #9 were open drafts in their existing stack when work started.

## Legacy behavior

qpeaks_find in src/qalgorithms_qpeaks.cpp returns -1 for null intensity,
coordinate or result pointers and -2 for an input shorter than MINLENGTH.
It allocates length+1 log entries and logs the entire input before segmentation.
The selection predicate is the strict raw comparison I > 2.7182818284590452,
with a comment tying the choice to log intensities greater than one.
The negated predicate skips a prefix of invalid samples and the positive
predicate finds the half-open run [rangeStart, pointIdx). Equality, zeros,
negative values and NaNs terminate or exclude runs. Positive infinity satisfies
the legacy comparison. There is no separate finite-value gate in that loop.
It resumes after the first separator. Runs reaching length include the final
sample. MINLENGTH is five. A selected run longer than 65,535 raises an error
because the legacy coordinate representation has a 16-bit limit.

Each run passes offset raw, log and coordinate pointers into regression.
Optional cumulative degrees of freedom are rebased by subtracting the value
immediately before the run. Appended regression span starts and unresolved group
coordinates are translated back by rangeStart. The existing five-point boundary
test in test/test_qpeaks.cpp exercises a run occupying the entire input.
Other fixtures cover peaks sharing a continuum above e. The legacy coordinates
are assumed equidistant by the kernel. The segmentation loop does not interpolate.

## Segmented search

Call HyperbolicSearchPlan::evaluateSegmented with full AlignedIntensityTrace
views, the expected coordinate step, its relative tolerance and optionally a
raw intensity threshold. Each view supplies the original intensity and
coordinate arrays and length. The caller owns those arrays. Output window centers,
scales and half-open span bounds directly index the original arrays.
A finite strictly increasing coordinate grid and exact agreement of coordinates
across traces establish alignment. An edge outside the declared spacing tolerance
breaks a run. The dimensionless fixedStep of the fit remains a separate input.

The model-independent segmentIntensities function scans the aligned arrays once.
A sample is eligible when every trace is finite and above the threshold.
Runs are the intersection of valid sample indices across traces. Empty
intersections leave every affected window explicitly unevaluated. Per-sample
diagnostic priority across traces is nonfinite, negative, then below threshold.
The result records short runs and coordinate breaks. The search reports one
status per original window and retains every window and the original family size.
A linear merge of sorted window centers and spans checks whether each complete
window belongs to an eligible span. Stage C then receives its original raw slice
and logs the selected positive values. One geometry cache serves all windows
in an evaluation. The existing linear fit, raw amplitude correction, analytic
area and covariance calculations supply the inference results.

The threshold is in the supplied raw units. Positive gain changes can move the
same signal across a fixed threshold. An explicit threshold can reflect units
or known acquisition settings. The default e is a legacy selection convention.
A stored zero is an excluded observation with unknown physical interpretation.
Original data and index spacing supply the fit directly.

Selection by observed intensities is exploratory even for a prespecified plan.
The effective selection mode is reported in search.options. Nominal statistics
and covariance retain their fixed-window interpretation. Bonferroni uses the
whole declared family. Post-selection coverage, the physical meaning of zeros,
actual instrument error and tail adequacy require separate evidence.

## Numerical and real-profile evidence

The new test covers complete positive input, all zeros, threshold equality and
the adjacent representable values around e, short and terminal runs, interior
zeros, malformed intensities, a coordinate break and misaligned traces. It
compares shared and local running coefficients on a retained common run with
the direct positive window. Stage C comparisons use the direct raw slice and
require exact RSS, adjusted statistics, apex, corrected amplitude and height,
analytic area, area standard error and full covariance.

The real-profile tool reproduces the 80 scans per injection from the earlier
audit, using the same checksum-verified Zenodo 20729183 inputs. It searches
all centers at half-widths 2 through 8, fixed u step .35, and uses the lower-third
sqrt(m/z) gap with .45 relative tolerance as the coordinate spacing policy.
The wider-window fit uses the original index basis and the original coordinates
are retained by index. Physical coordinate calibration remains unverified.

| Metric | Injection 1 | Injection 19 |
|---|---:|---:|
| Complete scans evaluated | 80 | 80 |
| Zeros | 72,279 | 73,945 |
| Positive samples up to and including e | 0 | 0 |
| Above-threshold runs after coordinate splits | 39,578 | 40,574 |
| Runs with at least five points | 2,503 | 2,499 |
| Short runs | 37,075 | 38,075 |
| Declared windows | 994,819 | 1,016,813 |
| Evaluated windows with finite identifiable coefficients | 45,649 | 47,130 |
| Mathematically integrable windows | 10,179 | 10,166 |

The previous full-array search rejected all 160 scans. Matching the earlier
audit by original scan, center and scale retains all 933 and 1,094 previously
eligible cohort records, with exact agreement of integrability and Stage C
availability. A window can be recorded in both audit cohorts. The larger full
search uses more centers and widths, so these counts represent fits rather than
distinct peaks.

### Interpreting area comparisons

For these TripleTOF profiles the true peak areas are unknown. An independent
numerical quadrature of the same fitted function over the same domain can check
the analytic integral. Agreement verifies that calculation for the fitted
shape in dimensionless u.

Quadrature of the observed intensities over their actual coordinates provides
a descriptive area over the observed window. Its value depends on the window,
background and handling of sampling gaps. A comparison with the model needs
the same interval, coordinate units and background convention. The whole-line
fit additionally includes extrapolated tails. Comparing the fitted integral
over the observed interval with the whole-line integral separates that tail
contribution within the model.

Simulated profiles with known integrals can quantify recovery under specified
sampling, noise and shape assumptions. For the real profiles, report agreement
over the observed interval and sensitivity to the window and tails. Absolute
area error on these spectra is unidentifiable from the available observations.

## Performance and reproducibility

Local Windows Release with GCC 15.1.0 and Intel i7-11850H. Times are medians
over the 80 different scans per injection. Full evaluation includes segmentation,
geometry preparation, window inference and collection, and excludes plan
construction, file decoding and evidence writing.

| Component | Injection 1, us/scan | Injection 19, us/scan |
|---|---:|---:|
| Segmentation | 12.95 | 12.40 |
| Preparing logs for the running benchmark | 3.20 | 3.30 |
| Running regression and coefficient classification | 6.05 | 5.60 |
| Full segmented search | 1194.00 | 1240.75 |

The workflow separately builds the immutable base and current branch in Debug
and Release. It records full CTest suites, exact legacy probe and 27-file CLI
comparisons, paired legacy benchmarks and executable code hashes. The required
real-data job records checksums, prior-audit window records, per-scan timings,
matched windows and summaries. Build the standalone collector with
QALGORITHMS_ENABLE_REAL_INTENSITY_SEGMENTS=ON and use its three file arguments
for input mzML, scan CSV and window CSV. summarize.py takes the two CSVs and
the matching earlier audit JSONL. The final CI and exact tested commit are
linked in the Draft PR.
