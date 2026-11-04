# mzML feature-list working branch

Branch: `work/mzml-featurelists` in `GeRe87/qAlgorithms`.
Starting commit: `ceec9eb6459253660a17caf5f80909f8f3ab143a`.

## Purpose and status

This branch establishes a reproducible command-line workflow from LC-HRMS
MS1 profile mzML to per-file feature CSVs. It is an experimental engineering
baseline. It is not a scientifically validated release.

The complete synthetic workflow detects the two known signals at their expected
m/z and retention time in positive, negative, and mixed-polarity files. A linear
null-model correction and conservative selection among existing valid fits now
also handle several synthetic overlapping-peak cases. Complete processing of
representative real profile data is still blocked by ambiguous EIC peak
selection. A diagnostic partial export is not a successful complete CLI run.

The legacy qPeaks suite retains unresolved cases (`pt_04`, expected three peaks,
and `pt_03`, expected no peaks). The reflected three-peak grouping regression is now fixed and included in
the regular numerical regressions. The full test suite must therefore
still be reported as failing. Test expectations have not been relaxed.

Unsupported unresolved apex groups still produce an explicit processing error. With
`--skip-error`, subsequent files are still processed and the batch returns a
nonzero status. These difficult profiles must be resolved before treating this
branch as a general-purpose working release.

### Peak-selection changes and their limits

The linear null-model residual now evaluates `intercept + slope * x`, matching
the slope returned by the linear regression. The previous sign inflated the
null-model residual for nonzero slopes. An exactly linear null model is valid,
and cannot justify an additional peak model with an equal or larger residual.

Only apex groups that the existing selector marks as conflicting enter the
new refinement. It associates already validated fits with observed local
maxima near their fitted apex, treats a flat maximum as one observation, and
compares alternatives on a shared raw-intensity SSE interval. Observed valleys
bound that comparison interval; they do not create a new fit. Fit coefficients,
original windows, degrees of freedom and uncertainties are retained. Previously
selected stable representatives are protected, and the same supported maximum
is emitted only once.

Nearby alternatives are combined only if their complete continuous fitted-apex
range is within the existing resolution limit and their windows overlap. A
final check uses the actual selected fitted apices and overlapping windows.
Missing observed support, missing separating valleys, and incompatible selected
apices remain explicit errors. Statistical acceptance thresholds are unchanged.
This is a bounded selection of peak descriptions, not quantitative deconvolution
or proof that two reported features correspond to two chemical compounds.

## Mandatory-group infeasibility regression

A local investigation of the first unresolved real EIC distinguishes failure of
the search from an infeasible grouping hypothesis. The search currently requires
one accepted description for each associated observed local maximum. For three
overlapping descriptions separated by more than the existing resolution limit,
the overall fitted-apex range must exceed twice that limit. No search strategy
or increased state budget can repair a hypothesis violating this condition.

A synthetic selector regression captures three mandatory groups with 4, 6 and 7
alternatives, spanning only 5.95 samples overall. Exhaustive enumeration of all
168 triples proves none feasible under the unchanged separation rule. The
selector must report explicit infeasibility instead of succeeding by dropping a
mandatory group. This originally established honest strict-mode error handling. The optional
collective output described below now handles this case without asserting
three resolved peaks; strict API calls retain the original error.

The corresponding measured trace, candidate coefficients, plots and exact
constraint proof are supplied as local analysis artifacts, not committed as
measurement data. Overlapping descriptions can be alternative accounts of one
signal region; the association with observed maxima does not itself
establish the number of separate physical peaks. The optional collective policy below preserves the established resolved-peak
tests. It does not infer a physical peak count or a new baseline rule from
this trace.

The new regression passes in native Debug and Release. Full CTest remains 6/10
with all four previously failing targets and original expectations retained.
That test-only checkpoint did not alter production output. The subsequent
collective implementation below extends the output contract explicitly.

## Explicit collective descriptions (native Windows continuation)

The user approved marked collective descriptions for statistically accepted
regions that cannot be represented as separate compatible peaks. Feature
construction opts into this policy; the default `qpeaks_find` API and all
centroiding remain strict. Statistical thresholds, accepted coefficients,
equidistant sample-index regression and RT/area transformation are unchanged.

Joint selection first seeks the original complete assignment. A failed search
can then return one **already statistically accepted** representative if it
covers the entire observed-maximum/fitted-apex envelope, keeps its apex in that
component's original apex range, and respects selected neighbours. Previously
accepted broad bridge alternatives may also represent this envelope. The
representative minimises raw-intensity SSE on that same envelope. This envelope
is not the union of unrelated alternative fit-window tails. No new fit is
created and no supported maximum in the replaced component is silently dropped:
all its original maximum anchors (first point for a plateau) are recorded.
Protected stable fits cannot be replaced. A missing full-envelope fit, an
unsupported maximum, or a protected-fit conflict remains an explicit error.

The main CSV appends five numeric columns to the original twenty:
`isUnresolved`, `unresolvedReason`, `unresolvedMaximaCount`,
`unresolvedLowerRT`, `unresolvedUpperRT`. Resolved rows contain zero in all five.
For a collective, reason `1` means the joint assignment was proven infeasible;
reason `2` means the bounded search stopped before proving feasibility or
infeasibility. These are not evidence for a particular number of chemical peaks.
The `*_unresolved_groups.csv` companion records `featureID,binID,maximumNumber,
retentionTime,reason`, one row per associated observed maximum, in output ID
order. It is created with a header even when no collective occurs. Successful
collective output emits an explicit CLI warning. **Area, height, position and
uncertainties in a collective row belong to its representative fit. They are
not summed analyte quantities, deconvolved values or validated total areas.**

Both output files are required. An unwritable companion returns a failure;
files already written before that failure are not a successful output pair.
`--skip-existing` preserves a complete pair and rejects an orphan main file or
orphan companion. It checks regular-file presence, not the scientific validity
or freshness of previously existing contents.

Two selection details are verified separately. After close maxima merge, the
valley comparison now uses the actual right boundary maximum of the left group,
not its possibly weaker first maximum. The previous condition rejected a
synthetic valid valley; reversing this one correction reproduces the failure.
Joint search also removes alternatives that have no compatible witness in
another remaining group. For the unchanged disjoint-window/apex-separation/
protected-fit predicate, five extremal witnesses are sufficient; tests compare
this reduction with exhaustive domain checks, including protected exceptions.
A shared-window/apex-range bound proves certain components infeasible without
exhausting the search budget. Neither optimisation changes the objective,
resolution threshold, accepted fits or budget.

Native GCC 15.1 Debug and Release both build. Full CTest remains **6/10**:
`test_qpeaks`, `test_qpeaks_pt03`, `test_qpeaks_pt04`, and `test_utils` retain
their original failures. The numerical regressions and all **15** generated
CLI scenarios pass in both configurations. New checks cover strict infeasibility,
marked collective output, an explicit search-limit reason, complete maximum
metadata and coordinate offsets, a validated bridge fit, rejection when no
full-cover candidate exists, protected neighbours, the merged-maximum valley,
output-pair failures and preservation of an existing pair. The earlier Linux
verification below is historical and is not the current native test result.

On the real `230222_SFC_Modelling_MixF_01_20.mzML` (SHA256
`b072c6bf1331bed1b1d6cb74ccaaad9bad9640a5e1d8b3f2988887a5d82a8da9`),
isolated positive EIC 119 now succeeds with 63 descriptions: 49 ordinary and
14 collective. Main and companion CSVs are byte-identical in Debug and Release.
The isolated helper uses local bin ID 0 for original EIC 119 and assigns DQSC=1;
its CSVs are diagnostic, not production whole-file exports. Its measured input
contains 708 points; maximum anchors and bounds use their actual EIC RT values.
The complete Debug and Release runs both import all 715 positive MS1 scans, produce 579351
centroids and 7017 EICs, then fail explicitly at EIC 180 because no accepted
fit covers its full unresolved envelope with an apex inside the component
hypothesis range. Neither reaches negative-mode
feature construction or produces a complete positive or negative
feature list. Earlier successful centroid-only runs covered 715 scans per
polarity; that does not establish feature-list completeness. No partial
successful-EIC export is promoted to a complete measurement result.

### Complete per-polarity audit of the collective policy

A private negative-only entry point changes only the polarity selection in the
CLI, and links the same verified core library. It completes all 9490 negative
EICs from 715 scans and 608532 centroids: 18640 descriptions, of which 1894 are
collective (1880 infeasible, 14 search-limited), with 21187 maximum-anchor rows.
All numeric fields, collective counts, bounds and companion references pass
the exported-CSV audit. Debug and Release main/companion CSVs are byte-identical across the entire
negative list. The private entry point reports 1307.683 s CPU time in Debug
and 93.249 s in Release; these are observed runs, not a controlled benchmark.
The state budget is not a wall-clock time guarantee.
This is complete processing under the current description policy, not a count
of resolved compounds or a quantitative validation.

A separate Release diagnostic evaluates every positive and negative EIC, using
the original bin ID in its output and retaining complete rejected-EIC manifests.
It accepts 7011/7017 positive EICs; IDs 180, 432, 1664, 2190, 3494 and 3902 fail.
Its positive `PARTIAL` list contains 19692 descriptions, including 2689
collectives (2672 infeasible, 17 search-limited) and 32548 maximum-anchor rows.
All six rejected traces are retained locally. Successful EICs without accepted
peaks still count as processed; the missing six EICs may contain an unknown
number of signals. The diagnostic returns status 1 for this partial result.
Its separate negative audit has zero rejected EICs and matches the negative-only
result by complete row contents and associated maxima. Equal-RT ties can receive
different IDs when EICs are sorted individually before the final non-stable
`std::sort`; this does not alter the descriptions or their metadata.

EIC180's first failing component covers indices 388..419, RT
456.527649..491.545532 seconds, with 12 associated maxima. Of 41814 accepted
candidates in that continuum, 9812 eligible fits cover the full region; none
has its fitted apex inside the component hypothesis. The closest covering
candidate has apex 428.653839. Thus coverage of a window alone is insufficient;
allowing an apex outside its component would require changing membership and
neighbour guarantees. The exception's short coverage wording includes these
additional eligibility conditions; the local proof and plot make them explicit.

A per-description audit retains sample start/length, smallest/largest actual EIC
RT interval and interval coefficient of variation (CV). All 19692 positive
partial rows and all 18640 negative rows have non-identical observed intervals.
Median CV is 0.05168335 positive and 0.04093320 negative; 6288 positive and 6225
negative windows have CV above 0.2. These are descriptive interval variations,
not estimated area errors or statistical acceptance thresholds. Missing EIC
scans and acquisition timing both contribute. No area correction was introduced,
and no row is labelled quantitatively area-validated.

RT interpolation and area corrections remain deferred. Irregular RT intervals
still prevent quantitative validation of areas and uncertainties. The next
selection issue is representation of a conflict region that lacks a single
accepted full-envelope model; the current implementation deliberately reports
that limitation instead of assembling unvalidated summed quantities.

## Statistical-cascade audit continuation

Retain the equidistant sample-index regression, its precomputed matrices, and
existing RT/area transformation. RT interpolation and area corrections are
explicitly deferred by the user.

Two concrete implementation errors are corrected without changing the cascade's
0.05 significance level:

- Log-domain residuals are captured from the original fitted coefficients,
  before the raw-intensity height correction changes the intercept. The latter
  optimises a different error measure and must not be treated as log-fit noise.
- `F_stat` no longer overwrites the requested alpha with 0.05. Its existing
  0.05 cache is used only for that level; other levels use their own calculation.

Exact zero log residuals are valid. A perfect complex fit is accepted against a
nonzero linear residual, while an equally exact line cannot justify additional
parameters. Zero-curvature/zero-residual candidates remain insignificant. No
artificial noise floor is introduced. The existing uncertainty calculation now
also permits a zero MSE; no area formula or RT correction is added.

Independent regression checks use a five-point residual orthogonal to all four
design columns: [1,-4,6,-4,1]*0.0111 on a quadratic log peak of curvature -0.5.
The exact RSS is 70*0.0111^2 and F is 202.9056..., above the unchanged 199.5
threshold. The previous height-corrected residual rejected this fixture.
Additional tests use the closed-form F(2,1) upper-tail quantile
(1/alpha^2-1)/2 for alpha 0.1, 0.05 and 0.01, including cache isolation.

Joint selection now forward-checks surviving alternatives, computes an
admissible compatible-error bound, and branches on the smallest remaining
domain. The existing unconstrained error bound is retained as a cheap initial
prune. Its objective, compatibility rules and explicit 200,000-state limit
are unchanged. A three-group test has exactly one feasible assignment among
27 combinations (apices 2.5,6,9.5). Removing its sole feasible left alternative
must raise an explicit infeasibility error rather than omit a group.

The original pt_03 and pt_04 expectations remain failing (2 instead of 0, and
2 instead of 3). The pt_04 candidate at center 3/scale 2 has original log RSS
0.0816954180, null RSS 16.8265746760, F=102.4835888111, and F-critical=199.5.
Its curvature test passes, but the unchanged F test rejects it. Lowering that
threshold simply to obtain the legacy count is not a verified correction.
The prior raw-height-corrected RSS was 0.1533269169. A candidate-level CSV audit
is supplied locally. Broad/background detections, statistical calibration over
many tested windows, and the covariance treatment of clipped fitting windows
remain unresolved; this is not complete validation of the cascade.

The current strict CLI still exits nonzero in the positive polarity at EIC 119
with no compatible accepted descriptions for all supported maxima. EIC IDs
cannot be compared directly with earlier revisions because the changed
acceptance already changes centroid and bin construction. Current complete
centroid-only processing succeeds: 579,351 positive and 608,532 negative
centroids over 715 scans per polarity. Debug/Release centroid CSVs are byte-equal.

The complete Release diagnostic retains 5,090 positive PARTIAL feature candidates
from 7,017 EICs (688 fail), and 8,099 negative candidates from 9,490 EICs
(638 fail). All output numeric fields are finite, areas/heights positive,
feature IDs unique, RTs sorted, bounds consistent with observations, and no
feature refers to a failed EIC. Failures are positive 541 infeasible/100
search-limit/47 missing-valley, negative 509/90/39 respectively. The diagnostic
iteration completed in 1,026.061 seconds including raw-EIC CSV output. This is
not a successful strict feature run or a speedup claim; the input to binning
and the accepted candidate sets differ from earlier revisions. The diagnostic's
inherited exit zero and success line refer only to its iteration, not completeness.
The full real-feature diagnostic was run in Release; Debug verification covers
the full test suite, centroiding and strict first-failure path, not all EICs.

The search-budget control uses 2,501 descriptions of an infeasible three-maximum
system. The old fixed order hits the unchanged 200,000-state limit; forward
checking proves infeasibility. A controlled single execution took 2.62 seconds
old versus 0.06 new. This does not establish a general real-data runtime gain.
The full native Debug and Release suites still pass 6/10 targets, including all
numerical regressions and the 11 synthetic CLI scenarios. The four original
failing targets remain failing with their original assertions.

## Joint fit selection continuation

Independently minimising fit error for each supported maximum can select
mutually incompatible apices, even when other already accepted fits satisfy
the existing resolution limit. The resolver now tries joint selection only
when those independent choices conflict. It minimises the sum of the existing
raw-intensity scoring-domain SSEs, retaining one unchanged accepted description
for every supported maximum/group and preserving stable representatives. Apex
ranges that cannot interact are solved independently. A search exceeding
200,000 visited states per interacting component remains an explicit processing
error; an infeasible component also remains an explicit error. Neither case
returns a reduced set of peaks as success. Acceptance thresholds are unchanged.

A reduced ten-point profile checks two retained maxima, unchanged separation
requirements, and reflection. A separate two-copy test checks four retained
maxima across independent components. The existing numerical and synthetic CLI
regressions pass in native Debug and Release. The full suite remains 6/10:
`pt_03`, `pt_04`, their containing legacy suite, and the exact beta-inverse
comparison remain failing with their original expectations.

The isolated original positive EIC 115 now completes with seven fits, and the
entire original negative spectrum 1917 completes with 900 centroids, in both
Debug and Release. These counts are reproducibility checks, not reference
chemical feature counts. The complete centroid-only CLI run now succeeds for
all 715 scans per polarity: 578,984 positive and 608,230 negative centroids.
The corresponding Debug and Release CSVs are byte-identical.

The strict full-file feature CLI still exits nonzero: its first unresolved
positive trace is now EIC 120 (no compatible accepted fits for all supported
maxima). The private positive diagnostic completes all 7,070 EICs, with 684
failed EICs and 5,006 retained PARTIAL feature candidates. This recovers 155
complete EICs relative to the previous snapshot. All 3,962 previously returned
positive feature descriptions are preserved byte-for-byte as CSV fields except
for sequential output IDs, and no previously successful EIC now fails.
All positive numeric fields
are finite and areas/heights positive, IDs unique, retention times sorted, and
no returned feature refers to a failed EIC. The 684 failures comprise 500
infeasible selections, 138 explicit search-limit errors, 45 missing separating
valleys, and one missing supported maximum. No failed EIC contributes partial
features. The negative diagnostic also completes all 9,516 EICs: 613 fail, and 8,138
PARTIAL feature candidates are retained. Its failures are 451 infeasible joint
selections, 130 explicit search-limit errors and 32 missing separating valleys.
All negative numeric fields are finite, areas/heights positive, IDs unique and
RTs sorted, and no feature refers to a failed EIC. In 5,275 negative feature
windows max(dt)/min(dt) exceeds 1.1. The diagnostic runtime was 414.105 seconds
for both polarities on this Windows machine; this includes writing raw EIC
observations and is not a strict successful processing benchmark. Its inherited
final success message and exit zero refer only to completing the diagnostic
iteration; the PARTIAL labels, failure manifests and strict CLI exit one remain
authoritative about completeness. The earlier counts below are historical snapshots.

In 3,468 positive feature windows max(dt)/min(dt) exceeds 1.1, using measured
points within CSV bounds with a 1e-4 second rounding tolerance. Quantitative
area/shape/uncertainty validity remains unresolved for all features. The joint
selector changes fit selection, not the even-spacing assumption or fit
coefficients. It is not a quantitative deconvolution method.

## Native Windows continuation, 2026-11-04

Source baseline: `c00e415c4a0520f6a87cc1af27aaacad11a27818`.
MSYS2 UCRT64 GCC 15.1.0, CMake 4.0.3, Ninja; default architecture flags.

Apex grouping now computes connected components using pairwise apex proximity
and overlapping fit windows. The old implementation kept the seed window fixed
while expanding the apex range, so the partition depended on candidate order.
A chain of overlapping windows is now checked in all 24 candidate permutations
and both reflections, plus singleton assignment. The original full-pipeline
reflection regression is included in the normal regression target and passes.
Fit acceptance thresholds and explicit unresolved-conflict errors are unchanged.

The vendored beta routines now use C++ `std::lgamma` for positive shape
parameters; their unused sign output required unavailable `lgamma_r` on UCRT64.
Test assertion output is flushed before the existing breakpoint trap, so captured
Windows test logs retain the actual failing assertion.

| Verification | Debug | Release |
| --- | --- | --- |
| Native configure/build | Pass | Pass |
| Full CTest suite, including separate legacy diagnostics | 6/10 pass | 6/10 pass |
| Numerical regressions including reflection/permutations | Pass | Pass |
| Synthetic mzML CLI scenarios | 11/11 pass | 11/11 pass |
| `pt_03`, expected zero peaks | FAIL: 2 | FAIL: 2 |
| `pt_04`, expected three peaks | FAIL: 2 | FAIL: 2 |
| `test_utils` exact beta-inverse equality | FAIL | FAIL |

The four failing targets are `test_qpeaks`, `test_qpeaks_pt03`,
`test_qpeaks_pt04`, and `test_utils`. Separate legacy targets ensure the first
failure cannot hide the other case. No tests are disabled or expected to fail.
On this Windows math library `incbi(0.5, 1, 0.5)` is
`0.24999999999999989`, rather than bit-exact `0.25`. The original exact assertion
is retained. A separate diagnostic checked 24 closed-form beta/inverse-beta
values with absolute error below 1e-12; this does not make `test_utils` pass.

Candidate-level inspection shows that `pt_04` has only two accepted fits before
grouping. The fit around sample 3.27324 (zero-based) is rejected by the existing
F-test; grouping cannot recover it. `pt_03` admits broad near-baseline fits and
selects two features. Resolving these acceptance problems requires further
statistical/model work rather than tuning thresholds to the two fixtures.

### Real sample now available

The user supplied the Windows location of
`230222_SFC_Modelling_MixF_01_20.mzML` (270,841,344 bytes), SHA-256
`b072c6bf1331bed1b1d6cb74ccaaad9bad9640a5e1d8b3f2988887a5d82a8da9`.
Independent XML inspection confirms 715 MS1 scans per polarity. Scan intervals
range from 0.618403 to 1.407683 seconds positive and 0.648982 to 1.485695
seconds negative; 180/181 of 714 intervals differ by over 10% from their
polarity median. Areas and uncertainties remain quantitatively unvalidated.

Strict Debug and Release stop at positive EIC 115. A private diagnostic, which
retains features only when an entire EIC succeeds, inspected all 7,070 positive
EICs: 839 failed and 3,962 features were retained. All numeric CSV fields are
finite and areas/heights positive. The CSV is explicitly PARTIAL. Failures are
793 unresolved selected-apex distances, 45 missing separating valleys, and one
missing supported maximum. Original EIC IDs are preserved in a separate failure
manifest; the inner singleton diagnostic error calls its local EIC zero.

The diagnostic cannot produce a negative feature list: negative centroiding
fails at original spectrum index 1917 with the unresolved selected-apex-distance
error. No spectrum is skipped. Historical negative feature counts therefore
are not reproduced by this revision and must not be reused as current results.

In 2,717 of the 3,962 positive feature RT windows, the largest observed interval
exceeds the smallest by more than 10%. This audit uses measured points between
the rounded CSV RT bounds (1e-4 second inclusion tolerance). It is a spacing
indicator, not an area-error estimate. No feature area is analytically validated.

The group traversal now indexes unassigned candidates by apex position, removes
each assigned candidate from that index, and checks only nearby candidates.
This preserves the pairwise graph definition while avoiding repeated full-list
scans on long real EICs. An independent all-pairs transitive-reachability oracle
checks 50 sets of 24 windows, in addition to the original permutation/reflection
tests. Native Debug/Release retain the same six passing and four failing targets.

For writable temporary test files in PowerShell, set `TEMP` and `TMP` to an
existing writable directory before CTest. The first attempt using the restricted
system temporary directory failed with access denied; reruns in the workspace
passed all 11 CLI scenarios. The historical Linux results below describe the
baseline and are not a new Linux verification of these changes.

## Build

Use a C++17 GCC toolchain, CMake, and the zlib development library. Python 3 is
needed for the generated mzML integration test. Vendored C++ dependencies are
included in the repository. The initial Linux build uses GCC 13.3 and zlib 1.3.

```sh
git clone --branch work/mzml-featurelists https://github.com/GeRe87/qAlgorithms.git
cd qAlgorithms
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

The final command includes the unresolved legacy test described above. To run
the new numerical regressions and complete mzML smoke test separately:

```sh
ctest --test-dir build --output-on-failure -R "test_qpeaks_regressions|test_mzml_cli"
```

Debug and Release configurations are respected. CPU-specific compilation is
opt-in (`-DQALGORITHMS_NATIVE_ARCH=ON`); the default build does not use
`-march=native`. Standard zlib is the default. Native zlib-ng can be selected
explicitly with `-DQALGORITHMS_USE_ZLIB_NG=ON` when its CMake package is installed.

### Windows build recipe

Use the **MSYS2 UCRT64** environment consistently. The following packages belong
to the same environment:

```sh
pacman -S --needed git mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-zlib mingw-w64-ucrt-x86_64-python
```

Then run the clone and CMake commands above in that UCRT64 shell, adding
`-G Ninja` to the CMake configuration command. There is no hard-coded MSYS2
installation path in the build anymore. This recipe still needs a native
Windows validation run; it is not evidence of a completed Windows test.

The resulting executable is `build/qAlgorithms.exe`. Run it in the UCRT64 shell,
which supplies the runtime DLL path. For a PowerShell session with the default
MSYS2 installation, set the path for that session first:

```powershell
$env:Path = "C:\msys64\ucrt64\bin;" + $env:Path
.\build\qAlgorithms.exe -i "C:\Data\sample.mzML" -o "C:\Data\features" --printfeatures
```

## Process data

```sh
./build/qAlgorithms -i "/data/sample.mzML" -o "/data/features" --printfeatures
./build/qAlgorithms -i "/data/samples" -o "/data/features" --printfeatures
```

Input directories are searched recursively. The mzML extension is matched
without case sensitivity. Repeated paths are deduplicated; distinct named
samples are retained even when their byte content is identical. Output
directories are created as needed. Supplying `-i` and `-o` without an output
flag defaults to feature CSV. The older `-printfeatures` spelling and `-pf`
remain accepted.

`--printall` additionally writes centroid and bin CSVs. `--skip-existing`
preserves existing outputs. Use a fresh output directory for a new analysis;
without that flag, existing matching output files can be replaced. Distinct
inputs that would produce the same output filename are rejected. A mixed
polarity file writes separate positive and negative results.

| Input | Feature output |
| --- | --- |
| Single-polarity `sample.mzML` | `sample_features.csv` |
| Mixed-polarity `sample.mzML` | `sample_positive_features.csv`, `sample_negative_features.csv` |

The feature CSV contains ID, bin ID, m/z, m/z uncertainty, retention time,
retention-time uncertainty, RT bounds, area, area uncertainty, height, height
uncertainty, scale, quality fields, and regression coefficients. Retention
times are in seconds. **DQSB is still the upstream placeholder `-1`; do not
interpret it as a computed bin quality score.**

A valid profile input without detected features writes a header-only CSV and
returns success. Malformed input, unsupported input, unresolved apex groups,
or failure to write the feature CSV returns a nonzero status. Under
`--skip-error`, any failed file still makes the final batch status nonzero.
A later failure in a mixed-polarity file may leave an already completed CSV
from its first polarity; a nonzero run must not be treated as complete.

## Supported input and limits

- MS1 **profile** spectra, positive, negative, or interleaved polarity.
- Plain `mzML` and `indexedmzML` with ordinary mzML XML element names.
- 32-bit and 64-bit floating-point binary arrays, selected separately for m/z
  and intensity, with zlib or no compression. Array order is determined by CV
  accession; shared referenceable parameter groups are supported.
- Scan start times in seconds or minutes; selected scans must have strictly
  increasing retention times.
- Peak regression still assumes locally even spacing in the supplied samples.
  The current EIC construction does not resample irregular retention times or
  insert missing scans. Consequently, irregular acquisition timing is an
  unresolved limitation for peak shape, area and uncertainty estimates. Correct
  CSV coordinates alone do not validate those quantitative properties.
- Centroid-only input and MS2-only input are rejected. MS2 processing, MS-Numpress,
  cross-sample feature alignment, compound identification and deconvolution
  are outside this working baseline.
- Current internal index widths limit one polarity to 65,535 selected MS1 scans
  and a contiguous above-threshold profile region to 65,535 points. Larger
  inputs are rejected explicitly instead of wrapping indices.

## Reproducible verification

Verified on 2026-11-04 on Linux x86-64 with GCC 13.3, CMake 4.4.4 and
system zlib 1.3:

| Check | Debug | Release |
| --- | --- | --- |
| Configure and build | Pass | Pass |
| Full CTest suite | 7 of 8 targets pass | 7 of 8 targets pass |
| Generated mzML/CLI integration scenarios | 11 of 11 pass | 11 of 11 pass |
| qPeaks coordinate, indexing, weighting, linear-null and selection regressions | Pass | Pass |

The integration scenarios are grouped in the single `test_mzml_cli` CTest
target. The sole failing target is the retained legacy `test_qpeaks` suite,
which stops at `pt_04`: three peaks are expected, but two are reported. It has
not been disabled or marked as an expected success. The full CTest command
therefore returns a failure status in both configurations.

The added tests include exact linear trends and orthogonal residuals, analytic
single and paired Gaussian signals, broad peaks, close unequal peaks and their
reflections, monotone traces, plateaus, and preservation of a stable neighbour.
Expected peak counts and positions were fixed before evaluating the changes.
Passing count/position checks does not establish quantitative accuracy: one
weaker peak in a close synthetic pair still has about 13.1% area error.

The formerly failing reflection example is also directly runnable:

```sh
./build/tests_bin/test_qpeaks --known-reflection-failure
```

This full-pipeline reflected three-peak regression now passes its original
count and half-sample position requirements in native Windows Debug and Release.
It is also executed by `test_qpeaks_regressions`; no expectation was relaxed.
This establishes the specified regression, not universal reflection invariance.

The Python integration test generates all its own inputs using a fixed random
seed and standard-library XML, binary encoding and compression. It requires
known peaks near m/z 100.125 / 150.375 and RT 22.25 / 31.4 seconds for positive
mode, with distinct negative-mode coordinates. All numeric CSV values must be
finite and detected peak areas/heights positive. It also covers valid profiles
without detected peaks, unsupported spectra, malformed data, missing option
values, output-name collisions, batch error handling and default CSV output.

```sh
python3 test/integration/test_mzml_cli.py ./build/qAlgorithms --work-dir ./smoke-evidence
```

With `--work-dir`, generated mzML inputs, feature CSVs and process logs are
retained. Without it, successful temporary runs are cleaned up. This is an
engineering test of the processing chain, not a validation of analytical
accuracy, detection limits, false-positive rates or uncertainty estimates.

### Real measurements without a reference feature list

Such files are useful for checking complete import, processing, finite outputs,
repeatability, and agreement of reported peak locations with the original
profile signals. Record the input hash, exact source revision, command, exit
status and logs alongside each result. Select any manually inspected peaks by
a stated rule, rather than choosing only visually convincing examples.

Without known reference signals, these checks do not establish how many true
features were missed, how many detections are false, compound identities, or
quantitative accuracy. A raw extracted ion trace used for visual comparison
also has different intensity/area semantics from the fitted qAlgorithms
feature area; its extraction window must be reported explicitly.

The strict CLI continues to return an error for an unresolved EIC. A private
diagnostic may isolate each EIC, retain features only after the whole EIC
succeeds, and record every rejected EIC with its original ID, polarity, RT
bounds and error. Such files must be labelled `PARTIAL`. Successful-EIC
fractions are not feature recall, because excluded EICs may contain an unknown
number of peaks. The partial diagnostic is not a catch-and-continue production
mode and does not change the CLI's failure status.

## Next acceptance work

1. Resolve ambiguous real-data EIC selection, the legacy `pt_03` / `pt_04`
   cases without weakening expected
   behavior or silently discarding ambiguous regions.
2. Account for irregular retention-time spacing before validating areas,
   shapes and uncertainties quantitatively.
3. Achieve a complete strict CLI run on representative profile mzML data;
   then compare feature locations and quantities with known signals or an
   agreed independent reference analysis.
4. Resolve the remaining native Windows test failures, and measure
   time and memory on successful realistic-file runs before broader use.
