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

Unresolved apex groups now produce an explicit processing error. With
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
