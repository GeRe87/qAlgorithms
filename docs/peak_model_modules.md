# Native peak-model modules

The host loads a separately built shared library from an explicit path.
The library exports the C symbol `qalgorithms_peak_model_v1` and returns
the descriptor defined in `include/qalgorithms_peak_model_abi.h`.
The host checks the ABI version and descriptor size, then dispatches through
its callbacks. The descriptor supplies the model name.

## Source layout

- `qalgorithms_peak_module.cpp` owns the library handle and adapts the native fit to the qPeaks validation code.
- `detail/qalgorithms_peak_cascade.h` orders the shared validation stages.
- `qalgorithms_peak_background.cpp` compares raw backgrounds and replays the same model fit on rotated observations.
- `qalgorithms_peak_builtin_models.cpp` connects the enum-based models to the validation code.

The library owns its opaque fit and supplies its destructor. The host
destroys each fit before starting another one in the same validation.
Callbacks borrow input and output buffers for the duration of the call.
Concurrent calls require support from the module. Use modules you trust
to execute as native code in the qAlgorithms process.

## Run a fixed-window diagnostic

```text
qAlgorithms --peak-model-module PATH STEP MIN_AREA MIN_COVERAGE I0 I1 I2 I3 I4 [...]
```

Supply an odd-length prespecified window with at least five intensities.
The window center is coordinate zero. The module receives the declared
sample step and reports area in intensity times that coordinate's unit.
The command prints the model name and each stage's status. The C++ entry
point is `NativePeakModel::validateCandidate` in `include/qalgorithms_peak_module.h`.

The module defines its shape, fit and statistical callbacks. The host
checks the structure and the plausibility of a passing result. Validation
of the statistical assumptions belongs to that model. The external entry
point currently covers a single prespecified window.

The `test_peak_module` integration test builds the fixture in a separate
CMake project after the host exists. It loads the resulting library into
that executable, checks the stages and rejects an incompatible descriptor.
The test compares the host's hash and modification time before and after.
The fixture demonstrates the ABI and uses test-specific margins.

## Run a declared window family

```text
qAlgorithms --peak-model-module-search PATH STEP MIN_AREA MIN_COVERAGE ALPHA ROTATIONS MIN_SCALE MAX_SCALE [--quadratic MIN_AREA MIN_COVERAGE] I0 I1 I2 I3 I4 [...]
```

The command enumerates every fitting center and half-width in the supplied
range. Add `--quadratic` to test the same windows with the built-in
quadratic model as well. Its area threshold is in intensity times sample
index. The module threshold uses intensity times the declared coordinate unit.

In C++, pass `NativeCandidateModelWindows` entries to `CandidateSearchPlan`.
The plan owns shared references to the loaded models. The C++ API also
accepts both built-in models and explicit, model-specific window lists.
Each native descriptor name must be unique within the family. Results use
the descriptor name for native windows and an optional model enum for built-in
windows. The command prefixes window labels with `native:` or `builtin:`
and percent-encodes delimiters in native names.

Declare the complete model and window family, thresholds and random seeds
before inspecting the intensities. The raw-background stage uses the
same Bonferroni allocation for all model-window pairs, including skipped
windows. The rotation budget must resolve that threshold. The guarantee
uses the fixed-window raw Gaussian null and the rotation assumptions.
Module-specific statistical stages keep their own assumptions.

The integration test uses the separately built module in both a native-only
family and a mixed family. It checks that adding the built-in model doubles
the denominator and that windows containing zero remain counted.
