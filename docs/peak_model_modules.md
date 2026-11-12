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
