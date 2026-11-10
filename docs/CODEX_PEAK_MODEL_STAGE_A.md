# Codex implementation briefing: qPeaks modular PeakModel contract (Stage A)

**Repository:** `GeRe87/qAlgorithms`  
**Your working branch:** `agent/peak-model-contract-stage-a`  
**Immutable starting point (before this briefing commit):** `8c6b51087205e11d9c09296e74b8a1f8ecf5ea59`  
**Pull request target:** `research/hyperbolic-running-regression` (stacked on Draft PR #4).  
**State:** only a briefing, no implementation yet.

## Objective

Refactor the original piecewise-quadratic qPeaks implementation behind a **minimal model interface**, preserving existing numerical results and performance. The model selector should ultimately support `PiecewiseQuadratic` and `Hyperbolic`, but **Stage A must only route the existing piecewise model**. An attempt to use hyperbolic in the production pipeline during Stage A must fail explicitly, not silently select the legacy model.

Do not directly modify `main`, `work/mzml-featurelists`, `feat/multitrace-shared-shape`, `feat/linear-peak-folding`, `feat/peak-folding-mzml-cli`, or `research/hyperbolic-running-regression`. Do not edit or merge PRs #1–#4. Work only on your isolated branch, open a **Draft** PR after implementation.

## Existing code / assumptions to inspect

- `src/qalgorithms_qpeaks.cpp`: `findCoefficients`, `regAt`, `regExp_fac`, `correctB0`, `calcRegSpan`, `validateRegressions`, `calcRegressionProperties`, `insignificantCoefficient`, `peakPosition`, `peakPositionUncert`, `peakHeightUncert`, `fullWidthHalfMax`, `peakArea`, `matProductReg`, `f_testRegression`, `findFeatures`.
- `include/qalgorithms_qpeaks.h`: `qpeaks_find`, the original and multitrace running regression APIs, and experimental folding.
- `include/qalgorithms_datatypes.h`: `RegCoeffs`, `RegressionGauss`, and `FeaturePeak` contain old-model-specific representation.
- The existing `qalgo_matInverse` array of six precomputed coefficients, linked into legacy running regression.
- Existing old-model centroiding and feature outputs, including uncertainty calculations and the specific raw-intensity `correctB0` stage.
- `experiments/hyperbolic_running/`: separate research-only C++ kernel, checks and benchmark from Draft PR #4. **Do not integrate into the production pipeline in Stage A**.

## Required architecture

1. Add `enum class PeakModelKind { PiecewiseQuadratic, Hyperbolic }` and a small `PeakFitOptions`, defaulting to `PiecewiseQuadratic`. Preserve the existing `qpeaks_find(...)` public symbol/semantics; expose opt-in selection via a compatible overload or entrypoint if reasonable.
2. Define a statically dispatched, strongly typed model contract (templates/traits or equivalent) with no virtual dispatch and no per-window allocation in the hot running kernel. **Select model once outside the inner loop**.
3. Separate the **model-specific math** from common orchestration where safe: running fit / cached geometry; log-value; amplitude-only intensity refit; apex; FWHM; integrability/shape constraints; closed analytic integral; covariance and property gradients. Some operations may be conditional and model-specific. Do not pretend all models have all properties.
4. Separate **shape validity** (e.g. interior maximum and integrability) from model-independent **candidate and peak validation** (observed support, statistical significance, baseline/noise assumptions, ambiguity/interference). Old F/t tests and degrees of freedom must not be mechanically reused across different models.
5. Introduce an internal or adapter-style model-neutral `PeakProperties` containing position, height, width, area, uncertainty, valid/invalid status and units. **Do not** replace existing `RegCoeffs`, `RegressionGauss`, `FeaturePeak` data layouts or persisted output formats in Stage A. Retain the legacy types through an explicit compatibility adapter.
6. Preserve all existing piecewise-quadratic equations and numerical evaluation order as much as possible. Avoid speculative cleanups. Keep the existing six-value inverse cache, raw-intensity amplitude postfit and legacy gradient calculations functioning.
7. Make the default/explicit `PiecewiseQuadratic` paths equivalent. On selecting `Hyperbolic` via the new production API during Stage A, produce a clear unsupported-model status/exception **before processing or writing outputs**.

## Scientific contract (for the separate Stage B / Stage C PRs)

The research model is `ln(I_j(u)) = a_j + b*u - C_j*cosh(u) - D_j*sinh(u)`, with a linearly fitted pooled `b` and individual `a_j,C_j,D_j`. It has precomputed geometry and running moments, a scalar raw-height correction, analytic apex and closed whole-line Bessel-K area in **dimensionless u** under integrable-tail constraints `C_j>|D_j|`.

Stage B will bind this kernel via `PeakModelKind::Hyperbolic`. Stage C will build model-specific statistical gates and full post-height-refit area-uncertainty propagation including shared coefficient covariance, and handle physical coordinate/Jacobian interpretation. **Neither is authorized for implementation in Stage A.** The existing Orbitrap/TripleTOF research reports do not establish certified reference-area accuracy.

## Acceptance checks (mandatory)

- Record original base SHA and baseline Debug/Release CMake/CTest results **before modifying production code**.
- Exact same output and ordering for legacy running coefficients and candidate selection; equal acceptance/rejection decisions, spans, positions, height after `correctB0`, areas, FWHM, uncertainties, centroid/EIC/feature outputs. Write differential tests against the baseline source/fixtures, documenting any tiny numerical tolerance.
- Default API and explicit `PiecewiseQuadratic` yield equivalent results; `Hyperbolic` fails safely and explicitly in Stage A.
- Preserve existing CLI behavior, file column names, persistence, exit codes and output naming.
- Include Debug and Release CI evidence; identify the **three documented existing legacy qPeaks test failures** separately rather than masking or quietly modifying expectations. Do not claim an entire test suite passes if baseline failures remain.
- Benchmark pre-/post-refactor running fits for equal windows/traces and same hardware/compiler; present median, variability and relevant compilation flags; any slowdown must be characterized and justified.
- Do not accidentally introduce a new runtime virtual call, repeated dynamic allocation, or matrix inversion in the inner running loop.
- No merge, release, destructive action or production approval. Publish a Draft PR with exact base/head SHA, changed files, evidence and open questions.

## Scope management

Prioritize a **small, reviewable vertical slice**. A clean contract with the legacy model working and tested is more important than ambitious simultaneous rewrites. If a full extraction would alter output or blow up scope, stop at an interface + adapters + differential tests and record which coupling remains.

## Agent handoff

1. Inspect repository HEAD, active stacked PRs and code paths; avoid disturbing concurrent worktrees.
2. Build and test baseline.
3. Implement Stage A on `agent/peak-model-contract-stage-a` (or a new Codex work branch based on its HEAD).
4. Verify against baseline and benchmark.
5. Open a Draft PR **targeting `research/hyperbolic-running-regression`**; include an architecture note, test evidence, regression differential summary, and Stage B/C plan.
6. Return the PR URL, source/head SHAs and unresolved blockers. **No merge**.
