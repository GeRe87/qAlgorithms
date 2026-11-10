# Research-only hyperbolic running regression for qPeaks

**Status: experimental.** This isolated C++17 kernel is built only with
`-DQALGORITHMS_ENABLE_HYPERBOLIC_RESEARCH=ON`. Production qPeaks,
qCentroids, the existing inverse-cache header, feature selection, and the
experimental mzML CLI are unchanged.

## Linear model and closed integral

For aligned, equally spaced, positive log-intensity traces with a
**fixed externally specified** dimensionless coordinate step
`u=(sample-index-window-center)*step`, fit

```text
y_j(u) = a_j + b*u - C_j*cosh(u) - D_j*sinh(u)
       = a_j + b*u - p_j*exp(u) - q_j*exp(-u)

p_j=(C_j+D_j)/2, q_j=(C_j-D_j)/2
```

All unknown coefficients are linear in OLS. Each of the `m` traces has
three independent parameters `a_j,C_j,D_j`; `b` is shared and is the
apex-centered skew coefficient. There is **no fitted lambda**, shape
grid search, nonlinear optimization, fixed left/right breakpoint, or
per-window generic matrix inversion.

When `p_j>0, q_j>0` the model has one strict maximum, and its **full area
in u** is

```text
u_apex = asinh(b/(2*sqrt(p_j*q_j))) + 0.5*log(q_j/p_j)

A_u = 2*exp(a_j)*(q_j/p_j)^(b/2) * K_b(2*sqrt(p_j*q_j))
```

where `K_b` is the modified Bessel K-function. `logAreaU` uses the C++17
`std::cyl_bessel_k` implementation and returns NaN if the peak cannot
be interpreted as an integrable whole-line shape. This calculation is
**outside the running regression hot loop**. Conversion to a mass-domain
area requires the corresponding coordinate Jacobian. The Bessel
function may underflow for extreme coefficient magnitudes; a
production implementation needs a scaled-log Bessel path.

## Orthogonalized least squares: no large matrix inversions

For a symmetric sample grid `u=-h*step,...,h*step`, define

```text
n = 2h+1
cbar = sum(cosh(u))/n
Vc = sum((cosh(u)-cbar)^2)
Ss = sum(sinh(u)^2)
kappa = sum(u*sinh(u))/Ss
R = sum((u-kappa*sinh(u))^2)
```

Let `Y0_j=sum(y_j)`, `YC_j=sum(cosh(u)*y_j)`,
`YS_j=sum(sinh(u)*y_j)`, and `YU_j=sum(u*y_j)`.
Because even and odd regressors are orthogonal, exact pooled OLS is

```text
b   = sum_j(YU_j-kappa*YS_j)/(m*R)
C_j = -(YC_j-cbar*Y0_j)/Vc
a_j = Y0_j/n + C_j*cbar
D_j = kappa*b - YS_j/Ss
```

The cache precomputes `cbar,1/n,1/Vc,1/Ss,kappa,1/R` and endpoint
`cosh/sinh`. Residual norms are accumulated directly in extended
precision to reduce cancellation. Only four data moments per trace are
needed; when expanding a centered window from half-width `h-1` to
`h`, update them using the new left and right observations. Work is
**O(m) per candidate window**, not O(m*n), with O(m+maxHalfWidth)
working memory and no allocations in the inner loop.

A window is marked non-identifiable if the relative even/odd
residual information becomes tiny. The research routine still emits
its unconstrained coefficients for speed/diagnostics; such results
must not be considered accepted peaks.

## Complete covariance structure (homoskedastic OLS only)

Multiply the following factors by the estimated *pooled* log-domain
residual MSE after selection:

```text
Var(b)       = sigma²/(m*R)
Var(C_j)     = sigma²/Vc
Var(a_j)     = sigma²*(1/n + cbar²/Vc)
Cov(a_j,C_j) = sigma²*cbar/Vc
Var(D_j)     = sigma²*(1/Ss + kappa²/(m*R))
Cov(b,D_j)   = sigma²*kappa/(m*R)
Cov(D_j,D_k) = sigma²*kappa²/(m*R), j != k
```

Remaining even-versus-odd covariances are zero under the equal-grid,
equal-weight OLS assumptions. The covariance factors are checked
against the inverse of an independently assembled **full**
`(3m+1)x(3m+1)` normal matrix for `m=1,3,10`. The hot kernel
intentionally does not compute RSS or coefficient-level area
uncertainty; these must be added before quantitative deployment.
Nonlinear area transformations require appropriate propagation
(Delta method or dedicated uncertainty assessment), not a claim of
exact linear uncertainty.

## Build and use

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DQALGORITHMS_ENABLE_HYPERBOLIC_RESEARCH=ON
cmake --build build --target test_hyperbolic_running_research \
  bench_hyperbolic_running_research --parallel 2
ctest --test-dir build -R '^test_hyperbolic_running_research$' \
  --output-on-failure
./build/bench_hyperbolic_running_research 9 8
```

The benchmark uses **the actual production**
`qAlgorithms::forEachCoefficientsMultiTrace` as its baseline, not an
approximation or Python port. Both fit the same 3969 windows from
129-sample traces, 1/4/16 traces per case, and report median
nanoseconds/window (nine paired, alternating-order trials, eight
full runs per trial).

Two deterministic shapes with Gaussian log noise (sigma 0.015) are
compared: a Gaussian and a smooth asymmetric hyperbolic peak.
At a selected center/half-width 24, independent, known noise-free
truth is used to report:

- RMSE of fitted vs true log-intensities on identical samples
- Mean absolute apex error (dimensionless u)
- Mean absolute percent error of whole-line area in u
- Count of valid integrated candidate profiles

The **legacy baseline has m+3 coefficients** (one intercept per trace,
three shared shape coefficients). The new model has **3m+1
coefficients** (three per trace, shared b). Timing measures the real
difference in parameterization; neither equal fit complexity nor
general superiority is assumed. Both include callback consumption of
their entire coefficient sets; cache construction and peak-area
integration are outside the timed hot loops.

The benchmark is *microbenchmark evidence*, not a feature-detector,
centroiding validation, mass-accuracy or instrument resolution study.
Do not compare fitted area in u to ion area in m/z without the
appropriate physical coordinate conversion.

## Boundaries / scientific gates

- Matched, complete, equally spaced trace grids are mandatory for the
  cached orthogonal formula; nonuniform or missing positions require
  a separate weighted/generalized reference solver.
- `step` **must be supplied independently**, not tuned on the measured
  peaks. It determines the physical tail shape even though all
  regression coefficients remain linear.
- Large shared `b` alone does not enforce exactly translated and
  identical shapes; `C_j` and `D_j` may still vary between traces.
- Near apex, `u` and `sinh u` become nearly collinear. The pooled
  asymmetry uncertainty can be enormous for narrow windows. Never
  silently use these fits to report confident peak areas.
- `C_j <= abs(D_j)` makes at least one tail non-integrable and must be
  rejected for quantification; this is **not** an extra constrained
  regression step.
- Selection, cross-trace noise correlation, baseline, peak overlap,
  uncertainty in step/axis and mass calibration are not accounted for.
- Final comparison on independently measured, suitable profile-mode
  TOF/Orbitrap data and benchmarks across hardware are still required.

**No merge or production integration is authorized by this experiment.**
