# Experimental linear peak folding (qPeaks)

This is a **coefficient-level research API**, deliberately separate from
the existing peak-picking / feature-identification pipeline. The implementation
lives in the existing qPeaks files:
- header: `include/qalgorithms_qpeaks.h`
- source: `src/qalgorithms_qpeaks.cpp`
- tests: `test/test_qpeaks_folding.cpp`

No changes are made to qCentroids, the mzML reader, the previous
`findCoefficients`/multi-trace regression, or the six-value matrix cache.

## Mathematical model

For preselected, measured log-intensity profiles j=1..m, at possibly
irregular local coordinates u_ji, fit

```
log I_j(u) = a_j + d_j*u + qL*u^2 for u < 0
log I_j(u) = a_j + d_j*u + qR*u^2 for u >= 0
```

The a_j and d_j are specific to the peak, while qL/qR are common across
all selected profiles. This is **linear in every coefficient**. In symmetric
mode constrain qL=qR=q, giving exactly translated quadratic peaks with
individual d_j and a_j. In asymmetric mode the curvature breakpoint is
**at u=0 for all profiles**: allowing individual slopes DOES NOT shift
this breakpoint. Thus a shifted asymmetric shape is only approximate
unless the coordinate origins have first been aligned to an externally
determined breakpoint.

The local u values should be in instrument sampling-bin units (typically
a calibrated transform such as sqrt(m/z) for TOF), with a local anchor
near each candidate peak; the peak-specific transformations are supplied
by the caller and **must use the same physical scale**. Arbitrary grids
and unequal numbers of available measurements per peak are supported.
There is no need to interpolate the intensities. Omitted measurements
are not silently converted into zeros.

Input is `std::vector<FoldingProfileView>`, each containing a pointer to
double-precision u, a pointer to finite log intensities and the sample
count. The output includes local intercepts and slopes, shared curvature(s),
full conditional covariance factors, residual sums, degrees of freedom,
and conditional estimated apex and its standard error if identifiable
and within the measured coordinate window.

## Exact Schur-complement computation

For each profile, let Z_j=[1,u] and let C_j be either the single column
u^2 (symmetric) or two disjoint columns u^2 * (u<0),
u^2 * (u>=0) (asymmetric). Denote

```
H_j = Z_j^T Z_j              [2x2]
K_j = Z_j^T C_j              [2xq]
G   = sum(C_j^T C_j)         [qxq]
r_j = Z_j^T y_j             [2]
s   = sum(C_j^T y_j)         [q]

S = G - sum(K_j^T H_j^-1 K_j)
t = s - sum(K_j^T H_j^-1 r_j)

beta = S^-1 t                  [q=1 or q=2]
alpha_j = H_j^-1 (r_j-K_j beta)
```

Only per-profile 2x2 matrices and one final 1x1 or 2x2 need inversion.
The runtime uses O(m) scratch space beyond the input and output; unlike
a full inverse it does not grow quadratically with m. Per-peak moment
accumulation is O(number of supplied points).

The complete inverse normal matrix blocks are available analytically:

```
Cov(beta) / sigma^2 = S^-1
Cov(alpha_j,beta) / sigma^2 = -H_j^-1 K_j S^-1
Cov(alpha_j) / sigma^2 =
    H_j^-1 + H_j^-1 K_j S^-1 K_j^T H_j^-1
```

All returned `*Factor` covariance fields must be multiplied by the pooled
residual MSE to obtain the conditional variances/covariances. There are
N-(2m+q) residual degrees of freedom, where N is the total number of
supplied observations and q is 1 (symmetric) or 2 (asymmetric).

Only when BOTH fitted half-curvatures are negative and the apex lies
inside the observed local coordinate range do we report an isolated
unimodal apex, x_apex=-d_j/(2*q_side). Its conditional standard
error uses the delta method and the **full covariance including
Cov(d_j,q_side)**:

```
dx/dd = -1/(2*q_side)
dx/dq = d_j/(2*q_side^2)
Var(x_apex) ~ (dx/dd)^2 Var(d_j)
            + (dx/dq)^2 Var(q_side)
            + 2 (dx/dd)(dx/dq) Cov(d_j,q_side)
```

In the asymmetric model, zero slope corresponds to a side switch;
the program does not claim a differentiable apex uncertainty there.

If the u-axis is defined as
`u=(sqrt(m/z)-sqrt(m/z_anchor))/samplingStep` for a calibrated TOF,
the inferred mass is `(sqrt(m/z_anchor)+samplingStep*u_apex)^2`.
The associated local mass standard error, **conditional on calibration**,
follows by multiplication with
`2*samplingStep*(sqrt(m/z_anchor)+samplingStep*u_apex)`.
Do not assume this transform without checking the instrument's actual
sampling law and missing-bin representation.

## Scientific boundaries

- This is *unweighted* least squares on **log intensities**. Zero/negative
  raw intensities must be excluded before the log transformation.
- Homoscedastic, uncorrelated log residuals and fixed coordinate anchors
  are assumptions of the reported classical covariance; correlated
  Fourier-profile bins require generalized/error-aware treatment.
- Profile selection, overlap, shared-peak-form hypothesis, ion identity
  and instrument calibration are **not** validated by this numerical kernel.
- The local apex covariance does NOT propagate uncertain selection,
  alignment, anchor positions or deviations from a common peak shape.
  The method should not be reported as a calibrated m/z uncertainty
  without external mass standards and independent replicate validation.
- A sum of unresolved physical ions is not a sum of log peaks. Deconvolution
  is outside scope, and overlapping peaks should be detected/excluded from
  training rather than silently absorbed.
- Per-profile individual slopes and the asymmetric breakpoint can jointly
  lead to a biased mass estimate when origins are misaligned. The symmetric
  version is the correct exact-translation baseline.

## Tests and integration status

`test/test_qpeaks_folding.cpp` constructs the full design/normal matrix
independently for single and multiple (1,2,4,10) profiles with irregular
u grids, compares fitted coefficients and every reported covariance
factor to a full inverse, and checks delta-method apex uncertainty
against the corresponding full covariance. It also verifies exact
symmetric peak translation, non-increase in conditional slope variance
factor with pooling, input validation and rank-deficiency handling.

Neither the previous seven-value multi-trace inverse nor the existing
six-value single-trace cache is required to run this irregular-grid
fitter. The API is intended for an experimental, explicit profile
selection/coordinate preprocessing path in future work. No automatic
replacement of qPeaks/qCentroids is performed.
