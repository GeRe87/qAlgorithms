# qPeaks statistical contract and analytic audit

Status: a conditional mathematical reference and an audit of the current code,
not a claim that the real measurement follows this noise model. The retained
pt_03/pt_04 expectations are controls, not a basis for choosing thresholds.
Equidistant fitting and the precomputed inverse matrices are retained.

## 1. State the random model before assigning a test distribution

For a **fixed, predetermined** symmetric window x=-s,...,s, define

    z_i = log(y_i) = X_i beta + epsilon_i,
    X_i = [1, x_i, x_i^2 1(x_i<0), x_i^2 1(x_i>0)].

The exact reference considered here assumes positive observations and
independent epsilon_i ~ N(0,sigma^2), with common unknown variance, fixed design,
and full-rank X. These are assumptions about errors in log intensity; taking
logs does not establish them. Additive raw-intensity noise, clipping, a fitted
background, serially correlated scans and aggregated/interpolated observations
need their own derivation. Counting independent points alone does not repair
correlations or unequal variances. Index equidistance is the declared working
model; irregular physical RT and quantitative areas remain separate limitations.

Under this reference, with n=2s+1 and four coefficients,

    beta_hat = M X' z, M=(X'X)^(-1), nu=n-4,
    RSS = ||z-X beta_hat||^2, sigma_hat^2=RSS/nu,
    Cov(beta_hat)=sigma^2 M.

For any **fixed** contrast c, (c'beta_hat-c'beta)/(sigma_hat sqrt(c'Mc))
has Student t_nu distribution. The log-linear null is beta_2=beta_3=0.
For both least-squares models on the identical fixed window,

    F = ((RSS_linear-RSS_full)/2) / (RSS_full/(n-4)) ~ F_(2,n-4)

under that null. This tests improvement over a log-linear trend, not correctness
of the peak shape, absence of instrument noise or chemical peak identity.
Failure to reject does not establish that no physical peak exists. A good fit
and a sufficiently informative fit are different questions.

The assumptions and residual-variance convention agree with the primary
[NIST least-squares reference](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd431.htm)
and [NIST assumptions discussion](https://www.itl.nist.gov/div898/handbook/pri/section2/pri245.htm).
The following matrix and selection arguments are derived directly for qPeaks.

## 2. Why the current trimmed-window residual is not a standard F denominator

Let H=X M X' be the full-window projection, R=I-H, and S a diagonal indicator
selecting m shape-window observations. Even if S were fixed in advance, the
residual from full-window coefficients evaluated on S is

    RSS_shape = z' Q z, Q=R S R.

Under a correctly specified full model,

    E(RSS_shape)=sigma^2 tr(Q)=sigma^2 [m-tr(SH)].

In general Q is not idempotent. Its quadratic form is a weighted sum of
independent chi-square variables, not chi-square_(m-4). Replacing m-4 with
tr(Q) matches a mean but does not restore the exact chi-square law.
For fixed S the full coefficient estimate remains independent of Rz; that
independence alone does not turn this residual into a Student-t denominator.
In production S is chosen from the fitted coefficients, adding selection.

For the linear fit on the shape window let Q0 be its residual projection,
embedded in the full window. The numerator is z'(Q0-Q)z. Q0-Q need not be a
projection or positive semidefinite; its product with Q need not vanish.
Consequently the usual two independent chi-square construction of F fails.
The guard against negative improvement prevents an invalid rejection; it
cannot supply the missing null distribution.

Independent QR calculations verify the expected full-window projector
identities and expose the trimmed-window violations:

| Full n / selected m | Current denominator df | E(RSS_shape)/sigma^2 | Smallest eigenvalue of Q0-Q |
| --- | ---: | ---: | ---: |
| 5 / 5 | 1 | 1 | numerical zero |
| 7 / 6 (x=-3..2) | 2 | 2.9141604 | -0.3422057 |
| 9 / 7 (x=-2..4) | 3 | 4.1356305 | -0.2223055 |
| 17 / 12 (x=-8..3) | 8 | 9.3828936 | -0.3440974 |

The first two trimmed supports occur in the retained pt_04 candidates.
The earlier private full-domain experiment restores this necessary fit-domain
consistency, but does not by itself validate selection or resolve both controls.

## 3. Why the data-selected curvature is not a prespecified one-sided t test

`insignificantCoefficient` chooses beta_2 when beta_hat_1<0 and beta_3 otherwise.
`calcRegSpan` requires the chosen apex-bearing curvature to be negative. These
are functions of the same data, and the precomputed covariance contains
Cov(beta_hat_1,beta_hat_2)=sigma^2 D>0 and
Cov(beta_hat_1,beta_hat_3)=-sigma^2 D. Choosing a side is therefore informative
about its estimated curvature; it is not an independent choice.

A counterexample already exists under the most favourable full-window,
homoscedastic Gaussian reference, with all true slope/curvatures zero.
Normalize the three coefficient estimates to (Z1,Z2,Z3), and write
rho=D/sqrt(C E)>0, U~chi-square_nu independently. With tcrit=t_(nu,1-alpha),
the **isolated signed curvature gate**, without the other cascade filters,
rejects when the selected negative curvature exceeds tcrit. By reflection,
its rejection probability is

    2 E_U integral_[tcrit sqrt(U/nu)]^infinity
        phi(z) Phi(rho z/sqrt(1-rho^2)) dz.

For z>0 the Phi term is strictly greater than 1/2 and less than 1.
Thus alpha < P(reject) < 2 alpha. This is a derivation, not a fitted
correction factor. It does not assert this rate for the complete cascade.
The F gate and the geometry/selection stages can restrict the final event.

Independent simulations use the exact normal coefficient law from an
independently inverted design Gram matrix, and independent chi-square noise.
There are one million draws per scale, with recorded reproducible seeds:

| Scale / nu | Fixed left negative-curvature rate | Data-selected side rate |
| --- | ---: | ---: |
| 2 / 1 | 0.050049 | 0.096377 |
| 3 / 3 | 0.050340 | 0.099383 |
| 8 / 13 | 0.050130 | 0.100043 |
| 25 / 47 | 0.050034 | 0.100227 |

Monte Carlo standard errors for the selected rates are about 0.00030;
estimates can fluctuate slightly above the analytic upper bound 0.10.
These simulations isolate a distributional claim, not instrument noise,
conditioning on geometry or real-file false-discovery rates.
Changing the side rule or doubling a critical probability without specifying
the scientific null is not adopted as a production correction.

## 4. A possible analytically justified peak-shape test, not yet adopted

For a predetermined support [a,b] with a<0<b, the endpoint derivatives are

    dL=beta_1+2 a beta_2, dR=beta_1+2 b beta_3.

Within the current continuously differentiable piecewise-quadratic log model,
dL>0 and dR<0 imply exactly one interior maximum on that support: each half's
derivative is linear, and the derivative is continuous at zero. The case of
beta_1>=0 puts the crossing on the right; beta_1<0 puts it on the left.
This permits a convex non-apex half if its valley lies outside the support.

The null 'dL<=0 OR dR>=0' is a union. Reject it only if **both** fixed one-sided
contrast tests reject at alpha. This intersection-union test has size at most
alpha under the stated fixed-window model, regardless of correlation between
the contrasts: at least one true component null bounds the probability of the
joint rejection. Equivalently p_peak=max(p_L,p_R) is a valid conservative
p-value for that declared null. No invented amplitude ratio, width threshold
or additional matrix inversion is required; contrast variances use c'Mc.

Important limits: this tests existence of a peak on the declared support within
the model, not general lack of fit. A data-chosen support and the multiscale
search still require a selection/family argument. Requiring signs at the full
window edges would reject some currently permitted fits whose convex tail
contains a valley inside that window. This is a real change of the scientific
target, so it is a candidate derivation rather than an unannounced replacement.
For the first pt_04 five-point fit, the endpoint evidence is asymmetric:
approximately t_left=0.5096 and t_right=13.1715 against t_(1,.95)=6.3138.
The ordinary full-window F statistic is 102.4836 against 199.5 (p about .0697).
A visually apparent peak need not be established by these five observations
under the stated four-coefficient, unknown-variance reference. These values
are not a justification to lower alpha or alter the retained expected count.

## 5. A dimensionally invalid area gate: corrected independently

The production check `area<=e` compared a dimensional quantity (intensity times
axis units) with a unitless constant. It is not a statistical test. For the
same stored intensities y_i=10 exp(-(i-40)^2/18), i=0..80, the original code
returned one peak with dx=0.1 seconds (area about 7.51988), but zero with
dx=0.1/60 minutes. This changes representation of the same physical axis.

This continuation removes that numeric cutoff. It retains a finite, strictly
positive representability requirement. It does not replace it with another
height/area significance threshold or change the area formula. The downstream
EIC builder's `area>1` assumption is replaced with an explicit finite-positive
runtime check, active in Release as well as Debug. This accepts valid subunit
areas while rejecting zero, negative, infinite and NaN values explicitly.
New tests preserve the same count and fitted coefficients in different axis
units, check analytic Gaussian area/apex scaling, retain subunit areas through
EIC construction, and enforce rejection of invalid values. Original test
assertions and pt_03/pt_04 expectations are unchanged.

This does **not** remove the separate input segmentation at intensity>e.
That floor is also tied to numeric units and is not a calibrated detection
limit. Its interaction with log-domain observations and sample selection must
be resolved separately. No claim of general intensity-unit invariance follows
from correcting the area gate. No interpolation or quantitative RT/area
correction is introduced here.

## 6. Contract for subsequent cascade changes

Every acceptance stage must state: (a) physical/mathematical question,
(b) null and alternative, (c) estimator and exact data support, (d) assumptions,
(e) statistic and reference law, (f) scope of its error guarantee and
(g) treatment of nonrepresentable or unsupported cases.

The current layers must not be conflated:

| Current layer | Actual role / outstanding requirement |
| --- | --- |
| Input intensity>e | Uncalibrated numerical segmentation; requires separate treatment |
| Geometric apex/span checks | Model/support eligibility; not p-values |
| Chosen-side t gate | Selection-dependent; fixed-coefficient interpretation invalid |
| Log-linear F comparison | Exact only under compatible fixed-domain fitting and noise assumptions |
| Raw-intensity b0 correction | Different estimator; does not calibrate log covariance or property coverage |
| Area positivity/finiteness | Numeric/model representability, no significance claim |
| Solitary-fit height versus adjacent samples | Existing geometric/data-dependent heuristic, not a statistical test |
| Apex grouping / collective fallback | Representation of accepted fits, not new statistical evidence |

A chosen alpha is a declared error tolerance, not a measured physical constant.
It may be configurable and needs an explicit scope; changing it to satisfy a
fixture is not allowed. Applying several filters to one window does not imply
independence or multiply their error probabilities. Searching many windows
does not inherit a local alpha guarantee. Bonferroni is justified only after
valid per-window p-values and a defined family exist; its sparse-peak power
cost in the preceding paired experiment remains relevant.

The next method implementation should first establish a coherent fixed-window
reference while retaining the precomputed equidistant design; then define
support/search selection and the physical noise interpretation. The original
control outcomes remain visible throughout. The whole cascade is not yet
statistically validated, and the two retained count failures remain unresolved.
