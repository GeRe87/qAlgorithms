# Research: closed-form log-Gamma peak folding with running regression

**Status: research only.** Isolated Python proof-of-concept based on the Draft
mzML folding branch. It neither replaces nor modifies qCentroids, qPeaks,
feature detection, or the existing symmetric/asymmetric folding solver.

## 1. Model and exact translation

In a dimensionless local acquisition coordinate `z`, for an orientation
`s ∈ {-1, +1}` and a **fixed** `lambda > 0`:

```text
log I_j(z) = a_j + d*z + c_j*exp(s*lambda*z)
                c_j < 0,  s*d > 0
```

The coefficient `d` is pooled over profiles, while `a_j, c_j` are local.
For fixed `lambda`, all unknown coefficients are linear in log-intensity
OLS. The same `d,lambda,s` imply **identical shape under exact translations**:

```text
p = d/(s*lambda) > 0
log I_j(z) = log H_j + p*(s*lambda*(z-mu_j) + 1
                         - exp(s*lambda*(z-mu_j)))
```

with local apex `mu_j` and height `H_j`. This is not a
piecewise-quadratic function and has no breakpoint.

The peak is strictly log-concave, with a unique maximum:

```text
z_apex,j = log(-d/(s*lambda*c_j)) / (s*lambda)
(log I)''(z_apex,j) = -s*d*lambda < 0
```

## 2. Closed area and coordinate units

The full positive line-shape integrates analytically in `z`:

```text
A_z,j = exp(a_j) * Gamma(p) / (lambda * (-c_j)^p),
p = d/(s*lambda) > 0
```

Compute `log(A_z)` using `lgamma` to avoid overflow. The implementation
returns log-area and does not exponentiate it automatically. A fixed physical
linear scale `x = x0 + h*z` multiplies the area by `h`; a nonlinear
transformation from native instrument coordinates to m/z requires the
appropriate Jacobian. **This area must not be reported as an m/z area
without a justified coordinate conversion.**

A useful special case is `z = log(m/m_ref)` and `s=+1`: the exact
**m/z-domain** area is

```text
A_m = m_ref * exp(a_j) * Gamma((d+1)/lambda)
      / (lambda * (-c_j)^((d+1)/lambda))
```

provided the full positive-mass model and its units are appropriate. This
multiplicative mass-coordinate model is generally *not* the same as a
constant-width translation in TOF or Orbitrap native sampling coordinates.

## 3. Pooled OLS with cached 2x2 inverses

For each profile construct columns `[1, E]` with `E=exp(s*lambda*z)`.
Let `G=L^T L` (2x2), `h=L^T z`, `g=L^T y`.

```text
S = z^T z - h^T G^-1 h
r = z^T y - h^T G^-1 g
d = sum_j r_j / sum_j S_j
[a_j, c_j]^T = G_j^-1 (g_j - h_j*d)
```

The inverse `G^-1`, `h`, and Schur denominator `S` depend only on
the coordinate geometry, lambda and orientation: they can be cached.
For identical uniform geometries the cached design occupies O(window length) memory and per-trace moments occupy O(number of profiles). Distinct nonuniform geometries require their own basis arrays (up to O(total samples)). The Schur solve is scalar.
Different sampling patterns and profile lengths remain valid but require
corresponding geometries.

## 4. True constant-time running sums on a uniform grid

For a window with fixed offsets `z_i = (i-(n-1)/2)*h` maintain:

```text
Y0 = sum y_i
Yz = sum z_i*y_i
Ye = sum exp(s*lambda*z_i)*y_i
Yy = sum y_i^2
```

All three regression RHS moments and the sum of squared observations
update in O(1) per shifted window; neither a large inverse nor a
per-window rescan is needed.

**Numerical stability:** right-shifts are preferable for `s=+1`,
left-shifts for `s=-1`. Opposite traversal gives an exponential
recursion multiplier greater than one and can amplify roundoff over many
windows. Long streaming runs should also periodically recompute sums
from scratch to control accumulated error. The `Geometry` class provides
the cached precomputation and the two update directions.

For nonuniform native coordinates, the fixed 2x2 inverse is generally
not reusable. An O(1) update of the *small Gram matrix and RHS moments*
using added/removed samples is still possible in a fixed coordinate frame,
but matrix factorizations and recentering require separate implementation.

## 5. Reproducible verification

Run from this folder:

```powershell
python -m unittest discover -s . -p "test_*.py" -v
```

The seven deterministic tests verify:
- Exact apex positions and Gamma total areas of translated peaks, both orientations
- Constant-time rolling moments against direct recomputation over 318 shifts
- Multiple profile lengths in one shared fit
- Nonuniform profiles and rejection of irregular-grid rolling updates
- Degenerate design and nonintegrable-fit rejection
- Independently assembled full dense OLS on noisy nonuniform profiles (both orientations)
- Analytical area versus independent Simpson quadrature (both orientations)

Locally all seven tests passed with Python 3. The tests are a proof of
mathematical/numerical feasibility, not scientific mass-calibration
or instrument validation. In particular, the exact-response tests are
deliberately generated from the assumed model.

## 6. Scientific restrictions and next decisions

- The **shape parameter lambda is fixed externally**. Selecting lambda
  from a discrete grid means multiple linear fits and model selection,
  not a single globally linear optimization.
- Shape invariance requires a shared `d` **and** the same dimensionless
  coordinate scale. If native peak widths change with m/z, group similar
  resolution regions or use externally determined scale normalization.
- Unconstrained OLS can yield `s*d <= 0` or `c_j >= 0`.
  This prototype correctly rejects such cases; a constrained solver
  is not implemented.
- The log-domain model assumes positive intensities; selection bias,
  baseline, censoring, dependent samples, peak overlap, heteroscedasticity,
  and mass-calibration uncertainties are not addressed.
- The left/right tail decays are intrinsically different. The model is
  exact for its own family but can bias apex locations for Gaussian or
  other near-symmetric experimental peak shapes.
- Comparison with the smooth log-cosh/Beta and hyperbolic/Bessel families
  is a separate required validation track. Both have closed total
  integrals but do not provide the same unconstrained, exactly pooled
  translation structure.
- No real QTOF mass accuracy claims are made.

The existing experimental mzML route and its reference Schur solver are
unchanged. No production or merge approval is implied.
