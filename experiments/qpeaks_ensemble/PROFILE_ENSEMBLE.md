# Profile ensemble pilot

This experiment reuses the v3 slope-field basins and scale graph. It keeps
junction flags and combines native log profiles within each existing branch.
Weights are the squared parabolic taper of each native window. The quotient
derivative includes the derivative of the weights. Strict slope changes
bracketed on the scan grid give geometric maxima. Opposite crossings inside
one scan interval can be missed. No additional fit of the observations is
performed and no four-coefficient export is produced.

Shape diagnostics use linear intensities, scan coordinates and the chord
through the sampled support boundaries, or through the minimum sampled scan
between adjacent maxima in a branch. The area is the signed trapezoidal
excess above that chord and the asymmetry is the ratio of the two half widths
relative to the raw profile apex. These definitions can yield negative areas
or missing shape estimates for poorly bounded candidates. They are neither
quantitative RT areas nor statistical acceptances.

Conditional uncertainty uses scale-cached linear least-squares kernels and
the gradient of each metric with respect to the original log observations.
The same observation contributes jointly through all windows that contain it.
The derivative checks include the apex, intensity, chord area and width
ratio. Conditional errors assume fixed membership, integration boundaries and
crossing identity. The noise experiment uses known iid Gaussian log noise
with sigma 0.05. It uses double-precision reference kernels and eight scales.
It does not provide a calibrated covariance for real EICs or uncertainty
intervals after branch selection.

## Results

On 128 repetitions of each exact single-peak model, the oracle membership
gives one maximum in every case. The nominal 95% conditional intervals cover
apex in 121/128 symmetric and 122/128 asymmetric cases. Chord area is covered
in 119/128 and 124/128, width asymmetry in 121/128 and 117/128, and intensity
in 122/128 and 113/128. Automatic grouping gives exactly one proposal in
105/128 and 76/128. These are small conditional coverage experiments.

The native paired cohort contains 280 cases and 878,865 regressions.
The profile has exactly the reference count on 16/32 wide single cases,
32/32 narrow single cases, 32/32 triples and 1/32 overlapping pairs.
The v3 refit returns two proposals on 17/32 overlapping pairs and the
strict cascade two descriptions on 32/32. The profile proposes 154 and
4 maxima on pt_03 and pt_04, with the same counts on their reflected
controls. Neither legacy expectation is met. Debug and Release input
banks give identical results excluding timings. The Python profile
stage takes longer than the v3 refit stage in this pilot.

On ten fixed real EICs the profile gives 17,694 geometric maxima,
the refit 5,163 and native selection 4,195. The cascade with collectives
enabled gives 293 descriptions and succeeds on all ten. Strict cascade
fails explicitly on eight. On an additional length-stratified sample of
32 EICs per polarity the profile gives 4,814 positive and 1,707 negative
proposals. The refit gives 1,197 and 529, native selection 1,149 and
478, and the cascade with collectives enabled 115 and 67 descriptions.
The sample selects eight EICs from each length quartile within each
polarity. Quartile boundaries use integer fractions of the complete
inventory sorted by length and bin ID. Sampling uses Python Random seed
698093. Real EICs have no reference labels. The results do not establish
accuracy, or provide complete feature lists for the mzML.

## Reproduction

Run against matching paired inputs, native fields and the two archived v3
comparisons. The comparator verifies the hashes before pairing cases.

```text
python experiments/qpeaks_ensemble/compare_profile_ensemble.py INPUT NATIVE OUT REFIT/case-results.json REPRESENTATIVE/case-results.json
python experiments/qpeaks_ensemble/study_profile_uncertainty.py UNCERTAINTY.json
python experiments/qpeaks_ensemble/test_profile_ensemble.py
```

NumPy and SciPy are required. No production code is used for the added
profile evaluation, and production detector behavior is unchanged.

## Apex corridor and valley audit

The audit replays the profile outputs against the exact same field bank.
Every maximum is classified by whether its apex belongs to a native root
corridor of its own branch, of another branch on the same positive block,
between its own corridors, or outside all corridors. The test uses closed
corridor bounds without a distance tolerance. A root corridor is a geometric
construct and can move under noise or differ from the ensemble apex.

Of 17,694 maxima in the fixed real diagnostics, 4,515 lie in their own
branch's corridors and 12,841 in another branch's corridor. In the
64-EIC stratified sample the corresponding counts are 2,081 and 4,174
out of 6,521 maxima. No real branch has more than one maximum in its own
corridors. A hypothetical hard own-corridor rule gives the expected count
in 29/32 narrow cases, 29/32 triples and 7/32 overlap cases, compared
with 32/32, 32/32 and 1/32 before restriction. It retains 78 maxima
for pt_03 and three for pt_04. It is an attribution diagnostic only.

The valley diagnostic measures both log-height drops to the shared sampled
separator between adjacent maxima of a single branch. The analytic standard
error under known iid log noise is computed from the difference between the
apex and separator influence vectors on the same original scans. The
selected locations and membership are held fixed for propagation. In 64
trials per separation, both drops exceed 1.96 conditional standard errors
in 0, 0, 1 and 56 cases at separations of 0, 4, 6 and 9 scans for equal
Gaussians of width 2.5 scans. The first two true profiles are unimodal and
the last two bimodal. The pointwise contrast is exploratory after searching
for peaks and separators. This audit does not measure a valley between
different branches, and no noise standard errors are assigned to real EICs.

```text
python experiments/qpeaks_ensemble/audit_profile_corridors.py NATIVE PROFILE AUDIT
python experiments/qpeaks_ensemble/study_valley_contrasts.py VALLEY.json
```
