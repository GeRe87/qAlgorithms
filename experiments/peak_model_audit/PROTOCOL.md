# Fixed real-profile audit protocol

Base: 03649958394fecc759de802ec429d4ecd5f27d95. Experimental audit only.
Two public TripleTOF serum pool injections from Zenodo 20729183, with published
MD5 verification. Same analytical run, not certified standards or independent
instruments. No known masses/areas/noise truth supplied.

80 positive profile MS1 scans, evenly stratified over scan order. Per scan,
64 fixed index-bin centers (or all available centers if fewer), windows of
9 and 17 observed samples. Separately: up to 64 evenly subsampled observed
local maxima >=100, using the same widths. The maxima cohort is data-selected.
Fixed index windows are chosen before reading intensities, but intensity/spacing
eligibility still censors the analyzed subset; all attempted windows are recorded.
No null-test or physical accuracy claim is based on this subset.

Fixed u step .35, matching prior research. No optimized step, interpolation,
zero filling or coordinate/Jacobian calibration. Eligibility: positive finite
raw samples and increasing sqrt(m/z), local gap/global lower-third gap ratio
in [.55,1.45], local gap CV<=.20. Rejection reasons are retained, not filtered
out of denominators. Full scans are also passed to the new search in default
exploratory mode, recording complete-plan input rejection versus evaluation;
no prespecified statistical acceptance is enabled on observed data.

Independent QR fits compare every retained kernel coefficient/prediction.
Residual lag-1 and Durbin-Watson-style ratio are compared descriptively to
10000 projected iid Gaussian residuals R*epsilon per geometry, seed 617281,
R=I-QQ^T. Residuals are correlated even for independent errors because of the
fit. No unadjusted zero-correlation baseline. Reference quantiles are not
selection-adjusted p-values and do not diagnose raw noise versus model mismatch.
Record zero/low residual energy, raw-refit validity, integrability/support,
window-width sensitivity and unsupported properties. No production routing.
