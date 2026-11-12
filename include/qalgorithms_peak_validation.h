#ifndef QALGORITHMS_PEAK_VALIDATION_H
#define QALGORITHMS_PEAK_VALIDATION_H

#include "qalgorithms_peak_fit.h"
#include <optional>
#include <cstdint>

namespace qAlgorithms
{
    namespace peak_model
    {
        enum class PeakValidationStatus { notEvaluated, failed, passed, unavailable };

        struct PeakGeometryValidation
        {
            PeakValidationStatus form = PeakValidationStatus::notEvaluated;
            PeakValidationStatus apexSupport = PeakValidationStatus::notEvaluated;
            double apex = std::numeric_limits<double>::quiet_NaN();
        };

        struct PeakFormNoise
        {
            double logMse = 0;
            size_t residualDf = 0;
            double alpha = .05;
        };

        struct PeakWindowCoverage
        {
            PeakValidationStatus status = PeakValidationStatus::unavailable;
            double fraction = std::numeric_limits<double>::quiet_NaN();
            double standardError = std::numeric_limits<double>::quiet_NaN();
            double lowerConfidenceBound = std::numeric_limits<double>::quiet_NaN();
        };

        // Form uses shape coefficients alone. Quadratic mixed tails refer to
        // the clipped peak span; a flat or rising flank at b1=0 has no strict apex.
        bool hasPeakShape(const RegCoeffs &coefficients);
        bool hasPeakShape(const HyperbolicModel::Coefficients &coefficients, double sharedSlope);

        // Apex is relative to the center, in samples for quadratic and in u
        // for hyperbolic. The support test counts strictly bracketing samples
        // on the original regression grid, before any quadratic flank clipping.
        PeakGeometryValidation validatePeakGeometry(const RegCoeffs &coefficients, size_t minimumSamplesPerSide = 2);
        PeakGeometryValidation validatePeakGeometry(const HyperbolicModel::Coefficients &coefficients,
            double sharedSlope, size_t scale, double step, size_t minimumSamplesPerSide = 2);

        // Fixed-window log-OLS diagnostics under independent homoskedastic
        // Gaussian errors. Supply the fitted log residual MSE and residual DF
        // of the same window. Quadratic expects 2*h-3 DF; hyperbolic expects
        // m*(2*h-2)-1 DF for m aligned traces.
        // Quadratic uses the full symmetric design at coefficients.scale.
        // These tests do not calibrate a data-selected window or a scan.
        PeakValidationStatus validatePeakFormConfidence(const RegCoeffs &coefficients,
            double logMse, size_t residualDf, double alpha = .05);
        PeakValidationStatus validatePeakFormConfidence(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, double logMse, size_t residualDf, double alpha = .05);

        // Fraction of fitted area between the first and last sample. The
        // quadratic denominator uses the legacy peak domain (positive
        // curvature flanks end at their continuous valley). Hyperbolic uses
        // the whole real line. An affine axis scale and the amplitude cancel.
        // Pass if the approximate one-sided, fixed-window delta-method lower
        // bound meets the caller's minimumFraction (0 < minimumFraction < 1,
        // 0 < alpha < .5). Same noise contract as form confidence. A saturated
        // fraction or unresolved derivative can leave the bound unavailable.
        PeakWindowCoverage validatePeakWindowCoverage(const RegCoeffs &coefficients,
            const PeakFormNoise &noise, double minimumFraction);
        PeakWindowCoverage validatePeakWindowCoverage(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, double minimumFraction);

        struct PeakAreaEvidence
        {
            PeakValidationStatus status = PeakValidationStatus::unavailable;
            bool shapeCertified = false;
            double logLowerBound = std::numeric_limits<double>::quiet_NaN();
            double left = std::numeric_limits<double>::quiet_NaN();
            double right = std::numeric_limits<double>::quiet_NaN();
        };

        // Test area of exp(E[log intensity]) against an explicit positive minimum, in
        // intensity * sample index (quadratic) or intensity * u (hyperbolic).
        // The lower bound integrates only within the observed window and the
        // certified peak domain. A valid but uncertified shape fails with a
        // zero area bound (logLowerBound = -infinity).
        // Supply ORIGINAL log-OLS coefficients, before raw amplitude refit,
        // and their same-window residual variance and DF. This fixed-window,
        // per-trace test assumes independent homoskedastic Gaussian log errors
        // and a prespecified model/window/threshold. The target is the modeled
        // curve, without arithmetic-mean bias correction or background subtraction.
        PeakAreaEvidence validatePeakAreaEvidence(const RegCoeffs &coefficients,
            const PeakFormNoise &noise, double minimumArea);
        PeakAreaEvidence validatePeakAreaEvidence(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise, double minimumArea);

        struct PeakCandidateValidation
        {
            PeakGeometryValidation geometry;
            PeakValidationStatus formConfidence = PeakValidationStatus::notEvaluated;
            PeakWindowCoverage windowCoverage{PeakValidationStatus::notEvaluated};
            PeakAreaEvidence areaEvidence{PeakValidationStatus::notEvaluated};
        };

        // Run the common ordered stages, stopping at the first non-pass.
        // The area evidence and subsequent coverage stages run only when
        // their respective minimum is explicitly supplied.
        PeakCandidateValidation validatePeakCandidate(const RegCoeffs &coefficients,
            const PeakFormNoise &noise, size_t minimumSamplesPerSide = 2,
            std::optional<double> minimumCoverage = std::nullopt, std::optional<double> minimumArea = std::nullopt);
        PeakCandidateValidation validatePeakCandidate(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise,
            size_t minimumSamplesPerSide = 2, std::optional<double> minimumCoverage = std::nullopt,
            std::optional<double> minimumArea = std::nullopt);

        struct PeakBackgroundOptions
        {
            double alpha = .05;
            size_t rotations = 999;
            uint64_t seed = 0;
            size_t minimumSamplesPerSide = 2;
        };

        struct PeakBackgroundComparison
        {
            double backgroundRSS = std::numeric_limits<double>::quiet_NaN();
            double explainedFraction = std::numeric_limits<double>::quiet_NaN();
            double pValue = std::numeric_limits<double>::quiet_NaN();
        };

        struct PeakBackgroundValidation
        {
            PeakValidationStatus status = PeakValidationStatus::unavailable;
            bool supportedPeakShape = false;
            PeakBackgroundComparison constant, linear;
        };

        // One prespecified model and equally spaced, odd-length window. The fitted
        // alternative adds a nonnegative peak amplitude to each unconstrained
        // background. Smallest reported p is 1/(rotations+1).
        // Rotation calibration assumes independent, homoskedastic Gaussian
        // errors on uncensored raw intensities. A window with any nonpositive
        // observation supplies no peak template; simulations obey the same rule.
        // fixedStep is 1 for quadratic and the model step in u for hyperbolic.
        // Pass requires rejecting both backgrounds. Adaptive model/window
        // selection, clipping, heteroskedasticity and serial correlation are
        // outside this calibration. Choose the simulation seed before testing.
        PeakBackgroundValidation validatePeakBackground(PeakModelKind model, const float *intensities,
            size_t length, double fixedStep, const PeakBackgroundOptions &options = {});
    }
}

#endif
