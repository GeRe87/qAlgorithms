#ifndef QALGORITHMS_PEAK_VALIDATION_H
#define QALGORITHMS_PEAK_VALIDATION_H

#include "qalgorithms_peak_fit.h"

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

        struct PeakCandidateValidation
        {
            PeakGeometryValidation geometry;
            PeakValidationStatus formConfidence = PeakValidationStatus::notEvaluated;
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

        // Run the common ordered stages, stopping at the first non-pass.
        PeakCandidateValidation validatePeakCandidate(const RegCoeffs &coefficients,
            const PeakFormNoise &noise, size_t minimumSamplesPerSide = 2);
        PeakCandidateValidation validatePeakCandidate(const HyperbolicModel::Fit &fit, size_t trace,
            const HyperbolicModel::Geometry &geometry, const PeakFormNoise &noise,
            size_t minimumSamplesPerSide = 2);
    }
}

#endif
