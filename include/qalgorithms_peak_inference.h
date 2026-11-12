#ifndef QALGORITHMS_PEAK_INFERENCE_H
#define QALGORITHMS_PEAK_INFERENCE_H

#include "qalgorithms_peak_fit.h"

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace qAlgorithms
{
    namespace peak_model
    {
        enum class WindowSelection
        {
            exploratory,
            prespecified
        };

        enum class ConditionalDecision
        {
            unavailable,
            rejected,
            passed
        };

        enum class InferenceStatus
        {
            unavailable,
            available,
            invalidFit,
            degenerateNoise,
            numericalFailure,
            screenedOut // exploratory geometry/support screen; properties were not evaluated
        };

        // The caller must explicitly declare fixed windows and the complete
        // tested family to enable a conditional gate. This does not calibrate
        // data-selected windows or correlated/non-Gaussian observations.
        struct HyperbolicInferenceOptions
        {
            WindowSelection selection = WindowSelection::exploratory;
            size_t familySize = 1;
            double familyAlpha = 0.05;
            size_t minimumSamplesPerSide = 2;
        };

        struct HyperbolicInference
        {
            InferenceStatus status = InferenceStatus::unavailable;
            ConditionalDecision groupDecision = ConditionalDecision::unavailable;
            bool observedSupportSufficient = false; // meaningful when status=available
            size_t residualDegreesOfFreedom = 0;
            size_t comparisonDegreesOfFreedom = 0;
            double residualSumSquares = 0;
            double nullResidualSumSquares = 0;
            double logNoiseVariance = std::numeric_limits<double>::quiet_NaN();
            double comparisonF = std::numeric_limits<double>::quiet_NaN();
            double nominalP = std::numeric_limits<double>::quiet_NaN();
            double adjustedP = std::numeric_limits<double>::quiet_NaN();
            double logNominalP = std::numeric_limits<double>::quiet_NaN();
            double logAdjustedP = std::numeric_limits<double>::quiet_NaN();
            std::vector<HyperbolicProperties> properties;
            // Natural-property covariance, row-major order:
            // [apexU_0, height_0, areaU_0, apexU_1, height_1, areaU_1, ...].
            // Includes cross-trace shared-shape and amplitude-refit dependence.
            // First-order delta method, conditional on fixed grid/window/model;
            // excludes selection, calibration and model-misspecification error.
            std::vector<double> covariance;
        };

        // Fit and infer one EXPLICITLY selected window from the same positive
        // raw observations (log samples are derived here, never independently
        // supplied). No scan-wide detection, ambiguity resolution or CLI change.
        // Assumes independent homoskedastic Gaussian errors in log intensity.
        HyperbolicInference inferHyperbolicWindow(
            const std::vector<const float *> &rawTraces,
            size_t length,
            size_t center,
            size_t scale,
            double fixedStep,
            const HyperbolicInferenceOptions &options = {});

        // Reuse geometry across windows with the same fixed step.
        HyperbolicInference inferHyperbolicWindow(
            const std::vector<const float *> &rawTraces, size_t length, size_t center, size_t scale,
            const HyperbolicModel::Geometry &cache, const HyperbolicInferenceOptions &options = {});

        struct AffinePeakCoordinates
        {
            double origin = 0;       // physical coordinate of window center
            double unitsPerU = 1;    // fixed, nonzero; negative axes permitted
            AxisUnit unit = AxisUnit::CallerDefined;
        };

        struct MappedHyperbolicProperties
        {
            AxisUnit axisUnit = AxisUnit::CallerDefined;
            // Same order as covariance: position, height, area per trace.
            std::vector<std::array<double, 3>> values;
            std::vector<double> covariance;
        };

        // Only an externally supplied affine map is supported. The area uses
        // abs(unitsPerU), position/cross-covariances retain its signed derivative.
        // No inferred instrument calibration or nonlinear Jacobian approximation.
        MappedHyperbolicProperties mapHyperbolicAffine(
            const HyperbolicInference &inference,
            const AffinePeakCoordinates &coordinates);
    }
}

#endif
