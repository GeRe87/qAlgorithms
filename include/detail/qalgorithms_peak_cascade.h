#ifndef QALGORITHMS_DETAIL_PEAK_CASCADE_H
#define QALGORITHMS_DETAIL_PEAK_CASCADE_H

#include "qalgorithms_peak_validation.h"
#include <cmath>

namespace qAlgorithms
{
    namespace peak_model
    {
        inline PeakGeometryValidation checkPeakGeometry(bool shape, double apex, size_t scale,
            double step, size_t minimumSamplesPerSide)
        {
            PeakGeometryValidation result;
            result.form = shape ? PeakValidationStatus::passed : PeakValidationStatus::failed;
            if (!shape)
            {
                return result;
            }
            result.apex = apex;
            result.apexSupport = PeakValidationStatus::unavailable;
            if (scale == 0 || !(step > 0 && std::isfinite(step)))
            {
                return result;
            }
            const double edge = static_cast<double>(scale) * step;
            if (!std::isfinite(edge))
            {
                return result;
            }
            const double side = minimumSamplesPerSide == 0 ? 0
                : (static_cast<double>(scale) - static_cast<double>(minimumSamplesPerSide - 1)) * step;
            const bool supported = std::isfinite(apex) && (minimumSamplesPerSide == 0
                ? std::abs(apex) <= edge : std::abs(apex) < side);
            result.apexSupport = supported ? PeakValidationStatus::passed : PeakValidationStatus::failed;
            return result;
        }

        // Model adapters supply the stage calculations. This routine owns
        // their order and leaves skipped stages marked notEvaluated.
        template<typename CheckConfidence, typename CheckArea, typename CheckCoverage>
        PeakCandidateValidation runPeakCascade(PeakGeometryValidation geometry, CheckConfidence checkConfidence,
            std::optional<double> minimumArea, CheckArea checkArea,
            std::optional<double> minimumCoverage, CheckCoverage checkCoverage)
        {
            PeakCandidateValidation result{geometry};
            if (geometry.apexSupport == PeakValidationStatus::passed)
            {
                result.formConfidence = checkConfidence();
            }
            if (result.formConfidence != PeakValidationStatus::passed)
            {
                return result;
            }
            if (minimumArea)
            {
                result.areaEvidence = checkArea(*minimumArea);
                if (result.areaEvidence.status != PeakValidationStatus::passed)
                {
                    return result;
                }
            }
            if (minimumCoverage)
            {
                result.windowCoverage = checkCoverage(*minimumCoverage);
            }
            return result;
        }

        inline PeakValidationStatus peakCascadeStatus(const PeakCandidateValidation &stages)
        {
            for (PeakValidationStatus stage : {stages.geometry.form, stages.geometry.apexSupport,
                stages.formConfidence, stages.areaEvidence.status, stages.windowCoverage.status})
            {
                if (stage != PeakValidationStatus::passed)
                {
                    return stage;
                }
            }
            return PeakValidationStatus::passed;
        }
    }
}

#endif
