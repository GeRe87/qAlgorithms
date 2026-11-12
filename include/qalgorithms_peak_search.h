#ifndef QALGORITHMS_PEAK_SEARCH_H
#define QALGORITHMS_PEAK_SEARCH_H

#include "qalgorithms_peak_inference.h"
#include "qalgorithms_intensity_segments.h"

#include <variant>

namespace qAlgorithms
{
    namespace peak_model
    {
        struct PeakSearchWindow
        {
            size_t center = 0;
            size_t scale = 0;
        };

        struct HyperbolicSearchOptions
        {
            WindowSelection selection = WindowSelection::exploratory;
            double familyAlpha = 0.05;
            size_t minimumSamplesPerSide = 2;
            // Zero uses the entire plan size. A larger declared family also
            // accounts for other trace pools/steps/plans tested by the caller.
            // A nonzero value smaller than the plan size is rejected.
            size_t totalFamilySize = 0;
        };

        struct PeakSearchResult
        {
            PeakSearchWindow window;
            HyperbolicInference inference;
            // Indices refer to the complete canonical window list, not a
            // filtered list of fits. No winning hypothesis is selected.
            size_t unresolvedGroup = std::numeric_limits<size_t>::max();
        };

        struct UnresolvedSearchGroup
        {
            size_t start = 0, end = 0; // inclusive observed sample envelope
            std::vector<size_t> windowIndices;
        };

        struct HyperbolicSearchResult
        {
            size_t length = 0, traceCount = 0;
            double fixedStep = 0;
            HyperbolicSearchOptions options;
            size_t familySize = 0;
            std::vector<PeakSearchResult> windows;
            std::vector<size_t> candidateIndices;
            // Connected components of overlapping candidate support intervals.
            // Multiple scales/windows may describe the SAME peak. A component
            // establishes neither interference nor a physical analyte count.
            std::vector<UnresolvedSearchGroup> unresolvedGroups;
        };

        enum class SegmentedWindowStatus { outsideEligibleSpan, evaluated };

        struct ObservedApexHypothesis
        {
            size_t first = 0, last = 0; // inclusive observed maximum plateau
            std::vector<size_t> windowIndices;
        };

        struct DiagnosticApexGroup
        {
            std::vector<ObservedApexHypothesis> maxima;
            std::vector<size_t> unassignedWindowIndices;
        };

        struct SegmentedHyperbolicSearchResult
        {
            HyperbolicSearchResult search;
            IntensitySegments segmentation;
            std::vector<SegmentedWindowStatus> windowStatus;
            // Exploratory fits with sufficient support and an integrable shape.
            // Indices address search.windows; groups are inclusive sample envelopes.
            std::vector<size_t> diagnosticIndices;
            std::vector<UnresolvedSearchGroup> diagnosticGroups;
            // Single-trace only. One entry per diagnosticGroups element, with
            // unsupported observed maxima and unassigned fits retained.
            std::vector<DiagnosticApexGroup> apexGroups;
        };

        // A data-independent, fixed family of windows, canonicalized by center
        // then scale. Configure before observing raw data; adaptive construction
        // cannot acquire calibration merely by declaring 'prespecified'.
        class HyperbolicSearchPlan
        {
        public:
            HyperbolicSearchPlan(size_t length,
                                 double fixedStep,
                                 std::vector<PeakSearchWindow> windows,
                                 const HyperbolicSearchOptions &options = {});

            const std::vector<PeakSearchWindow> &windows() const;
            size_t length() const;
            size_t familySize() const;

            // Validate ALL traces before fitting. Every declared window remains
            // in the result, including invalid/nonintegrable/numeric failures;
            // none are removed from the multiple-testing denominator.
            // Covariance is the fixed-window reference only: selecting a result
            // does not give post-selection confidence intervals or area coverage.
            HyperbolicSearchResult evaluate(const std::vector<const float *> &rawTraces) const;

            // Original indices address the supplied coordinate arrays. Threshold
            // selection yields exploratory inference over the complete family.
            SegmentedHyperbolicSearchResult evaluateSegmented(const std::vector<AlignedIntensityTrace> &traces,
                double coordinateStep, double relativeTolerance = 0,
                double threshold = 2.7182818284590452) const;

        private:
            size_t length_;
            double fixedStep_;
            std::vector<PeakSearchWindow> windows_;
            HyperbolicSearchOptions options_;
            size_t familySize_;
        };

        // Enumerate a uniform, data-independent family. Center grid starts at
        // minScale; each center includes every fitting scale in [min,max].
        HyperbolicSearchPlan makeUniformHyperbolicSearchPlan(
            size_t length,
            double fixedStep,
            size_t minScale,
            size_t maxScale,
            size_t centerStride = 1,
            const HyperbolicSearchOptions &options = {});

        struct PiecewiseQuadraticSearchResult
        {
            int32_t status = 0; // qpeaks_find return value
            std::vector<RegressionGauss> peaks;
            std::vector<UnresolvedPeak> unresolvedGroups;
        };

        struct PeakSearchRequest
        {
            PeakModelKind model = PeakModelKind::PiecewiseQuadratic;
            const float *intensity = nullptr;
            size_t length = 0;
            const float *coordinates = nullptr; // piecewise-quadratic only
            const uint16_t *degreesOfFreedom = nullptr; // cumulative, piecewise-quadratic only
            bool describeUnresolved = false; // piecewise-quadratic only
            const HyperbolicSearchPlan *hyperbolicPlan = nullptr;
            struct Segmentation
            {
                const double *coordinates = nullptr;
                double coordinateStep = 0;
                double relativeTolerance = 0;
                double threshold = 2.7182818284590452;
            };
            const Segmentation *segmentation = nullptr; // optional, Hyperbolic only
        };

        // Single-trace search. The quadratic branch uses the legacy selection.
        // Hyperbolic requires a matching plan. Without segmentation it requires
        // positive finite intensities. With segmentation the result is exploratory.
        // Windows use original sample indices; fit quantities use the plan's u
        // coordinates. An apex maps to center + apexU / fixedStep in sample space.
        using ModelPeakSearchResult = std::variant<PiecewiseQuadraticSearchResult, HyperbolicSearchResult,
            SegmentedHyperbolicSearchResult>;
        ModelPeakSearchResult findPeaks(const PeakSearchRequest &request);
    }
}

#endif
