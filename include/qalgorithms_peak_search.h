#ifndef QALGORITHMS_PEAK_SEARCH_H
#define QALGORITHMS_PEAK_SEARCH_H

#include "qalgorithms_peak_inference.h"
#include "qalgorithms_intensity_segments.h"
#include "qalgorithms_peak_validation.h"

#include <variant>
#include <optional>
#include <memory>
#include <string>

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
        // supportedWindows skips inference for shapes with no finite,
        // adequately supported apex. Retained windows get full covariance.
        enum class SegmentedInferenceDetail { allWindows, supportedWindows };

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

        enum class ApexSelectionStatus { selected, unsupported, insufficientCommonSupport, ambiguous, numericalFailure };

        struct ApexSelection
        {
            size_t first = 0, last = 0; // observed maximum plateau
            ApexSelectionStatus status = ApexSelectionStatus::unsupported;
            size_t windowIndex = std::numeric_limits<size_t>::max(); // valid when selected
        };
        using HyperbolicApexSelection = ApexSelection;

        struct AffineSampleAxis
        {
            double origin = 0; // measurement coordinate of sample index zero
            double unitsPerSample = 1; // signed and nonzero
            AxisUnit unit = AxisUnit::CallerDefined;
        };

        struct MeasuredHyperbolicPeak
        {
            HyperbolicApexSelection selection;
            double position, height, area;
            // Row-major order: position, height, area. Fixed-window
            // covariance, conditional on the supplied axis.
            std::array<double, 9> covariance;
        };

        struct MeasuredHyperbolicPeaks
        {
            AxisUnit axisUnit = AxisUnit::CallerDefined;
            std::vector<MeasuredHyperbolicPeak> peaks;
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
            // Opt-in, single-trace, ordered by observed maximum. Selection
            // minimizes raw-intensity SSE on the intersection of the fit
            // windows, with at least MINLENGTH samples and a sample on each
            // side of the plateau. Exact ties between different estimates
            // are ambiguous. Covariance in search.windows is fixed-window.
            std::vector<HyperbolicApexSelection> apexSelections;
            // Singleton fits without an assignable observed maximum. For
            // overlapping windows the same information lives in apexGroups.
            std::vector<size_t> unassignedSingletons;
            std::optional<MeasuredHyperbolicPeaks> measurementPeaks;
        };

        // Map selected representatives on an externally calibrated affine
        // sample axis. The scale must apply to the entire modeled peak;
        // area has units intensity times measurement coordinate.
        // Unresolved hypotheses stay in apexSelections. Each
        // covariance belongs to one selected fit; no cross-peak covariance,
        // calibration error, or selection uncertainty is inferred.
        MeasuredHyperbolicPeaks mapSelectedHyperbolicPeaks(
            const SegmentedHyperbolicSearchResult &result, const AffineSampleAxis &axis);

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
                double threshold = 2.7182818284590452,
                SegmentedInferenceDetail detail = SegmentedInferenceDetail::allWindows,
                bool selectRepresentatives = false) const;

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

        std::vector<PeakSearchWindow> makeUniformPeakSearchWindows(
            size_t length, size_t minScale, size_t maxScale, size_t centerStride = 1);

        class NativePeakModel;

        struct CandidateModelWindows
        {
            PeakModelKind model = PeakModelKind::PiecewiseQuadratic;
            double fixedStep = 1;
            FullPeakCandidateOptions options;
            std::vector<PeakSearchWindow> windows;
        };

        struct NativeCandidateModelWindows
        {
            // The plan retains its own reference to the loaded library.
            std::shared_ptr<const NativePeakModel> model;
            double fixedStep = 1;
            FullPeakCandidateOptions options;
            std::vector<PeakSearchWindow> windows;
        };

        struct CandidateSearchWindowResult
        {
            // Native windows carry their descriptor name and no built-in enum.
            std::optional<PeakModelKind> model;
            PeakSearchWindow window;
            FullPeakCandidateValidation validation;
            double adjustedBackgroundP = std::numeric_limits<double>::quiet_NaN();
            std::string nativeModelName;
        };

        struct CandidateSearchResult
        {
            size_t familySize = 0;
            double backgroundAlphaPerWindow = 0;
            std::vector<CandidateSearchWindowResult> windows;
            std::vector<size_t> candidateIndices;
            struct MeasuredProperties
            {
                double position, area;
                std::array<double, 4> covariance; // row-major position, area
            };
            struct Measurement
            {
                ApexSelection selection;
                // Missing for unresolved hypotheses and for selected fits
                // without representable joint fixed-window covariance.
                std::optional<MeasuredProperties> properties;
            };
            struct Measurements
            {
                AxisUnit axisUnit = AxisUnit::CallerDefined;
                std::vector<Measurement> peaks;
                std::vector<size_t> unassignedWindowIndices;
            };
            std::optional<Measurements> measurements;
        };

        // Prespecify all windows and models before seeing intensities. A window
        // tested with two models counts twice. Each built-in kind and each
        // native descriptor name must be unique within the plan.
        // Bonferroni applies to the raw-background stage across the full family,
        // even when windows overlap. Earlier gates can only remove rejections.
        // The guarantee assumes the fixed-window raw Gaussian null for each
        // true null and independent random rotations within each test. Choose
        // the seeds before inspecting the data. It does not calibrate log-stage
        // claims, adaptive plans or multiple calls treated as a single family.
        class CandidateSearchPlan
        {
        public:
            // Uses familyAlpha in place of each model's background.alpha.
            // The rotation budget per model must resolve familyAlpha/familySize.
            CandidateSearchPlan(size_t length, std::vector<CandidateModelWindows> models, double familyAlpha = .05,
                std::vector<NativeCandidateModelWindows> nativeModels = {});

            // A supplied axis requests one representative per observed apex,
            // across the accepted windows and models. The affine calibration
            // must cover the full modeled peak. Covariance is conditional
            // on the selected window, model and axis.
            CandidateSearchResult evaluate(const float *intensities, size_t length,
                const AffineSampleAxis *measurementAxis = nullptr) const;

        private:
            size_t length_, familySize_ = 0;
            double windowAlpha_ = 0;
            struct PlannedModel
            {
                std::optional<PeakModelKind> kind;
                std::shared_ptr<const NativePeakModel> module;
                std::string name;
                double fixedStep;
                FullPeakCandidateOptions options;
                std::vector<PeakSearchWindow> windows;
            };
            std::vector<PlannedModel> models_;
        };

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
                SegmentedInferenceDetail detail = SegmentedInferenceDetail::allWindows;
                bool selectRepresentatives = false;
                // Requests representative selection and mapped output.
                // Independent of coordinates used for segmentation.
                const AffineSampleAxis *measurementAxis = nullptr;
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
        CandidateSearchResult findPeaks(const CandidateSearchPlan &plan, const float *intensities, size_t length,
            const AffineSampleAxis *measurementAxis = nullptr);
    }
}

#endif
