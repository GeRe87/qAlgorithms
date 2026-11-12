#include "qalgorithms_peak_search.h"
#include "detail/qalgorithms_peak_candidates.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace qAlgorithms
{
    namespace peak_model
    {
        namespace
        {
            bool windowLess(const PeakSearchWindow &left, const PeakSearchWindow &right)
            {
                return left.center < right.center || (left.center == right.center && left.scale < right.scale);
            }

            std::vector<UnresolvedSearchGroup> collectUnresolvedGroups(
                const std::vector<PeakSearchResult> &windows, const std::vector<size_t> &indices)
            {
                std::vector<UnresolvedSearchGroup> groups;
                std::vector<size_t> order = indices;
                std::sort(order.begin(), order.end(), [&](size_t left, size_t right)
                {
                    const auto &a = windows[left].window;
                    const auto &b = windows[right].window;
                    const size_t startA = a.center - a.scale, startB = b.center - b.scale;
                    return startA < startB || (startA == startB && left < right);
                });
                UnresolvedSearchGroup group;
                auto appendGroup = [&]()
                {
                    if (group.windowIndices.size() > 1)
                    {
                        std::sort(group.windowIndices.begin(), group.windowIndices.end());
                        groups.push_back(std::move(group));
                    }
                    group = {};
                };
                for (size_t index : order)
                {
                    const auto &window = windows[index].window;
                    const size_t start = window.center - window.scale, end = window.center + window.scale;
                    if (!group.windowIndices.empty() && start > group.end)
                    {
                        appendGroup();
                    }
                    if (group.windowIndices.empty())
                    {
                        group.start = start;
                        group.end = end;
                    }
                    else
                    {
                        group.end = std::max(group.end, end);
                    }
                    group.windowIndices.push_back(index);
                }
                appendGroup();
                std::sort(groups.begin(), groups.end(),
                    [](const UnresolvedSearchGroup &left, const UnresolvedSearchGroup &right)
                {
                    return left.windowIndices.front() < right.windowIndices.front();
                });
                return groups;
            }

            DiagnosticApexGroup classifyApices(const float *raw, const UnresolvedSearchGroup &group,
                const HyperbolicSearchResult &search)
            {
                DiagnosticApexGroup result;
                for (size_t first = group.start; first <= group.end;)
                {
                    size_t last = first;
                    while (last < group.end && raw[last + 1] == raw[first])
                    {
                        ++last;
                    }
                    if (first > group.start && last < group.end
                        && raw[first - 1] < raw[first] && raw[last + 1] < raw[last])
                    {
                        result.maxima.push_back({first, last, {}});
                    }
                    first = last + 1;
                }
                for (size_t index : group.windowIndices)
                {
                    const auto &entry = search.windows[index];
                    const auto &window = entry.window;
                    const double apex = static_cast<double>(window.center)
                        + entry.inference.properties.front().apexU.value / search.fixedStep;
                    const auto left = static_cast<size_t>(std::floor(apex));
                    const auto right = static_cast<size_t>(std::ceil(apex));
                    // As in the quadratic selector, an observed maximum must
                    // touch a sample bracketing the fitted apex.
                    const auto maximum = std::lower_bound(result.maxima.begin(), result.maxima.end(), left,
                        [](const ObservedApexHypothesis &candidate, size_t sample) { return candidate.last < sample; });
                    if (maximum != result.maxima.end() && maximum->first <= right
                        && maximum->first >= window.center - window.scale && maximum->last <= window.center + window.scale)
                    {
                        maximum->windowIndices.push_back(index);
                    }
                    else
                    {
                        result.unassignedWindowIndices.push_back(index);
                    }
                }
                return result;
            }
        }

        HyperbolicSearchPlan::HyperbolicSearchPlan(
            size_t length,
            double fixedStep,
            std::vector<PeakSearchWindow> windows,
            const HyperbolicSearchOptions &options)
            : length_(length), fixedStep_(fixedStep), windows_(std::move(windows)), options_(options), familySize_(0)
        {
            if (length < 5 || windows_.empty() || !(options.familyAlpha > 0 && options.familyAlpha < 1)
                || !std::isfinite(options.familyAlpha)
                || (options.selection != WindowSelection::exploratory && options.selection != WindowSelection::prespecified))
            {
                throw std::invalid_argument("Invalid declared search plan");
            }
            size_t maxScale = 2;
            for (const auto &window : windows_)
            {
                if (window.scale < 2 || window.center >= length || window.scale > window.center
                    || window.scale > length - 1 - window.center)
                {
                    throw std::invalid_argument("Declared search window outside trace");
                }
                maxScale = std::max(maxScale, window.scale);
            }
            // Validate all geometry independently of observations. The existing
            // selected-window API owns its own caches during evaluation.
            const HyperbolicModel::Geometry geometry(maxScale, fixedStep);
            (void)geometry;
            std::sort(windows_.begin(), windows_.end(), windowLess);
            for (size_t i = 1; i < windows_.size(); ++i)
            {
                if (windows_[i - 1].center == windows_[i].center && windows_[i - 1].scale == windows_[i].scale)
                {
                    throw std::invalid_argument("Duplicate window in declared search family");
                }
            }
            if (options.totalFamilySize != 0 && options.totalFamilySize < windows_.size())
            {
                throw std::invalid_argument("Declared family cannot be smaller than complete search plan");
            }
            familySize_ = options.totalFamilySize == 0 ? windows_.size() : options.totalFamilySize;
        }

        const std::vector<PeakSearchWindow> &HyperbolicSearchPlan::windows() const
        {
            return windows_;
        }

        size_t HyperbolicSearchPlan::familySize() const
        {
            return familySize_;
        }

        size_t HyperbolicSearchPlan::length() const
        {
            return length_;
        }

        HyperbolicSearchResult HyperbolicSearchPlan::evaluate(const std::vector<const float *> &rawTraces) const
        {
            if (rawTraces.empty())
            {
                throw std::invalid_argument("Search requires at least one aligned raw trace");
            }
            for (const float *raw : rawTraces)
            {
                if (raw == nullptr)
                {
                    throw std::invalid_argument("Null raw search trace");
                }
                for (size_t i = 0; i < length_; ++i)
                {
                    if (!(raw[i] > 0 && std::isfinite(raw[i])))
                    {
                        throw std::invalid_argument("Search requires positive finite uncensored observations");
                    }
                }
            }
            HyperbolicInferenceOptions inferenceOptions;
            inferenceOptions.selection = options_.selection;
            inferenceOptions.familyAlpha = options_.familyAlpha;
            inferenceOptions.minimumSamplesPerSide = options_.minimumSamplesPerSide;
            inferenceOptions.familySize = familySize_;
            HyperbolicSearchResult result;
            result.length = length_;
            result.traceCount = rawTraces.size();
            result.fixedStep = fixedStep_;
            result.options = options_;
            result.familySize = familySize_;
            result.windows.reserve(windows_.size());
            detail::findValidRegressions<PeakSearchResult>(windows_,
                [&](const PeakSearchWindow &window, PeakSearchResult *entry)
            {
                entry->window = window;
                entry->inference = inferHyperbolicWindow(rawTraces, length_, window.center, window.scale, fixedStep_, inferenceOptions);
                return entry->inference.groupDecision == ConditionalDecision::passed;
            }, [&](PeakSearchResult &entry, bool valid)
            {
                if (valid)
                {
                    result.candidateIndices.push_back(result.windows.size());
                }
                result.windows.push_back(std::move(entry));
            });
            result.unresolvedGroups = collectUnresolvedGroups(result.windows, result.candidateIndices);
            for (size_t i = 0; i < result.unresolvedGroups.size(); ++i)
            {
                for (size_t index : result.unresolvedGroups[i].windowIndices)
                {
                    result.windows[index].unresolvedGroup = i;
                }
            }
            return result;
        }

        SegmentedHyperbolicSearchResult HyperbolicSearchPlan::evaluateSegmented(
            const std::vector<AlignedIntensityTrace> &traces, double coordinateStep,
            double relativeTolerance, double threshold) const
        {
            for (const auto &trace : traces)
            {
                if (trace.length != length_)
                {
                    throw std::invalid_argument("Trace length does not match declared search plan");
                }
            }
            SegmentedHyperbolicSearchResult result;
            result.segmentation = segmentIntensities(traces, coordinateStep, relativeTolerance, threshold, MINLENGTH);
            auto &search = result.search;
            search.length = length_;
            search.traceCount = traces.size();
            search.fixedStep = fixedStep_;
            search.options = options_;
            search.options.selection = WindowSelection::exploratory;
            search.familySize = familySize_;
            search.windows.reserve(windows_.size());
            result.windowStatus.reserve(windows_.size());

            const HyperbolicInferenceOptions inferenceOptions{
                WindowSelection::exploratory, familySize_, options_.familyAlpha, options_.minimumSamplesPerSide};
            std::vector<const float *> raw;
            for (const auto &trace : traces)
            {
                raw.push_back(trace.intensity);
            }
            size_t maxScale = 2;
            for (const auto &window : windows_)
            {
                maxScale = std::max(maxScale, window.scale);
            }
            const HyperbolicModel::Geometry cache(maxScale, fixedStep_);
            auto span = result.segmentation.spans.begin();
            for (const auto &window : windows_)
            {
                while (span != result.segmentation.spans.end() && span->end <= window.center)
                {
                    ++span;
                }
                PeakSearchResult entry;
                entry.window = window;
                auto status = SegmentedWindowStatus::outsideEligibleSpan;
                if (span != result.segmentation.spans.end() && span->status == IntensitySpanStatus::eligible
                    && window.center - window.scale >= span->start && window.center + window.scale < span->end)
                {
                    entry.inference = inferHyperbolicWindow(raw, length_, window.center, window.scale, cache, inferenceOptions);
                    status = SegmentedWindowStatus::evaluated;
                    if (entry.inference.status == InferenceStatus::available && entry.inference.observedSupportSufficient)
                    {
                        result.diagnosticIndices.push_back(search.windows.size());
                    }
                }
                search.windows.push_back(std::move(entry));
                result.windowStatus.push_back(status);
            }
            result.diagnosticGroups = collectUnresolvedGroups(search.windows, result.diagnosticIndices);
            if (traces.size() == 1)
            {
                result.apexGroups.reserve(result.diagnosticGroups.size());
                for (const auto &group : result.diagnosticGroups)
                {
                    result.apexGroups.push_back(classifyApices(raw.front(), group, search));
                }
            }
            return result;
        }

        HyperbolicSearchPlan makeUniformHyperbolicSearchPlan(
            size_t length,
            double fixedStep,
            size_t minScale,
            size_t maxScale,
            size_t centerStride,
            const HyperbolicSearchOptions &options)
        {
            if (length < 5 || minScale < 2 || minScale > maxScale || maxScale > (length - 1) / 2 || centerStride == 0)
            {
                throw std::invalid_argument("Invalid uniform search geometry");
            }
            std::vector<PeakSearchWindow> windows;
            const size_t lastCenter = length - 1 - minScale;
            for (size_t center = minScale; center <= lastCenter;)
            {
                const size_t largest = std::min(maxScale, std::min(center, length - 1 - center));
                for (size_t scale = minScale; scale <= largest; ++scale)
                {
                    windows.push_back({center, scale});
                }
                if (centerStride > lastCenter - center)
                {
                    break;
                }
                center += centerStride;
            }
            return HyperbolicSearchPlan(length, fixedStep, std::move(windows), options);
        }

        ModelPeakSearchResult findPeaks(const PeakSearchRequest &request)
        {
            switch (request.model)
            {
            case PeakModelKind::PiecewiseQuadratic:
            {
                PiecewiseQuadraticSearchResult result;
                result.status = qpeaks_find(request.intensity, request.coordinates, request.degreesOfFreedom,
                    request.length, &result.peaks, request.describeUnresolved ? &result.unresolvedGroups : nullptr);
                return result;
            }
            case PeakModelKind::Hyperbolic:
                if (request.hyperbolicPlan == nullptr || request.hyperbolicPlan->length() != request.length)
                {
                    throw std::invalid_argument("Hyperbolic search requires a plan matching the trace length");
                }
                if (request.segmentation != nullptr)
                {
                    const auto &segmentation = *request.segmentation;
                    return request.hyperbolicPlan->evaluateSegmented(
                        {{request.intensity, segmentation.coordinates, request.length}},
                        segmentation.coordinateStep, segmentation.relativeTolerance, segmentation.threshold);
                }
                return request.hyperbolicPlan->evaluate({request.intensity});
            }
            throw std::invalid_argument("Unknown peak model");
        }
    }
}
