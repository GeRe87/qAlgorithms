#include "qalgorithms_peak_search.h"
#include "detail/qalgorithms_peak_candidates.h"
#include "detail/qalgorithms_peak_inference.h"

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
                const std::vector<PeakSearchResult> &windows, const std::vector<size_t> &indices,
                bool includeSingletons = false)
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
                    if (group.windowIndices.size() > 1 || (includeSingletons && !group.windowIndices.empty()))
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

            struct SelectionFit
            {
                size_t index;
                HyperbolicModel::Coefficients coefficients;
                double slope;
            };

            HyperbolicApexSelection selectRepresentative(const float *raw, const ObservedApexHypothesis &maximum,
                const HyperbolicSearchResult &search, const std::vector<SelectionFit> &fits)
            {
                HyperbolicApexSelection selection{maximum.first, maximum.last};
                const auto &indices = maximum.windowIndices;
                if (indices.empty())
                {
                    return selection;
                }
                size_t start = 0, end = search.length - 1;
                for (size_t index : indices)
                {
                    const auto &window = search.windows[index].window;
                    start = std::max(start, window.center - window.scale);
                    end = std::min(end, window.center + window.scale);
                }
                selection.status = ApexSelectionStatus::insufficientCommonSupport;
                if (start >= maximum.first || end <= maximum.last || end - start + 1 < MINLENGTH)
                {
                    return selection;
                }
                if (indices.size() == 1)
                {
                    selection.status = ApexSelectionStatus::selected;
                    selection.windowIndex = indices.front();
                    return selection;
                }
                double best = std::numeric_limits<double>::infinity();
                size_t bestIndex = indices.front();
                bool ambiguous = false;
                for (size_t index : indices)
                {
                    const auto fit = std::lower_bound(fits.begin(), fits.end(), index,
                        [](const SelectionFit &candidate, size_t i) { return candidate.index < i; });
                    const auto &window = search.windows[index].window;
                    double sse = 0;
                    for (size_t i = start; i <= end; ++i)
                    {
                        const double u = search.fixedStep * (static_cast<double>(i) - static_cast<double>(window.center));
                        const double error = raw[i] - std::exp(detail::hyperbolic::predict(fit->coefficients, fit->slope, u));
                        sse += error * error;
                    }
                    if (!std::isfinite(sse))
                    {
                        selection.status = ApexSelectionStatus::numericalFailure;
                        return selection;
                    }
                    if (sse < best)
                    {
                        best = sse;
                        bestIndex = index;
                        ambiguous = false;
                    }
                    else if (sse == best)
                    {
                        const auto &a = search.windows[bestIndex].inference;
                        const auto &b = search.windows[index].inference;
                        const auto &p = a.properties.front();
                        const auto &q = b.properties.front();
                        const double apexA = static_cast<double>(search.windows[bestIndex].window.center) + p.apexU.value / search.fixedStep;
                        const double apexB = static_cast<double>(window.center) + q.apexU.value / search.fixedStep;
                        ambiguous = ambiguous || apexA != apexB || p.height.value != q.height.value
                            || p.areaU.value != q.areaU.value || a.covariance != b.covariance;
                    }
                }
                selection.status = ambiguous ? ApexSelectionStatus::ambiguous : ApexSelectionStatus::selected;
                if (!ambiguous)
                {
                    selection.windowIndex = bestIndex;
                }
                return selection;
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
            double relativeTolerance, double threshold, SegmentedInferenceDetail inferenceDetail,
            bool selectRepresentatives) const
        {
            if (selectRepresentatives && traces.size() != 1)
            {
                throw std::invalid_argument("Representative selection requires a single trace");
            }
            if (inferenceDetail != SegmentedInferenceDetail::allWindows
                && inferenceDetail != SegmentedInferenceDetail::supportedWindows)
            {
                throw std::invalid_argument("Unknown segmented inference detail");
            }
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
            search.windows.resize(windows_.size());
            result.windowStatus.assign(windows_.size(), SegmentedWindowStatus::outsideEligibleSpan);
            size_t maxScale = 2;
            for (size_t i = 0; i < windows_.size(); ++i)
            {
                search.windows[i].window = windows_[i];
                maxScale = std::max(maxScale, windows_[i].scale);
            }

            const HyperbolicInferenceOptions inferenceOptions{
                WindowSelection::exploratory, familySize_, options_.familyAlpha, options_.minimumSamplesPerSide};
            const HyperbolicModel::Geometry cache(maxScale, fixedStep_);
            std::vector<std::vector<float>> logs(traces.size());
            std::vector<const float *> raw(traces.size()), logPointers(traces.size());
            std::vector<SelectionFit> selectionFits;
            auto requested = windows_.begin();
            for (const auto &span : result.segmentation.spans)
            {
                if (span.status != IntensitySpanStatus::eligible)
                {
                    continue;
                }
                while (requested != windows_.end() && requested->center < span.start)
                {
                    ++requested;
                }
                if (requested == windows_.end())
                {
                    break;
                }
                if (requested->center >= span.end)
                {
                    continue;
                }
                const size_t length = span.end - span.start;
                for (size_t j = 0; j < traces.size(); ++j)
                {
                    raw[j] = traces[j].intensity + span.start;
                    logs[j].resize(length);
                    for (size_t i = 0; i < length; ++i)
                    {
                        logs[j][i] = static_cast<float>(std::log(static_cast<double>(raw[j][i])));
                    }
                    logPointers[j] = logs[j].data();
                }
                HyperbolicModel::runningFit(logPointers, length, cache, [&](const auto &fit)
                {
                    const PeakSearchWindow window{span.start + fit.center, fit.scale};
                    while (requested != windows_.end() && windowLess(*requested, window))
                    {
                        ++requested;
                    }
                    if (requested == windows_.end() || windowLess(window, *requested))
                    {
                        return;
                    }
                    const size_t index = static_cast<size_t>(requested - windows_.begin());
                    ++requested;
                    auto &entry = search.windows[index];
                    entry.inference = inferFittedHyperbolicWindow(raw, logPointers, length, fit, cache,
                        inferenceOptions, inferenceDetail == SegmentedInferenceDetail::supportedWindows);
                    result.windowStatus[index] = SegmentedWindowStatus::evaluated;
                    if (entry.inference.status == InferenceStatus::available && entry.inference.observedSupportSufficient)
                    {
                        result.diagnosticIndices.push_back(index);
                        if (selectRepresentatives)
                        {
                            auto coefficients = fit.locals.front();
                            coefficients.intercept = entry.inference.properties.front().correctedIntercept;
                            selectionFits.push_back({index, coefficients, fit.sharedSlope});
                        }
                    }
                });
            }
            auto groups = collectUnresolvedGroups(search.windows, result.diagnosticIndices, selectRepresentatives);
            if (!selectRepresentatives)
            {
                result.diagnosticGroups = std::move(groups);
            }
            if (traces.size() == 1)
            {
                auto &classifiedGroups = selectRepresentatives ? groups : result.diagnosticGroups;
                result.apexGroups.reserve(classifiedGroups.size());
                for (auto &group : classifiedGroups)
                {
                    auto apices = classifyApices(traces.front().intensity, group, search);
                    if (selectRepresentatives)
                    {
                        for (const auto &maximum : apices.maxima)
                        {
                            result.apexSelections.push_back(selectRepresentative(traces.front().intensity, maximum,
                                search, selectionFits));
                        }
                        if (group.windowIndices.size() == 1)
                        {
                            if (!apices.unassignedWindowIndices.empty())
                            {
                                result.unassignedSingletons.push_back(group.windowIndices.front());
                            }
                            continue;
                        }
                        result.diagnosticGroups.push_back(std::move(group));
                    }
                    result.apexGroups.push_back(std::move(apices));
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

        CandidateSearchPlan::CandidateSearchPlan(size_t length, std::vector<CandidateModelWindows> models, double familyAlpha)
            : length_(length), models_(std::move(models))
        {
            if (length < 5 || models_.empty() || !(familyAlpha > 0 && familyAlpha < 1))
            {
                throw std::invalid_argument("Invalid candidate search family");
            }
            std::sort(models_.begin(), models_.end(), [](const auto &a, const auto &b) { return a.model < b.model; });
            for (size_t i = 0; i < models_.size(); ++i)
            {
                auto &model = models_[i];
                const auto &options = model.options;
                if ((model.model != PeakModelKind::PiecewiseQuadratic && model.model != PeakModelKind::Hyperbolic)
                    || (i > 0 && models_[i - 1].model == model.model) || model.windows.empty()
                    || !(model.fixedStep > 0 && std::isfinite(model.fixedStep))
                    || (model.model == PeakModelKind::PiecewiseQuadratic && model.fixedStep != 1)
                    || !(options.minimumArea > 0 && std::isfinite(options.minimumArea))
                    || !(options.minimumCoverage > 0 && options.minimumCoverage < 1)
                    || !(options.logAlpha > 0 && options.logAlpha < .5))
                {
                    throw std::invalid_argument("Invalid model in candidate search family");
                }
                std::sort(model.windows.begin(), model.windows.end(), windowLess);
                for (size_t j = 0; j < model.windows.size(); ++j)
                {
                    const auto &window = model.windows[j];
                    if (window.scale < 2 || window.center >= length || window.scale > window.center
                        || window.scale > length - 1 - window.center
                        || (model.model == PeakModelKind::PiecewiseQuadratic && window.scale > maxscale_global)
                        || (j > 0 && !windowLess(model.windows[j - 1], window)))
                    {
                        throw std::invalid_argument("Invalid or duplicate candidate search window");
                    }
                }
                if (model.windows.size() > std::numeric_limits<size_t>::max() - familySize_)
                {
                    throw std::invalid_argument("Candidate search family size overflow");
                }
                familySize_ += model.windows.size();
            }
            windowAlpha_ = familyAlpha / static_cast<double>(familySize_);
            for (auto &model : models_)
            {
                const double minimumP = 1 / (1 + static_cast<double>(model.options.background.rotations));
                if (minimumP > windowAlpha_)
                {
                    throw std::invalid_argument("Rotation budget cannot resolve the Bonferroni window alpha");
                }
                model.options.background.alpha = windowAlpha_;
            }
        }

        CandidateSearchResult CandidateSearchPlan::evaluate(const float *intensities, size_t length) const
        {
            if (length != length_ || intensities == nullptr || !std::all_of(intensities, intensities + length_,
                [](float value) { return std::isfinite(value); }))
            {
                throw std::invalid_argument("Candidate search requires a matching length and finite raw intensities");
            }
            CandidateSearchResult result;
            result.familySize = familySize_;
            result.backgroundAlphaPerWindow = windowAlpha_;
            result.windows.reserve(familySize_);
            for (const auto &model : models_)
            {
                for (const auto &window : model.windows)
                {
                    CandidateSearchWindowResult entry{model.model, window, {}};
                    entry.validation = validatePeakCandidate(model.model, intensities + window.center - window.scale,
                        2 * window.scale + 1, model.fixedStep, model.options);
                    const auto &background = entry.validation.background;
                    if (std::isfinite(background.constant.pValue) && std::isfinite(background.linear.pValue))
                    {
                        entry.adjustedBackgroundP = std::min(1.0, static_cast<double>(familySize_)
                            * std::max(background.constant.pValue, background.linear.pValue));
                    }
                    if (entry.validation.status == PeakValidationStatus::passed)
                    {
                        result.candidateIndices.push_back(result.windows.size());
                    }
                    result.windows.push_back(std::move(entry));
                }
            }
            return result;
        }

        CandidateSearchResult findPeaks(const CandidateSearchPlan &plan, const float *intensities, size_t length)
        {
            return plan.evaluate(intensities, length);
        }

        MeasuredHyperbolicPeaks mapSelectedHyperbolicPeaks(
            const SegmentedHyperbolicSearchResult &result, const AffineSampleAxis &axis)
        {
            const auto &search = result.search;
            if (search.traceCount != 1 || !(search.fixedStep > 0 && std::isfinite(search.fixedStep))
                || !std::isfinite(axis.origin) || !std::isfinite(axis.unitsPerSample) || axis.unitsPerSample == 0)
            {
                throw std::invalid_argument("Measurement output requires one trace and a finite affine sample axis");
            }
            const double unitsPerU = axis.unitsPerSample / search.fixedStep;
            if (!std::isfinite(unitsPerU) || unitsPerU == 0)
            {
                throw std::overflow_error("Measurement step is outside numeric range");
            }
            MeasuredHyperbolicPeaks output;
            output.axisUnit = axis.unit;
            for (const auto &selection : result.apexSelections)
            {
                if (selection.status != ApexSelectionStatus::selected)
                {
                    continue;
                }
                const auto &entry = search.windows.at(selection.windowIndex);
                if (entry.inference.properties.size() != 1 || entry.window.center >= search.length)
                {
                    throw std::invalid_argument("Selected fit does not match the single-trace search");
                }
                const double center = std::fma(static_cast<double>(entry.window.center), axis.unitsPerSample, axis.origin);
                if (!std::isfinite(center))
                {
                    throw std::overflow_error("Measurement center is outside numeric range");
                }
                const auto mapped = mapHyperbolicAffine(entry.inference, {center, unitsPerU, axis.unit});
                const auto &values = mapped.values.front();
                MeasuredHyperbolicPeak peak{selection, values[0], values[1], values[2], {}};
                if (peak.area == 0 && entry.inference.properties.front().areaU.value != 0)
                {
                    throw std::underflow_error("Measurement area underflowed");
                }
                for (size_t i = 0; i < peak.covariance.size(); ++i)
                {
                    peak.covariance[i] = mapped.covariance[i];
                    if (peak.covariance[i] == 0 && entry.inference.covariance[i] != 0)
                    {
                        throw std::underflow_error("Measurement covariance underflowed");
                    }
                }
                output.peaks.push_back(peak);
            }
            return output;
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
                    auto result = request.hyperbolicPlan->evaluateSegmented(
                        {{request.intensity, segmentation.coordinates, request.length}},
                        segmentation.coordinateStep, segmentation.relativeTolerance, segmentation.threshold,
                        segmentation.detail, segmentation.selectRepresentatives || segmentation.measurementAxis != nullptr);
                    if (segmentation.measurementAxis != nullptr)
                    {
                        result.measurementPeaks = mapSelectedHyperbolicPeaks(result, *segmentation.measurementAxis);
                    }
                    return result;
                }
                return request.hyperbolicPlan->evaluate({request.intensity});
            }
            throw std::invalid_argument("Unknown peak model");
        }
    }
}
