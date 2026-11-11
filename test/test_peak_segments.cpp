#include "qalgorithms_peak_search.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <numeric>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

static void require(bool condition, const char *message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

int main()
{
    try
    {
        constexpr double e = 2.7182818284590452;
        const float below = static_cast<float>(e);
        const float above = std::nextafter(below, std::numeric_limits<float>::infinity());
        const std::vector<float> raw = {above, above, above, above, above, 0, above, below, -1,
            std::numeric_limits<float>::quiet_NaN(), above, above, above, above, above};
        std::vector<double> x(23);
        std::iota(x.begin(), x.end(), 0.0);
        const auto spans = segmentIntensities({{raw.data(), x.data(), raw.size()}}, 1, 0, e, 5).spans;
        require(spans.size() == 7 && spans.front().end == 5 && spans.back().end == raw.size(), "Complete leading and trailing spans");
        require(spans[2].status == IntensitySpanStatus::tooShort && spans[3].status == IntensitySpanStatus::belowThreshold
            && spans[4].status == IntensitySpanStatus::negative && spans[5].status == IntensitySpanStatus::nonfinite,
            "Short, below-threshold and malformed inputs remain distinguishable");
        const float exact[] = {3, 4};
        require(segmentIntensities({{exact, x.data(), 2}}, 1, 0, 3, 1).spans[0].status
            == IntensitySpanStatus::belowThreshold, "Threshold equality excluded");

        std::vector<float> a(23, 0), b(23, 0);
        HyperbolicSearchOptions options;
        options.selection = WindowSelection::prespecified;
        options.totalFamilySize = 9;
        const HyperbolicSearchPlan plan(23, .2, {{6, 2}, {11, 4}, {17, 2}}, options);
        const std::vector<AlignedIntensityTrace> traces = {{a.data(), x.data(), 23}, {b.data(), x.data(), 23}};
        auto result = plan.evaluateSegmented(traces, 1);
        require(result.search.windows.size() == 3 && result.search.familySize == 9
            && result.segmentation.spans.size() == 1 && result.windowStatus[1] == SegmentedWindowStatus::outsideEligibleSpan,
            "All-zero traces retain declared family");
        for (size_t i = 5; i <= 17; ++i)
        {
            const double u = .2 * (static_cast<double>(i) - 11);
            const double noise = .01 * std::sin(1.7 * i);
            if (i <= 16)
            {
                a[i] = static_cast<float>(std::exp(8 + .3 * u - 2 * std::cosh(u) + noise));
            }
            if (i >= 7)
            {
                b[i] = static_cast<float>(std::exp(9 + .3 * u - 2.2 * std::cosh(u) - .15 * std::sinh(u) + noise));
            }
        }
        result = plan.evaluateSegmented(traces, 1);
        require(result.windowStatus[0] == SegmentedWindowStatus::outsideEligibleSpan
            && result.windowStatus[1] == SegmentedWindowStatus::evaluated
            && result.windowStatus[2] == SegmentedWindowStatus::outsideEligibleSpan
            && result.search.candidateIndices.empty(), "Joint intersection and exploratory decision");
        HyperbolicInferenceOptions referenceOptions;
        referenceOptions.familySize = 9;
        const auto expected = inferHyperbolicWindow({a.data() + 7, b.data() + 7}, 10, 4, 4, .2, referenceOptions);
        const auto &actual = result.search.windows[1].inference;
        require(actual.status == InferenceStatus::available && actual.covariance == expected.covariance
            && actual.logAdjustedP == expected.logAdjustedP && actual.residualSumSquares == expected.residualSumSquares,
            "Exact direct-slice covariance and statistics");
        for (size_t j = 0; j < 2; ++j)
        {
            const auto &p = actual.properties[j];
            const auto &q = expected.properties[j];
            require(p.apexU.value == q.apexU.value && p.height.value == q.height.value
                && p.correctedIntercept == q.correctedIntercept && p.areaU.value == q.areaU.value
                && p.areaUncertainty.value == q.areaUncertainty.value, "Exact apex, refit, area and SE");
        }
        // Independently compare the unchanged running fit on the retained
        // common span to a direct fit over just the requested window.
        const auto span = result.segmentation.spans[1];
        require(span.start == 7 && span.end == 17, "Original common-span coordinates");
        std::vector<float> logsA, logsB;
        for (size_t i = span.start; i < span.end; ++i)
        {
            logsA.push_back(static_cast<float>(std::log(static_cast<double>(a[i]))));
            logsB.push_back(static_cast<float>(std::log(static_cast<double>(b[i]))));
        }
        HyperbolicModel::Fit retained, sliced;
        const HyperbolicModel::Geometry geometry(4, .2);
        HyperbolicModel::runningFit({logsA.data(), logsB.data()}, logsA.size(), geometry, [&](const auto &fit)
        {
            if (fit.center == 4 && fit.scale == 4)
            {
                retained = fit;
            }
        });
        HyperbolicModel::runningFit({logsA.data(), logsB.data()}, 9, geometry, [&](const auto &fit)
        {
            if (fit.center == 4 && fit.scale == 4)
            {
                sliced = fit;
            }
        });
        require(retained.locals.size() == 2 && sliced.locals.size() == 2
            && retained.sharedSlope == sliced.sharedSlope, "Exact shared coefficient");
        for (size_t j = 0; j < 2; ++j)
        {
            require(retained.locals[j].intercept == sliced.locals[j].intercept
                && retained.locals[j].curvature == sliced.locals[j].curvature
                && retained.locals[j].skew == sliced.locals[j].skew, "Exact trace-specific coefficients");
        }
        a[11] = 0;
        require(plan.evaluateSegmented(traces, 1).windowStatus[1] == SegmentedWindowStatus::outsideEligibleSpan,
            "Internal zero splits the window");
        auto shifted = x;
        shifted[7] += .1;
        bool rejected = false;
        try
        {
            (void)plan.evaluateSegmented({traces[0], {b.data(), shifted.data(), 23}}, 1);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        require(rejected, "Misaligned grid rejected");
        const std::vector<float> positive(23, 100);
        const auto full = plan.evaluateSegmented({{positive.data(), x.data(), 23}}, 1);
        for (auto status : full.windowStatus)
        {
            require(status == SegmentedWindowStatus::evaluated, "All positive windows evaluated");
        }
        for (size_t i = 11; i < x.size(); ++i)
        {
            x[i] += 2;
        }
        const auto gap = plan.evaluateSegmented({{positive.data(), x.data(), 23}}, 1);
        require(gap.segmentation.coordinateBreaks == std::vector<size_t>{11}
            && gap.windowStatus[1] == SegmentedWindowStatus::outsideEligibleSpan, "No fit across a coordinate gap");

        // Exercise the common entry point with an offset peak, zeros and a
        // coordinate break inside the positive span.
        std::iota(x.begin(), x.end(), 0.0);
        for (size_t i = 17; i < x.size(); ++i)
        {
            x[i] += 2;
        }
        const HyperbolicSearchPlan dispatchPlan(23, .2, {{11, 4}, {11, 5}, {15, 3}}, options);
        const PeakSearchRequest::Segmentation segmentation{x.data(), 1};
        PeakSearchRequest request;
        request.model = PeakModelKind::Hyperbolic;
        request.intensity = b.data();
        request.length = b.size();
        request.hyperbolicPlan = &dispatchPlan;
        request.segmentation = &segmentation;
        // The wider window starts at sample 6. Restore that observation so
        // both fits describe the same positive span before the break.
        b[6] = static_cast<float>(std::exp(9 - .3 - 2.2 * std::cosh(1.0) + .15 * std::sinh(1.0)));
        const auto dispatched = std::get<SegmentedHyperbolicSearchResult>(findPeaks(request));
        require(dispatched.search.options.selection == WindowSelection::exploratory
            && dispatched.search.candidateIndices.empty() && dispatched.search.unresolvedGroups.empty()
            && dispatched.search.familySize == 9, "Adaptive segmentation retains exploratory inference");
        require(dispatched.windowStatus[2] == SegmentedWindowStatus::outsideEligibleSpan
            && dispatched.diagnosticIndices == std::vector<size_t>{0, 1}
            && dispatched.diagnosticGroups.size() == 1, "Supported overlapping fits are retained as diagnostics");
        const auto &group = dispatched.diagnosticGroups.front();
        require(group.start == 6 && group.end == 16 && group.windowIndices == std::vector<size_t>{0, 1},
            "Diagnostic group refers to original samples and canonical windows");
        for (size_t index : group.windowIndices)
        {
            const auto &entry = dispatched.search.windows[index];
            const auto &apex = entry.inference.properties.front().apexU;
            require(apex.status == PropertyStatus::Available && std::abs(entry.window.center + apex.value / .2 - 11) < 1,
                "Apex can be mapped back to the original sample grid");
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
