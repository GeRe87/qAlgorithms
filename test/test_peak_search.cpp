#include "qalgorithms_peak_search.h"

#include "cephes.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

namespace
{
    void require(bool condition, const char *message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }

    void testModelDispatch()
    {
        constexpr size_t length = 41;
        std::vector<float> raw(length), coordinates(length);
        for (size_t i = 0; i < length; ++i)
        {
            const double u = .2 * (static_cast<double>(i) - 20);
            raw[i] = static_cast<float>(5 + 800 * std::exp(-u * u / 2) * (1 + .01 * std::sin(2 * u)));
            coordinates[i] = static_cast<float>(100 + .03 * i);
        }
        PeakSearchRequest request;
        request.intensity = raw.data();
        request.coordinates = coordinates.data();
        request.length = length;
        for (bool collective : {false, true})
        {
            request.describeUnresolved = collective;
            std::vector<RegressionGauss> peaks;
            std::vector<UnresolvedPeak> groups;
            const auto count = qpeaks_find(raw.data(), coordinates.data(), nullptr, length, &peaks,
                collective ? &groups : nullptr);
            const auto result = std::get<PiecewiseQuadraticSearchResult>(findPeaks(request));
            require(count > 0 && count == result.status && peaks.size() == result.peaks.size(),
                "Default model uses the quadratic detector");
            require(groups.size() == result.unresolvedGroups.size(), "Legacy group reporting retained");
            for (size_t i = 0; i < peaks.size(); ++i)
            {
                const auto &a = peaks[i];
                const auto &b = result.peaks[i];
                require(a.position == b.position && a.area == b.area && a.height == b.height
                    && a.coeffs.x0 == b.coeffs.x0 && a.coeffs.scale == b.coeffs.scale && a.isValid == b.isValid,
                    "Selected legacy results preserved in order");
            }
        }
        HyperbolicSearchOptions options;
        options.selection = WindowSelection::prespecified;
        const HyperbolicSearchPlan plan(length, .2, {{20, 7}, {20, 8}}, options);
        request.model = PeakModelKind::Hyperbolic;
        request.hyperbolicPlan = &plan;
        const auto reference = plan.evaluate({raw.data()});
        const auto result = std::get<HyperbolicSearchResult>(findPeaks(request));
        require(result.familySize == 2 && result.windows.size() == reference.windows.size()
            && result.candidateIndices == reference.candidateIndices
            && result.unresolvedGroups.size() == reference.unresolvedGroups.size(),
            "Hyperbolic dispatch preserves the declared search result");
        for (size_t i = 0; i < result.windows.size(); ++i)
        {
            const auto &a = result.windows[i];
            const auto &b = reference.windows[i];
            require(a.window.center == b.window.center && a.window.scale == b.window.scale
                && a.inference.status == b.inference.status && a.inference.groupDecision == b.inference.groupDecision
                && a.unresolvedGroup == b.unresolvedGroup, "Hyperbolic window order and decisions retained");
        }
        --request.length;
        bool threw = false;
        try
        {
            (void)findPeaks(request);
        }
        catch (const std::invalid_argument &)
        {
            threw = true;
        }
        require(threw, "Mismatched hyperbolic plan rejected");
        request.model = static_cast<PeakModelKind>(99);
        threw = false;
        try
        {
            (void)findPeaks(request);
        }
        catch (const std::invalid_argument &)
        {
            threw = true;
        }
        require(threw, "Unknown peak model rejected");
    }

    std::vector<std::vector<size_t>> referenceGroups(const HyperbolicSearchResult &result)
    {
        // Independent all-pairs overlap graph + traversal, rather than the
        // production interval sweep. Compare transitive components explicitly.
        std::vector<bool> visited(result.windows.size(), false);
        std::vector<std::vector<size_t>> groups;
        for (size_t seed : result.candidateIndices)
        {
            if (visited[seed])
            {
                continue;
            }
            std::vector<size_t> members = {seed};
            visited[seed] = true;
            for (size_t k = 0; k < members.size(); ++k)
            {
                const auto &a = result.windows[members[k]].window;
                for (size_t other : result.candidateIndices)
                {
                    const auto &b = result.windows[other].window;
                    if (!visited[other] && a.center - a.scale <= b.center + b.scale
                        && b.center - b.scale <= a.center + a.scale)
                    {
                        visited[other] = true;
                        members.push_back(other);
                    }
                }
            }
            if (members.size() > 1)
            {
                std::sort(members.begin(), members.end());
                groups.push_back(std::move(members));
            }
        }
        return groups;
    }

    void testFamilyAndAmbiguity()
    {
        constexpr size_t length = 77;
        std::vector<float> raw(length);
        for (size_t i = 0; i < length; ++i)
        {
            // Two well separated identical shapes, no fitted baseline.
            const double u = .12 * (static_cast<double>(i) - (i < 38 ? 18 : 58));
            raw[i] = static_cast<float>(std::exp(7 + .3 * u - 2 * std::cosh(u) + .01 * std::sin(i * 1.7)));
        }
        HyperbolicSearchOptions options;
        options.selection = WindowSelection::prespecified;
        const std::vector<PeakSearchWindow> declared = {{58, 8}, {18, 10}, {17, 8}, {18, 8}, {10, 2}, {59, 8}, {3, 2}};
        HyperbolicSearchPlan plan(length, .12, declared, options);
        const auto result = plan.evaluate({raw.data()});
        require(result.familySize == declared.size() && result.windows.size() == declared.size(), "Complete declared denominator and all windows retained");
        for (size_t i = 0; i < result.windows.size(); ++i)
        {
            const auto &entry = result.windows[i];
            HyperbolicInferenceOptions single;
            single.selection = options.selection;
            single.familySize = declared.size();
            const auto reference = inferHyperbolicWindow({raw.data()}, length, entry.window.center, entry.window.scale, .12, single);
            require(entry.inference.status == reference.status && entry.inference.groupDecision == reference.groupDecision,
                "Batch and fixed-window status equivalence");
            require(entry.inference.covariance == reference.covariance, "Batch preserves complete fixed-window covariance");
            if (std::isfinite(reference.logNominalP))
            {
                require(entry.inference.logAdjustedP == reference.logAdjustedP, "Complete-family adjustment retained for invalid as well as valid shapes");
            }
        }
        require(result.candidateIndices.size() >= 4, "Several alternative strong fits retained");
        const auto expected = referenceGroups(result);
        require(expected.size() == result.unresolvedGroups.size() && expected.size() == 2, "Separate unresolved components");
        for (size_t group = 0; group < expected.size(); ++group)
        {
            require(expected[group] == result.unresolvedGroups[group].windowIndices, "Sweep vs independent overlap graph");
            for (size_t index : expected[group])
            {
                require(result.windows[index].unresolvedGroup == group, "Unresolved membership indexes complete window list");
            }
        }
        auto reversed = declared;
        std::reverse(reversed.begin(), reversed.end());
        const auto reordered = HyperbolicSearchPlan(length, .12, reversed, options).evaluate({raw.data()});
        require(reordered.candidateIndices == result.candidateIndices && referenceGroups(reordered) == expected,
            "Input order does not choose a winning window");
        options.totalFamilySize = 100;
        const auto larger = HyperbolicSearchPlan(length, .12, declared, options).evaluate({raw.data()});
        require(larger.familySize == 100, "External steps/trace pools can extend the declared family");
        for (size_t i = 0; i < result.windows.size(); ++i)
        {
            if (std::isfinite(result.windows[i].inference.logNominalP))
            {
                require(larger.windows[i].inference.logAdjustedP >= result.windows[i].inference.logAdjustedP,
                    "Adding declared tests cannot increase significance");
            }
        }
        const auto exploratory = HyperbolicSearchPlan(length, .12, declared).evaluate({raw.data()});
        require(exploratory.candidateIndices.empty() && exploratory.unresolvedGroups.empty(), "Exploratory search has no conditional detections");
    }

    void testTransitiveAndTouchingOverlap()
    {
        HyperbolicSearchOptions options;
        options.selection = WindowSelection::prespecified;
        for (size_t spacing : {10u, 14u})
        {
            const size_t scale = spacing / 2;
            const size_t length = 2 * scale + 2 * spacing + 1;
            std::vector<float> raw(length);
            std::vector<PeakSearchWindow> windows;
            for (size_t j = 0; j < 3; ++j)
            {
                windows.push_back({scale + j * spacing, scale});
            }
            for (size_t i = 0; i < length; ++i)
            {
                size_t nearest = windows.front().center;
                for (const auto &window : windows)
                {
                    if (std::abs(static_cast<double>(i) - window.center) < std::abs(static_cast<double>(i) - nearest))
                    {
                        nearest = window.center;
                    }
                }
                const double u = .12 * (static_cast<double>(i) - nearest);
                raw[i] = static_cast<float>(std::exp(7 - 2 * std::cosh(u) + .001 * std::sin(i * 1.7)));
            }
            const auto result = HyperbolicSearchPlan(length, .12, windows, options).evaluate({raw.data()});
            if (result.candidateIndices.size() != 3)
            {
                std::cerr << "Chain spacing " << spacing << ": ";
                for (const auto &entry : result.windows)
                {
                    std::cerr << static_cast<int>(entry.inference.status) << "/" << static_cast<int>(entry.inference.groupDecision)
                              << "/" << entry.inference.adjustedP << " ";
                }
                std::cerr << '\n';
            }
            require(result.candidateIndices.size() == 3, "Three distinct retained hypotheses exercise a chain");
            require(result.unresolvedGroups.size() == 1 && result.unresolvedGroups[0].windowIndices.size() == 3,
                "Transitive/touching interval component retains every hypothesis");
            require(windows.front().center + scale < windows.back().center - scale, "First and last support do not overlap directly");
            require(referenceGroups(result)[0] == result.unresolvedGroups[0].windowIndices, "Chain vs independent graph");
            const auto singleton = HyperbolicSearchPlan(length, .12, {windows.front()}, options).evaluate({raw.data()});
            require(singleton.candidateIndices.size() == 1 && singleton.unresolvedGroups.empty(), "Single hypothesis is not assigned an overlap group");
        }
    }

    void testPlanValidationAndGrid()
    {
        const auto plan = makeUniformHyperbolicSearchPlan(17, .12, 2, 5, 3);
        std::vector<PeakSearchWindow> expected;
        for (size_t center = 2; center < 15; center += 3)
        {
            for (size_t scale = 2; scale <= 5; ++scale)
            {
                if (scale <= center && scale < 17 - center)
                {
                    expected.push_back({center, scale});
                }
            }
        }
        require(expected.size() == plan.familySize(), "Uniform family count independently enumerated");
        for (size_t i = 0; i < expected.size(); ++i)
        {
            require(expected[i].center == plan.windows()[i].center && expected[i].scale == plan.windows()[i].scale, "Uniform center/scale order");
        }
        auto rejected = [](auto operation)
        {
            bool caught = false;
            try
            {
                operation();
            }
            catch (const std::invalid_argument &)
            {
                caught = true;
            }
            require(caught, "Invalid complete search request rejected");
        };
        rejected([] { (void)HyperbolicSearchPlan(17, .1, {{8, 2}, {8, 2}}); });
        rejected([] { (void)HyperbolicSearchPlan(17, .1, {}); });
        rejected([] { (void)HyperbolicSearchPlan(17, .1, {{0, 2}}); });
        rejected([] { (void)makeUniformHyperbolicSearchPlan(17, .1, 2, 5, 0); });
        HyperbolicSearchOptions options;
        options.totalFamilySize = 1;
        rejected([&] { (void)HyperbolicSearchPlan(17, .1, {{8, 2}, {8, 3}}, options); });
        std::vector<float> raw(17, 100);
        raw.back() = std::numeric_limits<float>::quiet_NaN();
        // Invalid input is outside the declared windows but still detected by
        // the whole-trace prescan, before any fit can be returned.
        rejected([&] { (void)HyperbolicSearchPlan(17, .1, {{8, 2}}).evaluate({raw.data()}); });
        rejected([&] { (void)plan.evaluate({}); });
        rejected([&] { (void)plan.evaluate({nullptr}); });
    }

    void testOverlappingNullFamily()
    {
        // The same independent observations are reused by overlapping windows;
        // the tests are dependent. Compare the family event to an independent
        // incomplete-beta calculation on the reported fixed-window F values.
        HyperbolicSearchOptions options;
        options.selection = WindowSelection::prespecified;
        const auto plan = makeUniformHyperbolicSearchPlan(25, .16, 3, 6, 2, options);
        std::mt19937 random(219831);
        std::normal_distribution<double> noise(0, .04);
        constexpr size_t trials = 800;
        size_t anyNominal = 0, anyCorrected = 0, anyCandidate = 0;
        for (size_t trial = 0; trial < trials; ++trial)
        {
            std::vector<float> raw(25);
            for (size_t i = 0; i < raw.size(); ++i)
            {
                raw[i] = static_cast<float>(std::exp(6 + .1 * i + noise(random)));
            }
            const auto result = plan.evaluate({raw.data()});
            bool nominal = false, corrected = false;
            for (const auto &entry : result.windows)
            {
                const auto &fit = entry.inference;
                if (!std::isfinite(fit.comparisonF))
                {
                    continue;
                }
                const double df = static_cast<double>(fit.residualDegreesOfFreedom);
                const double x = df / (df + 2 * fit.comparisonF);
                const double p = cephes::incbeta(df / 2, 1, x);
                nominal = nominal || p <= .05;
                corrected = corrected || p <= .05 / plan.familySize();
                if (std::isfinite(fit.adjustedP))
                {
                    require(std::abs(fit.adjustedP - std::min(1.0, p * plan.familySize())) < 2e-6,
                        "Search correction vs independent incomplete-beta reference");
                }
            }
            anyNominal += nominal;
            anyCorrected += corrected;
            anyCandidate += !result.candidateIndices.empty();
            require(result.candidateIndices.empty() || corrected, "Shape gates cannot create a family rejection without corrected evidence");
        }
        const double correctedRate = static_cast<double>(anyCorrected) / trials;
        require(correctedRate < .08 && anyNominal > anyCorrected, "Dependent-window null family check");
        std::cout << "Overlapping null family: windows=" << plan.familySize() << ", trials=" << trials
                  << ", unadjusted=" << static_cast<double>(anyNominal) / trials
                  << ", Bonferroni=" << correctedRate
                  << ", shaped candidates=" << static_cast<double>(anyCandidate) / trials << '\n';
    }
}

int main()
{
    try
    {
        testModelDispatch();
        testFamilyAndAmbiguity();
        testTransitiveAndTouchingOverlap();
        testPlanValidationAndGrid();
        testOverlappingNullFamily();
        std::cout << "Declared search and unresolved-group checks passed\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
