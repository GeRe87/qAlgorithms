#include "qalgorithms_datatypes.h"
#include "qalgorithms_qpeaks.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvariadic-macro-arguments-omitted"

#include "../src/qalgorithms_qpeaks.cpp" // NOLINT
#include "common_test_utils.hpp"

using namespace qAlgorithms;

// @todo test these specific values with qpeaks:

// std::vector<float> intensity = {16649.7441, 34118.1367, 46703.2617, 41025.5352, 53125.8086, 140134.125, 402735.812, 689393.625, 907356.938, 728596.75, 432838.844, 156655.438, 43841.4922, 11888.1475, 3812.54126};
// std::vector<float> intensity = {32, 475, 711, 472, 207, 132, 57, 14};
// std::vector<float> intensity = {157.883072, 1325.722046, 2188.603760, 5415.137695, 12294.484375, 16239.560547, 13575.218750, 9618.787109, 7654.178223, 20002.025391, 69383.062500, 147876.296875, 233001.171875, 244286.796875, 162216.375000, 82337.882812, 23717.978516, 5968.742676, 2921.130859, 945.386047};
// std::vector<float> intensity = {1274.10596, 5653.32959, 19341.3398, 51021.3789, 103777.039, 162754.797, 179871.859, 128027.453, 58688.5547, 17326.6309, 3294.46802};

// NOLINTBEGIN
struct PeakTest
{
    std::vector<float> intensity = {0};
    size_t expect_count = 0;
};

PeakTest pt_01 = {
    // Note: this is a double peak system in need of deconvolution which i included because
    // the old design for selecting the best entry of a set of regressions eliminated the
    // second, smaller peak despite the apexes being far enough apart. Interesting, at the
    // time of writing, no regression with positive coefficients was taken into consideration
    {8862.04883, 17619.8887, 23784.2598, 22516.0684, 17171.5332, 14893.7227, 16483.5371, 16239.6406, 11697.918, 5697.0332},
    2};

PeakTest pt_02 = {
    // Peak with slightly increased baseline on the right tail that proved difficult to detect
    // in beta versions of the program. This was due to the Chi-square filter not functioning for
    // the data after reversing the log transform
    {16427.9434, 34097.4414, 42639.7266, 102040.758, 264363.312, 486761.125, 909302, 783670.375, 404226.625, 194167.844, 81423.1406, 63319.4688, 62018.1602, 27416.2754},
    1};

PeakTest pt_03 = {

    // Data for a malformed bin which should not contain any features, but still allowed for
    // peaks to be constructed: (multi-line)
    {223.15509, 236.694443, 242.834946, 206.72879, 244.555389, 232.282364, 233.314102, 276.77417, 236.010162, 234.881836, 213.153366, 215.19165, 274.305542,
     227.700317, 240.399643, 279.177734, 267.526672, 254.195343, 258.700653, 261.078522, 270.642853, 261.583618, 254.256668, 250.748398, 256.782806, 274.030823,
     232.25827, 214.949051, 245.469711, 284.085175, 230.861053, 260.483185, 266.808777, 247.738373, 230.255249, 209.010635, 234.587158, 253.927856, 209.396133,
     259.263336, 254.414474, 285.63858, 227.385788, 217.065872, 274.531311, 229.829956, 246.372513, 230.911835, 215.58461, 250.233902, 232.271439, 217.880798,
     216.245712, 237.993469, 217.751038, 222.585983, 305.338898, 209.944839, 233.039291, 246.02948, 222.083221, 241.314194, 233.13353, 275.777863, 236.0522,
     236.623123, 269.473999, 227.634735, 233.807373, 215.515015, 292.07547, 241.878036, 247.750275, 246.052521, 226.036819, 233.290924, 191.461838, 226.269302,
     249.507401, 242.797043, 183.200607, 280.93985, 254.288147, 229.075333, 258.677795, 218.803909, 225.112, 216.144104, 248.790497, 236.507645, 228.74855,
     261.456085, 217.647079, 254.565613, 240.973572, 266.561127, 234.288925, 231.6884, 282.011963},
    0};

PeakTest pt_04 = {

    // This data contains three peaks, but only the two less intense ones were found.
    // The new apex grouping function marked them as conflicting, which should not be the case.
    {28913.7539, 15335.7256, 31011.6035, 59799.4961, 30304.0547, 101.776535, 7791.06396, 16737.2637, 7326.88721, 5033.8501, 11136.9746, 12789.5098, 7660.45898, 2526.31592},
    3};

PeakTest pt_05 = {
    // this peak lead to unresolvable regression conflicts in a version of the new regression elimination
    {1566.37317, 6721.91357, 10333.0713, 12963.7285, 18588.1035, 23554.9043, 26401.1895, 27222.6621, 27128.6387, 18343.3789, 9999.44727, 823.102112},
    1};

PeakTest pt_06 = {
    // This peak system produced groups that did not correctly group apexes before eliminating regressions.
    // Some of the offending regressions were poorly fitting and could probably be removed by a goodness-of-fit filter.
    // This problem only occurs for one of the three peaks present in the data. The test case was
    // taken from a failing check of group coherence, but the group is generally correct - the only
    // problem is that no real peaks are describved by the two conflicting apexes.
    {629.24292, 1937.10291, 3530.98169, 4416.47412, 3860.57593, 2283.37793, 1097.28284, 1384.93604,
     2879.35815, 4468.28662, 5404.12451, 5281.40186, 2488.4541, 234.448563, 522.505859, 3366.85596,
     10939.0303, 18904.9043, 24883.1719, 19762.5625, 11631.5967, 4173.91699, 858.461609},
    3};
// NOLINTEND

static int test_qpeaks_find(const PeakTest *test)
{
    const size_t len = test->intensity.size();
    std::vector<RegressionGauss> result;
    std::vector<float> data_x(len, 0);
    for (size_t i = 0; i < len; i++)
        data_x[i] = ((float)(i + 1) * 0.6F);

    size_t count = qpeaks_find(test->intensity.data(),
                               data_x.data(),
                               nullptr,
                               len,
                               &result);
    size_t expect = test->expect_count;
    assert(count == expect, "Expected %zu peaks, found %zu", expect, count);
    return 0;
}

static int test_qpeaks_set(void)
{
    test_qpeaks_find(&pt_06);
    test_qpeaks_find(&pt_05);
    test_qpeaks_find(&pt_04);
    test_qpeaks_find(&pt_03);
    test_qpeaks_find(&pt_02);
    // test_qpeaks_find(&pt_01); // removed until deconvolution is introduced, second apex is only four distinct points
    return 0;
}

static float peakVal_gauss(double x, double apex, double height, double sdev)
{
    return float(height * exp(-(x - apex) * (x - apex) / (2 * sdev * sdev)));
}

static double area_gauss(double height, double sdev)
{
    return height * sdev * sqrt(2 * M_PI);
}

static double fwhm_gauss(double sdev)
{
    // 0.5 * height = height * exp(-(x-apex)^2 / (2 * sdev^2))
    // 0.5 = exp(-(x-apex)^2 / (2 * sdev^2))
    // ln(0.5) * (2 * sdev^2) = -(x-apex)^2
    // sqrt(-ln(0.5) * 2 * sdev^2) = x - apex
    // x = apex + sqrt(-ln(0.5) * 2 * sdev^2)

    // since the gaussian peak is symmetrical, we only need the distance from the apex * 2
    // this means that the apex is not required for knowing the FWHM
    return 2 * sqrt(-log(0.5) * 2 * sdev * sdev);
}

// note on simulating points: it is not possible to just generate arbitrary points
// since the core assumption of qPeaks is that all points are evenly spaced. Furthermore,
// it only makes sense to characterise the algorithm with data where the full peak is captured.
// To make cases where two peaks are in one system or where noise is added easier, all
// peak functions are additive. It is the responsibility of the user to zero them before first use.
static void simulate_gauss(
    const std::vector<float> *xvals,
    double apex, double height, double sdev,
    std::vector<float> *simulated)
{
    assert(!simulated->empty(), "wrong usage of test function\n", NULL);
    assert(xvals->size() == simulated->size(), "wrong usage of test function\n", NULL);

    for (size_t i = 0; i < simulated->size(); i++)
    {
        double x = xvals->at(i);
        // @todo add the option to randomise here
        simulated->at(i) += peakVal_gauss(x, apex, height, sdev);
    }
}

static void control_sim_gauss()
{
    // generate data using a standard gaussian on an equidistant x axis
    float x_start = 100;
    float x_step = 1;
    float apex = 115;
    uint16_t length = 30;
    double sdev = 2.5;
    double height = 1000;

    double fwhm = fwhm_gauss(sdev);
    double area = area_gauss(height, sdev);

    std::vector<float> xvals(length, 0);
    std::vector<float> yvals(length, 0);
    for (uint16_t i = 0; i < length; i++)
    {
        xvals[i] = x_start + x_step * (float)i;
    }

    simulate_gauss(&xvals, apex, height, sdev, &yvals);

    std::vector<RegressionGauss> ret;
    qpeaks_find(yvals.data(), xvals.data(), nullptr, length, &ret);

    assert(!ret.empty(), "Peak not found\n", NULL);
    assert(ret.size() == 1, "Too many peaks found\n", NULL);

    RegressionGauss reg = ret.front();

    float apex_p = reg.position;
    float height_p = reg.height;
    float fwhm_p = reg.fwhm;
    float area_p = reg.area; // wrong result

    float area_e = (float)area_empiric(&xvals, &yvals);
    // RegCoeffs c = reg.coeffs;
    // float area_c = peakArea(c.b0, c.b1, c.b2, c.b3, x_step);

    assert(flt_equal(apex, apex_p, FLT_EPSILON), "inaccurate position\n", NULL);
    assert(flt_equal(height, height_p, reg.height_unc), "inaccurate height\n", NULL);
    assert(flt_equal(fwhm, fwhm_p, 10e-4), "inaccurate width\n", NULL);
    assert(flt_equal(area, area_p, 0.01), "inaccurate area (%f vs. %f), empiric %f\n", area, area_p, area_e);
}

static void test_apex_conversion_nonunit_spacing()
{
    // A fitted apex at sample index 2.25 must be interpolated in the x-axis
    // units. These two axes have steps smaller and larger than one.
    const float axes[2][5] = {{10, 10.25F, 10.5F, 10.75F, 11},
                              {10, 12, 14, 16, 18}};
    const float expected_positions[2] = {10.5625F, 14.5F};
    const float expected_uncertainties[2] = {0.05F, 0.4F};
    for (size_t i = 0; i < 2; ++i)
    {
        RegressionGauss peak;
        // log y = log(100) - (x - 0.25)^2, with x relative to sample 2.
        peak.coeffs = {log(100.0) - 0.0625, 0.5, -1, -1, 2, 2};
        peak.position = 2.25F;
        peak.position_unc = 0.2F;
        peak.height = 100;
        peak.area = 100;
        peak.area_unc = 1;
        peak.isValid = true;

        adjustRegression(&peak, axes[i]);

        assert(std::abs(peak.position - expected_positions[i]) < 1e-6F,
               "Apex conversion failed: expected %f, found %f\n",
               expected_positions[i], peak.position);
        assert(std::abs(peak.position_unc - expected_uncertainties[i]) < 1e-6F,
               "Position uncertainty must use the same x-axis units\n", NULL);
    }
}

static void test_weighted_mean_ignores_interpolated_values()
{
    // Zero is the missing-value marker. Its large interpolated intensity must
    // not dominate the measured masses: (1*100 + 3*104)/(1+3) = 103.
    const float weights[] = {1, 1000, 3};
    const float masses[] = {100, 0, 104};
    float uncertainty = 0;
    const float mean = weightedMeanAndVariance_EIC(weights, masses, 3, &uncertainty);

    assert(std::abs(mean - 103.0F) < 1e-6F,
           "Incorrect weighted mass: expected 103, found %f\n", mean);
    // The function reports sqrt(weighted MSD / measured point count):
    // sqrt((1*9 + 3*1) / (4*2)) = sqrt(1.5).
    assert(std::abs(uncertainty - sqrt(1.5)) < 1e-6,
           "Interpolated values changed the reported uncertainty: %f\n", uncertainty);
    const float mean_only = weightedMeanAndVariance_EIC(weights, masses, 3, nullptr);
    assert(std::abs(mean_only - 103.0F) < 1e-6F,
           "Mean-only path must use the same measured points\n", NULL);
}

static void test_peak_at_region_boundaries_and_cumulative_df()
{
    // Exactly five points require the final sample to participate. Both ends
    // are above e, so there is no below-threshold sentinel inside the region.
    const size_t length = 5;
    std::vector<float> x = {20, 20.5F, 21, 21.5F, 22};
    std::vector<float> y(length, 0);
    simulate_gauss(&x, 21.05, 1000, 0.4, &y);
    assert(y.front() > 3 && y.back() > 3,
           "Boundary regression fixture must start and end above the threshold\n", NULL);

    std::vector<RegressionGauss> peaks;
    const int32_t count = qpeaks_find(y.data(), x.data(), nullptr, length, &peaks);
    assert(count == 1, "Expected one complete boundary peak, found %d\n", count);
    assert(peaks.front().span.startIdx == 0 && peaks.front().span.length == 5,
           "The complete five-point region must be fitted\n", NULL);
    assert(std::abs(peaks.front().position - 21.05F) < 1e-3F,
           "Incorrect physical apex for boundary peak: %f\n", peaks.front().position);
    assert(peaks.front().df == 1, "Five measured points minus four fitted coefficients give one DF\n", NULL);

    // Prefixing three low observations must shift only the source span. Their
    // cumulative counts must not inflate the fitted peak's degrees of freedom.
    std::vector<float> prefixed_x = {18.5F, 19, 19.5F};
    std::vector<float> prefixed_y = {1, 1, 1};
    prefixed_x.insert(prefixed_x.end(), x.begin(), x.end());
    prefixed_y.insert(prefixed_y.end(), y.begin(), y.end());
    const uint16_t cumulative_df[] = {1, 2, 3, 4, 5, 6, 7, 8};
    peaks.clear();
    const int32_t prefixed_count = qpeaks_find(prefixed_y.data(), prefixed_x.data(),
                                               cumulative_df, prefixed_y.size(), &peaks);
    assert(prefixed_count == 1, "Expected the same peak after a prefix, found %d\n", prefixed_count);
    assert(peaks.front().span.startIdx == 3 && peaks.front().span.length == 5,
           "Peak span must refer to the original input indices\n", NULL);
    assert(peaks.front().df == 1,
           "Prefix observations must not add DF to the fitted region: %u\n", peaks.front().df);
    assert(std::abs(peaks.front().position - 21.05F) < 1e-3F,
           "Prefix observations shifted the physical peak position\n", NULL);
}

static void test_apex_grouping_more_than_128_candidates()
{
    std::vector<RegressionGauss> candidates(300);
    for (RegressionGauss &candidate : candidates)
    {
        candidate.position = 7.5F;
        candidate.span = {2, 11};
    }
    std::vector<int16_t> groups(candidates.size(), 77);
    const size_t count = groupRegsByApex(&candidates, groups.data());
    assert(count == 1, "Coincident candidates must form exactly one apex group\n", NULL);
    for (size_t i = 0; i < groups.size(); ++i)
        assert(groups[i] == 0, "Candidate %zu was not assigned to the common group\n", i);
}

static void test_apex_grouping_order_and_reflection()
{
    // A overlaps B, B overlaps C; A does not overlap C. D is separate.
    // Every permutation and its reflection must preserve this partition.
    std::vector<RegressionGauss> original(4);
    original[0].position = 4; original[0].span = {0, 6};
    original[1].position = 6; original[1].span = {3, 6};
    original[2].position = 8; original[2].span = {6, 6};
    original[3].position = 20; original[3].span = {17, 7};
    for (const bool reflected : {false, true})
    {
        std::vector<size_t> order = {0, 1, 2, 3};
        do
        {
            std::vector<RegressionGauss> candidates;
            for (const size_t index : order)
            {
                auto candidate = original[index];
                if (reflected)
                {
                    candidate.position = 30 - candidate.position;
                    candidate.span.startIdx = 31 - candidate.span.startIdx - candidate.span.length;
                }
                candidates.push_back(candidate);
            }
            int16_t groups[4];
            assert(groupRegsByApex(&candidates, groups) == 2,
                   "Connected windows must form two groups in every order", NULL);
            for (size_t i = 0; i < 4; ++i)
                for (size_t j = 0; j < 4; ++j)
                    assert((groups[i] == groups[j]) == ((order[i] == 3) == (order[j] == 3)),
                           "Apex partition changed under permutation/reflection", NULL);
        } while (std::next_permutation(order.begin(), order.end()));
    }
    std::vector<RegressionGauss> singleton(1);
    int16_t group = -1;
    assert(groupRegsByApex(&singleton, &group) == 1 && group == 0,
           "A singleton must receive a valid group ID", NULL);
}

static void test_indexed_grouping_matches_pairwise_graph()
{
    // Independent all-pairs reachability oracle, including disjoint windows.
    for (size_t trial = 0; trial < 50; ++trial)
    {
        const size_t n = 24;
        std::vector<RegressionGauss> candidates(n);
        bool reachable[n][n] = {};
        for (size_t i = 0; i < n; ++i)
        {
            candidates[i].span = {int32_t((i * 17 + trial * 7) % 50), int32_t(5 + (i * 11 + trial) % 16)};
            candidates[i].position = float(candidates[i].span.startIdx) + float(candidates[i].span.length - 1) / 2;
        }
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                reachable[i][j] = i == j ||
                    (std::abs(candidates[i].position - candidates[j].position) < 2 * GLOBAL_MINSCALE - 1 + 2 * FLT_EPSILON &&
                     candidates[i].span.startIdx <= candidates[j].span.endIdx() &&
                     candidates[j].span.startIdx <= candidates[i].span.endIdx());
        for (size_t k = 0; k < n; ++k)
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    reachable[i][j] = reachable[i][j] || (reachable[i][k] && reachable[k][j]);
        int16_t groups[n];
        groupRegsByApex(&candidates, groups);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                assert((groups[i] == groups[j]) == reachable[i][j],
                       "Indexed grouping differs from pairwise graph", NULL);
    }
}

static void test_joint_selection_keeps_both_supported_maxima()
{
    // Reduced ten-point profile: the independently best fits have apices
    // 3.4154 and 6.3443, below the existing three-sample resolution limit.
    // The already accepted alternative at 2.6431 allows both maxima to remain.
    std::vector<float> y = {32892.6719F,63924.832F,86557.8672F,87380.4922F,
        76967.9922F,76014.4922F,84940.7188F,84153.8281F,63490.9297F,34147.457F};
    for (int reflected = 0; reflected < 2; ++reflected)
    {
        std::vector<float> x(y.size());
        for (size_t i = 0; i < x.size(); ++i) x[i] = float(i);
        std::vector<RegressionGauss> peaks;
        const auto count = qpeaks_find(y.data(), x.data(), nullptr, y.size(), &peaks);
        assert(count == 2, "Joint selection must retain both supported maxima", NULL);
        std::sort(peaks.begin(), peaks.end(), [](const auto &a,const auto &b){return a.position < b.position;});
        assert(std::abs(peaks[0].position - 3) < 0.5 && std::abs(peaks[1].position - 6) < 0.5,
               "Selected fits must remain near both observed maxima", NULL);
        assert(peaks[1].position - peaks[0].position > 3 + 2 * FLT_EPSILON,
               "Jointly chosen fits must satisfy the unchanged resolution limit", NULL);
        std::reverse(y.begin(), y.end());
    }
    // Two separated copies must be solved independently without losing a maximum.
    std::vector<float> x(10), logs(10);
    for (size_t i = 0; i < 10; ++i) { x[i] = float(i); logs[i] = std::log(y[i]); }
    std::vector<RegCoeffs> coefficients;
    findCoefficients(logs.data(), y.size(), &coefficients);
    std::vector<RegressionGauss> candidates;
    validateRegressions(y.data(), logs.data(), x.data(), nullptr, &coefficients, y.size(), &candidates);
    const size_t originalCount = candidates.size();
    for (size_t i = 0; i < originalCount; ++i)
    {
        auto copy = candidates[i];
        copy.position += 30; copy.coeffs.x0 += 30; copy.span.startIdx += 30;
        candidates.push_back(copy);
    }
    std::vector<float> combined(40, 5);
    std::copy(y.begin(), y.end(), combined.begin());
    std::copy(y.begin(), y.end(), combined.begin() + 30);
    std::vector<int16_t> groups(candidates.size());
    const size_t groupCount = groupRegsByApex(&candidates, groups.data());
    std::vector<bool> conflicting(groupCount, true);
    resolveApexGroups(&candidates, combined.data(), combined.size(), groups.data(), conflicting);
    size_t retained = 0;
    for (const auto &candidate : candidates) retained += candidate.isValid;
    assert(retained == 4, "Independent joint-selection components must retain all four maxima", NULL);

}

static void test_log_fit_acceptance_precedes_height_correction()
{
    // This residual is orthogonal to all four design columns on [-2,2].
    // Its exact log-domain RSS is 70*e^2. The nested linear RSS adds 14*a^2.
    // Hence F=a^2/(10*e^2)=202.906..., above F(2,1;0.95)=199.5.
    // Raw-height correction is not an OLS log intercept and must not inflate
    // that residual or make this statistically accepted peak disappear.
    const double a = 0.5, e = 0.0111;
    const double residual[] = {1,-4,6,-4,1};
    std::vector<float> x(5), y(5), logs(5);
    for (size_t i = 0; i < 5; ++i)
    {
        x[i] = float(i);
        y[i] = float(1000 * std::exp(-a * (double(i)-2)*(double(i)-2) + e*residual[i]));
        logs[i] = std::log(y[i]);
    }
    std::vector<RegCoeffs> coefficients;
    findCoefficients(logs.data(), 5, &coefficients);
    assert(coefficients.size() == 1, "Five points must yield one four-coefficient fit", NULL);
    RegressionGauss original;
    original.coeffs = coefficients.front(); original.span = {0,5};
    const double rss = calcRSS_log(&original, logs.data());
    assert(std::abs(rss - 70*e*e) < 1e-6, "Log residual must match its independent orthogonal construction", NULL);
    assert(f_testRegression(logs.data(), rss, original.span), "The unchanged F test must accept the original log fit", NULL);
    std::vector<float> prediction(5);
    correctB0(y.data(), original.span, prediction.data(), &original.coeffs);
    assert(calcRSS_log(&original, logs.data()) > rss,
           "Fixture must distinguish raw-height correction from OLS log residual", NULL);
    std::vector<RegressionGauss> peaks;
    assert(qpeaks_find(y.data(), x.data(), nullptr, 5, &peaks) == 1,
           "Height correction must not overturn acceptance of this log-domain fit", NULL);
}

static void test_f_threshold_respects_significance_level()
{
    // F(2,1) has survival probability 1/sqrt(1+2F), so the
    // upper-tail alpha quantile is exactly (1/alpha^2 - 1)/2.
    // Repeated 0.05 calls bracket other levels to exercise cache isolation.
    for (const double alpha : {0.05, 0.1, 0.01, 0.05})
    {
        const double expected = (1/(alpha*alpha)-1)/2;
        const double actual = F_stat(alpha, 4, 2, 5);
        assert(std::abs(actual-expected) < 1e-6,
               "F threshold must use the requested alpha: %.9g versus %.9g", actual, expected);
    }
}

static void test_joint_search_propagates_constraints_without_dropping_groups()
{
    const std::vector<float> y = {5,10,50,100,50,45,90,45,40,80,40,10,5};
    std::vector<RegressionGauss> candidates;
    for (size_t group = 0; group < 3; ++group)
        for (const double offset : {-0.5, 0.0, 0.5})
        {
            const double position = 3 + 3*double(group) + offset;
            RegressionGauss reg;
            reg.position = float(position);
            reg.span = {0,13}; reg.coeffs.scale = 6; reg.coeffs.x0 = 6;
            reg.coeffs.b0 = std::log(100-10*double(group)) - 0.2*(position-6)*(position-6);
            reg.coeffs.b1 = 0.4*(position-6);
            reg.coeffs.b2 = reg.coeffs.b3 = -0.2;
            candidates.push_back(reg);
        }
    const auto original = candidates;
    std::vector<int16_t> groups(candidates.size(), 0);
    const std::vector<bool> conflicting = {true};
    resolveApexGroups(&candidates, y.data(), y.size(), groups.data(), conflicting);
    std::vector<float> selected;
    for (const auto &reg : candidates) if (reg.isValid) selected.push_back(reg.position);
    std::sort(selected.begin(), selected.end());
    // Enumerating the 3^3 alternatives under strict >3 separation leaves
    // exactly this combination. Per-maximum greedy choices cannot guarantee it.
    assert(selected == std::vector<float>({2.5F,6,9.5F}),
           "Joint search must find the unique feasible three-group assignment", NULL);
    candidates = original;
    candidates.erase(candidates.begin()); // Remove the sole feasible left option.
    groups.assign(candidates.size(), 0);
    bool rejected = false;
    try { resolveApexGroups(&candidates, y.data(), y.size(), groups.data(), conflicting); }
    catch (const std::runtime_error &error)
    { rejected = std::string(error.what()).find("No compatible accepted fits") != std::string::npos; }
    assert(rejected, "An infeasible assignment must fail instead of discarding a maximum", NULL);

    // Duplicating alternatives does not change feasibility. A fixed left-to-
    // right search visits over 200,000 states before proving this impossible.
    // The single right option constrains all middle choices immediately.
    candidates.clear();
    for (const size_t i : {1,2,3,4,5})
        for (size_t copy = 0; copy < 500; ++copy) candidates.push_back(original[i]);
    candidates.push_back(original[8]);
    groups.assign(candidates.size(), 0);
    rejected = false;
    try { resolveApexGroups(&candidates, y.data(), y.size(), groups.data(), conflicting); }
    catch (const std::runtime_error &error)
    { rejected = std::string(error.what()).find("No compatible accepted fits") != std::string::npos; }
    assert(rejected, "Constraint propagation must prove infeasibility within the unchanged search budget", NULL);
}

static void test_exact_log_fit_allows_zero_residual()
{
    RegressionGauss reg;
    reg.coeffs.b0 = 8; reg.coeffs.b1 = 0; reg.coeffs.b2 = reg.coeffs.b3 = -1;
    reg.coeffs.x0 = 2; reg.coeffs.scale = 2; reg.span = {0,5};
    const float logPeak[] = {4,7,8,7,4};
    const float logLine[] = {4,5,6,7,8};
    assert(calcRSS_log(&reg, logPeak) == 0, "An exact log fit has zero residual", NULL);
    assert(f_testRegression(logPeak, 0, reg.span), "Exact peak must improve on a non-exact line", NULL);
    assert(!f_testRegression(logLine, 0, reg.span), "An exact line cannot justify extra parameters", NULL);
    assert(!insignificantCoefficient(&reg.coeffs, 0, 1), "Nonzero curvature with zero residual remains defined", NULL);
    reg.coeffs.b2 = reg.coeffs.b3 = 0;
    assert(insignificantCoefficient(&reg.coeffs, 0, 1), "Zero curvature does not define a peak", NULL);
}

static void test_three_mandatory_groups_inside_two_resolution_intervals()
{
    // Synthetic analogue of an observed EIC conflict; no private measurement
    // data is embedded here. Accepted alternatives are associated with three
    // observed maxima, but all apex positions span only 5.95 samples.
    const std::vector<float> y = {5,10,30,100,60,50,90,70,95,65,40,10,5};
    const std::vector<std::vector<float>> positions = {
        {2.60F,2.71F,3.58F,3.95F},
        {5.15F,5.19F,5.84F,5.97F,6.72F,6.91F},
        {7.53F,7.54F,7.82F,8.02F,8.06F,8.43F,8.55F}};
    const double resolution = 3 + 2*FLT_EPSILON;
    size_t feasible = 0;
    for (float a : positions[0]) for (float b : positions[1]) for (float c : positions[2])
        feasible += b-a > resolution && c-b > resolution;
    assert(feasible == 0 && positions[2].back()-positions[0].front() < 2*resolution,
           "Three mandatory overlapping descriptions cannot fit this apex range", NULL);
    std::vector<RegressionGauss> candidates;
    for (const auto &group : positions)
        for (float position : group)
        {
            RegressionGauss reg;
            reg.position = position; reg.span = {0,13};
            reg.coeffs.x0 = 6; reg.coeffs.scale = 6;
            reg.coeffs.b0 = std::log(100.0)-0.15*(position-6)*(position-6);
            reg.coeffs.b1 = 0.3*(position-6);
            reg.coeffs.b2 = reg.coeffs.b3 = -0.15;
            candidates.push_back(reg);
        }
    std::vector<int16_t> groups(candidates.size(), 0);
    const std::vector<bool> conflicting = {true};
    bool rejected = false;
    try { resolveApexGroups(&candidates, y.data(), y.size(), groups.data(), conflicting); }
    catch (const std::runtime_error &error)
    { rejected = std::string(error.what()).find("No compatible accepted fits") != std::string::npos; }
    assert(rejected, "Impossible mandatory groups must not succeed by silently omitting a maximum", NULL);
}

static void test_linear_null_model_residuals()
{
    // These residuals are orthogonal to both columns of the linear model:
    // sum(e_i) = sum(i*e_i) = 0; sum(e_i^2) = 14/16 = 0.875.
    const float residuals[] = {0.25F, -0.5F, 0.5F, -0.5F, 0.25F};
    for (const float slope : {2.0F, -2.0F})
    {
        std::vector<float> exact(5);
        std::vector<float> perturbed(5);
        for (size_t i = 0; i < exact.size(); ++i)
        {
            exact[i] = 20.0F + slope * (float)i;
            perturbed[i] = exact[i] + residuals[i];
        }
        const double residual_sum = calcRSS_H0_cf2(perturbed.data(), {0, 5});
        assert(std::abs(residual_sum - 0.875) < 1e-12,
               "Linear null-model RSS for slope %f must be 0.875, found %f\n",
               slope, residual_sum);
        const double exact_sum = calcRSS_H0_cf2(exact.data(), {0, 5});
        assert(std::abs(exact_sum) < 1e-12,
               "An exact line with slope %f must have zero residual sum, found %f\n",
               slope, exact_sum);
    }
}

struct SyntheticGaussian
{
    double apex_index;
    double height;
    double sigma;
};

static std::vector<float> synthetic_gaussians(size_t length,
                                              const std::vector<SyntheticGaussian> &components)
{
    std::vector<float> y(length);
    for (size_t i = 0; i < length; ++i)
    {
        double intensity = 0;
        for (const auto &component : components)
        {
            const double z = ((double)i - component.apex_index) / component.sigma;
            intensity += component.height * std::exp(-0.5 * z * z);
        }
        y[i] = (float)intensity;
    }
    return y;
}

static std::vector<RegressionGauss> check_synthetic_peak_selection(
    const char *name, const std::vector<float> &y, const std::vector<SyntheticGaussian> &expected)
{
    // All cases below have analytically distinct modes. Their count and the
    // half-sample position tolerance were fixed before running the algorithm.
    const double x_start = 100;
    const double x_step = 0.2;
    std::vector<float> x(y.size());
    for (size_t i = 0; i < x.size(); ++i)
        x[i] = 100.0F + 0.2F * (float)i;

    std::vector<RegressionGauss> peaks;
    int32_t count = -1;
    try
    {
        count = qpeaks_find(y.data(), x.data(), nullptr, y.size(), &peaks);
    }
    catch (const std::exception &error)
    {
        (void)fprintf(stderr, "Synthetic case %s threw: %s\n", name, error.what());
        assert(false, "Synthetic peak selection must not reject the known profile: %s\n", name);
    }
    assert(count == (int32_t)expected.size() && peaks.size() == expected.size(),
           "%s: expected %zu peaks, found %d\n", name, expected.size(), count);
    std::sort(peaks.begin(), peaks.end(), [](const RegressionGauss &left, const RegressionGauss &right)
              { return left.position < right.position; });
    for (size_t i = 0; i < peaks.size(); ++i)
    {
        const RegressionGauss &peak = peaks[i];
        const double apex = x_start + x_step * expected[i].apex_index;
        assert(std::isfinite(peak.position) && std::abs(peak.position - apex) <= 0.5 * x_step,
               "%s: expected an apex near %f, found %f\n", name, apex, peak.position);
        assert(std::isfinite(peak.height) && peak.height > 0 &&
                   std::isfinite(peak.area) && peak.area > 0,
               "%s: returned peak properties must be finite and positive\n", name);

        // Diagnostics only: finding two maxima does not establish unbiased
        // quantitative deconvolution of their overlapping component areas.
        const double area = expected[i].height * expected[i].sigma * x_step * std::sqrt(2 * std::acos(-1.0));
        printf("Synthetic %s peak %zu: height expected %.8f, found %.8f, relative error %.6f%%; "
               "area expected %.8f, found %.8f, relative error %.6f%%\n",
               name, i, expected[i].height, peak.height,
               100 * (peak.height / expected[i].height - 1),
               area, peak.area, 100 * (peak.area / area - 1));
        (void)fflush(stdout);
    }
    return peaks;
}

static void test_analytic_peak_conflicts()
{
    // The valleys between these Gaussian pairs stay above e. Both signals
    // therefore enter one continuum, unlike separately thresholded peaks.
    for (const int separation : {8, 6})
    {
        const std::vector<SyntheticGaussian> components = {
            {12, 1000, 1.2}, {12.0 + separation, 700, 1.2}};
        std::vector<float> y = synthetic_gaussians((size_t)(25 + separation), components);
        const char *name = separation == 8 ? "two_gaussians_separation_8" : "two_gaussians_separation_6";
        check_synthetic_peak_selection(name, y, components);
        if (separation == 6)
        {
            std::reverse(y.begin(), y.end());
            check_synthetic_peak_selection("two_gaussians_separation_6_mirrored", y,
                                           {{12, 700, 1.2}, {18, 1000, 1.2}});
        }
    }
    const std::vector<SyntheticGaussian> unequal = {{12.25, 400, 1.2}, {17.75, 1000, 1.2}};
    std::vector<float> close_pair = synthetic_gaussians(31, unequal);
    check_synthetic_peak_selection("two_close_unequal_gaussians", close_pair, unequal);
    std::reverse(close_pair.begin(), close_pair.end());
    check_synthetic_peak_selection("two_close_unequal_gaussians_mirrored", close_pair,
                                   {{12.25, 1000, 1.2}, {17.75, 400, 1.2}});

    // This broad, plainly bimodal profile reproduced the unresolved-group
    // exception in 1b992a5. It is generated analytically, not from a user file.
    const std::vector<SyntheticGaussian> broad_pair = {{14, 1000, 3}, {26, 700, 3}};
    check_synthetic_peak_selection("two_broad_gaussians", synthetic_gaussians(41, broad_pair), broad_pair);

    // Resolving an ambiguous group must not split one broad maximum into
    // several peaks or invent peaks on a strictly monotone positive trace.
    const std::vector<SyntheticGaussian> single = {{30.25, 1000, 6}};
    check_synthetic_peak_selection("broad_single_gaussian", synthetic_gaussians(61, single), single);
    std::vector<float> monotone(31);
    for (size_t i = 0; i < monotone.size(); ++i)
        monotone[i] = (float)(20.0 * std::exp(0.05 * (double)i));
    check_synthetic_peak_selection("monotone_increasing", monotone, {});
    std::reverse(monotone.begin(), monotone.end());
    check_synthetic_peak_selection("monotone_decreasing", monotone, {});
}

static bool same_peak_description(const RegressionGauss &a, const RegressionGauss &b)
{
    return a.coeffs.b0 == b.coeffs.b0 && a.coeffs.b1 == b.coeffs.b1 &&
           a.coeffs.b2 == b.coeffs.b2 && a.coeffs.b3 == b.coeffs.b3 &&
           a.coeffs.scale == b.coeffs.scale && a.coeffs.x0 == b.coeffs.x0 &&
           a.span.startIdx == b.span.startIdx && a.span.length == b.span.length &&
           a.position == b.position && a.height == b.height && a.area == b.area &&
           a.position_unc == b.position_unc && a.height_unc == b.height_unc &&
           a.area_unc == b.area_unc && a.fwhm == b.fwhm && a.jaccard == b.jaccard &&
           a.dqs == b.dqs && a.df == b.df && a.numCompetitors == b.numCompetitors &&
           a.apexGroup == b.apexGroup && a.isValid == b.isValid;
}

static void test_conflict_preserves_stable_neighbor()
{
    // Fixed before testing: the broad pair is joined by one >e continuum to a
    // smaller, independent peak. The selected fit of that stable group must
    // survive a repair elsewhere in the continuum without refitting.
    const std::vector<SyntheticGaussian> components = {{14, 1000, 3}, {26, 700, 3}, {40, 100, 2}};
    const std::vector<float> y = synthetic_gaussians(51, components);
    std::vector<float> x(y.size());
    for (size_t i = 0; i < x.size(); ++i)
        x[i] = 100.0F + 0.2F * (float)i;
    size_t start = 0;
    size_t end = y.size();
    while (start < end && y[start] <= minIntensity_global)
        ++start;
    while (end > start && y[end - 1] <= minIntensity_global)
        --end;
    const size_t length = end - start;
    std::vector<float> log_y(length);
    for (size_t i = 0; i < length; ++i)
    {
        assert(y[start + i] > minIntensity_global,
               "The three analytical peaks must share one continuum\n", NULL);
        log_y[i] = std::log(y[start + i]);
    }
    std::vector<RegCoeffs> coefficients;
    findCoefficients(log_y.data(), length, &coefficients);
    std::vector<RegressionGauss> candidates;
    validateRegressions(y.data() + start, log_y.data(), x.data() + start, nullptr,
                        &coefficients, length, &candidates);
    std::vector<int16_t> groups(candidates.size(), -1);
    const size_t group_count = groupRegsByApex(&candidates, groups.data());
    std::vector<bool> conflicting(group_count, false);
    size_t stable_neighbor = candidates.size();
    bool has_conflict = false;
    for (size_t group = 0; group < group_count; ++group)
    {
        if (groupApexIsStable(&candidates, groups.data(), (int16_t)group) != -1)
        {
            conflicting[group] = true;
            has_conflict = true;
            continue;
        }
        const size_t chosen = selectFromGroup(&candidates, y.data() + start, groups.data(), (int16_t)group);
        candidates[chosen].isValid = true;
        if (std::abs(candidates[chosen].position + (double)start - 40.0) <= 0.5)
            stable_neighbor = chosen;
    }
    assert(has_conflict && stable_neighbor < candidates.size(),
           "This fixture must contain both a conflict and a stable smaller neighbor\n", NULL);

    RegressionGauss reference = candidates[stable_neighbor];
    adjustRegression(&reference, x.data() + start);
    reference.span.startIdx += (int32_t)start;
    const auto peaks = check_synthetic_peak_selection("conflict_with_stable_neighbor", y, components);
    assert(same_peak_description(peaks.back(), reference),
           "Conflict resolution changed the previously selected stable neighbor\n", NULL);

    // Reflect the same candidates and fixed group membership to isolate the
    // repair's symmetry. Regrouping the reflected raw profile is a separate
    // operation: the pre-existing grouping heuristic is not reflection-invariant.
    std::vector<float> reflected_y(y.begin() + start, y.begin() + end);
    std::reverse(reflected_y.begin(), reflected_y.end());
    for (RegressionGauss &candidate : candidates)
    {
        candidate.position = (float)(length - 1) - candidate.position;
        candidate.span.startIdx = (int32_t)length - candidate.span.startIdx - candidate.span.length;
        candidate.coeffs.x0 = (uint16_t)(length - 1 - candidate.coeffs.x0);
        candidate.coeffs.b1 = -candidate.coeffs.b1;
        std::swap(candidate.coeffs.b2, candidate.coeffs.b3);
    }
    reference = candidates[stable_neighbor];
    resolveApexGroups(&candidates, reflected_y.data(), length, groups.data(), conflicting);
    assert(same_peak_description(candidates[stable_neighbor], reference),
           "Reflected conflict resolution changed the stable neighbor\n", NULL);
    std::vector<float> selected_positions;
    const double reflected_start = (double)(y.size() - end);
    for (const RegressionGauss &candidate : candidates)
        if (candidate.isValid)
            selected_positions.push_back(candidate.position + (float)reflected_start);
    std::sort(selected_positions.begin(), selected_positions.end());
    assert(selected_positions.size() == 3,
           "Reflected conflict resolution must select three maxima without duplicates\n", NULL);
    const double reflected_centers[] = {10, 24, 36};
    for (size_t i = 0; i < selected_positions.size(); ++i)
        assert(std::abs(selected_positions[i] - reflected_centers[i]) <= 0.5,
               "Reflected conflict resolution misplaced peak %zu: %f\n", i, selected_positions[i]);
}

static void test_flat_maximum_is_one_group()
{
    // [4,8] is one connected maximum plateau. Candidate apexes on its two
    // ends must not be interpreted as two observed maxima or need a valley.
    const float y[] = {5, 10, 25, 70, 100, 100, 100, 100, 100, 70, 25, 10, 5};
    std::vector<RegressionGauss> candidates(2);
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        const uint16_t apex = i == 0 ? 4 : 8;
        candidates[i].position = (float)apex;
        candidates[i].span = {2, 9};
        candidates[i].coeffs = {std::log(100.0), 0, -0.125, -0.125, 4, apex};
    }
    const auto originals = candidates;
    const int16_t groups[] = {0, 0};
    resolveApexGroups(&candidates, y, sizeof(y) / sizeof(y[0]), groups, {true});
    size_t selected = 0;
    for (size_t i = 0; i < candidates.size(); ++i)
    {
        if (!candidates[i].isValid)
            continue;
        ++selected;
        RegressionGauss reference = originals[i];
        reference.isValid = true;
        assert(same_peak_description(candidates[i], reference),
               "Plateau grouping must select an unchanged candidate\n", NULL);
    }
    assert(selected == 1, "The two plateau ends must yield one selected maximum, found %zu\n", selected);
}

static void test_close_maxima_share_one_selected_fit()
{
    // Both distances are within the existing minimum resolution. The selector
    // must retain one existing description of the close pair, not claim two
    // resolved peaks, while preserving the independent stable neighbor.
    for (const int separation : {2, 3})
    {
        std::vector<float> y = {5, 10, 25, 60, 90, 100, 50, 80, 70, 50, 20,
                                10, 8, 20, 35, 40, 35, 20, 8, 5, 4};
        if (separation == 3)
            y.insert(y.begin() + 6, 50);
        const uint16_t apexes[] = {5, (uint16_t)(5 + separation), (uint16_t)(13 + separation)};
        const float heights[] = {100, 80, 40};
        std::vector<RegressionGauss> candidates(3);
        for (size_t i = 0; i < candidates.size(); ++i)
        {
            const uint16_t scale = i == 2 ? 3 : 4;
            RegressionGauss &candidate = candidates[i];
            candidate.position = (float)apexes[i];
            candidate.span = {(int32_t)apexes[i] - scale, 2 * scale + 1};
            candidate.coeffs = {std::log((double)heights[i]), 0, -0.125, -0.125, scale, apexes[i]};
            candidate.height = heights[i];
            candidate.area = (float)(heights[i] * 2 * std::sqrt(2 * std::acos(-1.0)));
            candidate.position_unc = 0.1F;
            candidate.height_unc = 1;
            candidate.area_unc = 2;
            candidate.df = (uint16_t)(candidate.span.length - 4);
            candidate.isValid = i == 2;
        }
        const auto originals = candidates;
        const int16_t groups[] = {0, 0, 1};
        resolveApexGroups(&candidates, y.data(), y.size(), groups, {true, false});
        assert(candidates.size() == originals.size(),
               "Selecting unresolved maxima must not create a new fitted candidate\n", NULL);
        assert(same_peak_description(candidates[2], originals[2]),
               "A nearby unresolved pair changed the independent stable neighbor\n", NULL);
        size_t selected_close = 0;
        for (size_t i = 0; i < 2; ++i)
        {
            if (!candidates[i].isValid)
                continue;
            ++selected_close;
            RegressionGauss reference = originals[i];
            reference.isValid = true;
            assert(same_peak_description(candidates[i], reference),
                   "The chosen unresolved-peak description must retain its original fit\n", NULL);
        }
        assert(selected_close == 1,
               "Maxima %d samples apart must select one unchanged fit, found %zu\n",
               separation, selected_close);
    }
}

static void test_known_reflected_three_peak_grouping()
{
    // Full-pipeline regression: preserve the original three-peak count and
    // half-sample location expectations. Also callable independently.
    std::vector<float> y = synthetic_gaussians(51, {{14, 1000, 3}, {26, 700, 3}, {40, 100, 2}});
    std::reverse(y.begin(), y.end());
    check_synthetic_peak_selection("known_reflected_three_peak_grouping", y,
                                   {{10, 100, 2}, {24, 700, 3}, {36, 1000, 3}});
}

static void test_features_follow_eic_retention_times()
{
    EIC eic;
    eic.scanNumbers = {4, 5, 7, 8, 10};
    eic.RT = {20, 20.5F, 21, 21.5F, 22};
    eic.df = {1, 2, 3, 4, 5};
    eic.mz.assign(5, 100);
    eic.DQSC.assign(5, 0.9F);
    eic.ints_area.assign(5, 0);
    simulate_gauss(&eic.RT, 21.05, 1000, 0.4, &eic.ints_area);

    // Other scans lie between this ion's observations. Taking a contiguous
    // slice of the full file's RT array would assign the wrong times to it.
    const std::vector<float> all_scan_times = {0, 1, 2, 3, 20, 20.5F, 20.75F, 21, 21.5F, 21.75F, 22};
    const std::vector<EIC> eics = {eic};
    std::vector<FeaturePeak> features;
    const size_t count = findFeatures(&eics, &all_scan_times, &features);
    assert(count == 1, "Expected one feature in the synthetic ion trace\n", NULL);
    assert(std::abs(features.front().retentionTime - 21.05F) < 1e-3F,
           "Feature apex must follow the EIC's own retention times: %f\n", features.front().retentionTime);
    assert(features.front().lowerRT == 20 && features.front().upperRT == 22,
           "Feature bounds must refer to its measured EIC region\n", NULL);
}

static void test_centroid_id_above_uint16()
{
    RegressionGauss peak;
    peak.position = 100;
    const CentroidPeak centroid = regToCen(&peak, 70000U, 42);
    assert(centroid.ID == 70000U, "Centroid ID was truncated: %u\n", centroid.ID);
    assert(centroid.number_MS1 == 42, "Centroid scan index changed during conversion\n", NULL);
}

#if 0

static double peakVal_EMG(double x, double apex, double height, double sdev, double tau)
{
    // implements the exponentially modified gaussian fit as a generative model

    // naive model: https://doi.org/10.1021/ac970481d, Eq. B1
    // y(x) = A/t * exp(s^2 / 2t^2 + (x_0 - x) / s)    *    (1 - erf( (x - x_0) / (sqrt(2) * s) - s / (sqrt(2) * t) ) )
    // A = area, t = tau (tailing), s = sdev, x_0 = apex
    // problem: the first term gets very large while the second term gets very small

    // alternative equation: https://doi.org/10.1002/cem.1343 Eq. 6
    // z = 1/sqrt(2) * ( (x_0 - x) / s + s / t )
    // y(x) = h * exp( -(x_0 - x)^2 / (2 * s^2) ) * (s/t) * sqrt(pi/2) * exp(z^2) * erfc(z)

    double z = (1 / sqrt(2)) * ((apex - x) / sdev + sdev / tau);
    double a = height * exp(-(apex - x) * (apex - x) / (2 * sdev * sdev)) * (sdev / tau) * sqrt(M_PI_2);
    double y = a * exp(z * z) * erfc(z);
    if (!(y < INFINITY))
        return -1;
    return (float)y;
}

static void simulate_EMG(
    float x_start, float x_step,
    double apex, double height, double sdev, double tau,
    std::vector<float> *xvals, std::vector<float> *yvals)
{
    assert(yvals->empty(), "wrong usage of test function\n", NULL);
    assert(xvals->empty(), "wrong usage of test function\n", NULL);

    float x = x_start;
    bool peakDone = false;
    for (size_t i = 0; i < 1000000; i++)
    {
        float y = (float)peakVal_EMG(x, apex, height, sdev, tau);
        if (y > 0.01)
        {
            yvals->push_back(y);
            xvals->push_back(x);
            peakDone = true;
        }
        else if (peakDone)
        {
            break;
        }
        x += x_step;
    }
}

// void print_regFit(const RegCoeffs *coeff, const std::vector<float> *x, const float delta_x)
// {
//     // only prints data for now
//     printf("    Predicted values:\n");
//     float x0 = x->at(coeff->x0);
//     for (size_t i = 0; i < x->size(); i++)
//     {
//         float xval = (x->at(i) - x0) / delta_x;
//         double y = regExpAt(coeff, xval);
//         printf("%f,", y);
//     }
//     printf("\n");
// }

struct ErrorEMG // NOLINT
{
    float r_tau;
    float d_area_rel, d_area_abs;
    float r_sdev, d_fwhm_rel, d_fwhm_abs;
    float r_apex, d_apex_abs;
    float r_height, d_height_rel, d_height_abs;
    float dqs, jaccard;
    float uncert_diff_position, uncert_diff_area, uncert_diff_height;
    bool negativeB23 = true;
};

static void control_sim_EMG(float x_start, float x_step, ErrorEMG *in_out)
{
    // generate data using an exponentially modified gaussian on an equidistant x axis

    float apex = in_out->r_apex;
    float sdev = in_out->r_sdev;
    float height = in_out->r_height;
    float tau = in_out->r_tau;

    std::vector<float> xvals;
    std::vector<float> yvals;

    simulate_EMG(x_start, x_step, apex, height, sdev, tau, &xvals, &yvals);
    x_start = xvals.front();
    // float x_end = xvals.back();
    size_t length = xvals.size();
    in_out->jaccard = 0;

    if (length < 5)
    {
        in_out->dqs = -2;
        return;
    }

    // take n equidistant samples from the generated distribution, starting at y_min.
    size_t idx_start = 0;
    for (; idx_start < length; idx_start++)
    {
        // @todo
    }

    std::vector<RegressionGauss> ret;
    qpeaks_find(yvals.data(), xvals.data(), nullptr, length, &ret);

    // printf("    Observed values:\n");
    // for (size_t i = 0; i < length; i++)
    // {
    //     printf("%f,", yvals[i]);
    // }
    // printf("\n");

    // if (ret.size() == 1)
    // {
    //     print_regFit(&ret.front().coeffs, &xvals, x_step);
    // }
    // else
    // {
    //     assert(ret.empty(), "error in test design!\n");
    //     printf("    No prediction possible, regression unfit\n####\n");
    // }

    xvals.clear();
    yvals.clear();
    simulate_EMG(x_start, x_step / 500, apex, height, sdev, tau, &xvals, &yvals);
    float fwhm_e = (float)fwhm_empiric(&xvals, &yvals);
    float area_e = (float)area_empiric(&xvals, &yvals);
    float apex_e = (float)position_empiric(&xvals, &yvals);
    float height_e = *maxVal(yvals.data(), yvals.size());

    in_out->r_height = height_e;
    in_out->r_apex = apex_e;

    // assert(ret.size() != 0, "Peak not found\n");
    // this check is removed since some configurations of parameters result in no peaks
    if (ret.empty())
    {
        in_out->dqs = -1;
        return;
    }

    assert(ret.size() == 1, "Too many peaks found\n", NULL);

    RegressionGauss reg = ret.front();

    float apex_p = reg.position;
    float height_p = reg.height;
    float fwhm_p = reg.fwhm;
    // float area_p = reg.area;

    RegCoeffs c = reg.coeffs;

    in_out->negativeB23 = (c.b2 < 0) && (c.b3 < 0);

    float area_c = (float)peakArea(&c, x_step, 0, nullptr);

    // empiric estimation of peak parameters

    // assert(abs(apex_e - apex_p) < 0.01, "inaccurate position (%f vs. %f)\n", apex_e, apex_p);
    // assert(abs(height_e - height_p) < height_e * 0.01, "inaccurate height (%f vs. %f)\n", height_e, height_p);
    // assert(abs(fwhm_e - fwhm_p) < fwhm_e * 0.05, "inaccurate width (%f vs. %f)\n", fwhm_e, fwhm_p);
    // assert(abs(area_e - area_p) < reg.uncert_area, "inaccurate area (%f vs. %f)\n", area_e, area_p);

    in_out->d_area_rel = abs(area_e - area_c) / area_c;
    in_out->d_area_abs = abs(area_e - area_c);
    in_out->d_fwhm_rel = abs(fwhm_e - fwhm_p) / fwhm_p;
    in_out->d_fwhm_abs = abs(fwhm_e - fwhm_p);
    in_out->d_apex_abs = abs(apex_e - apex_p);
    in_out->d_height_rel = abs(height_e - height_p) / height;
    in_out->d_height_abs = abs(height_e - height_p);
    in_out->dqs = reg.dqs;
    in_out->jaccard = reg.jaccard;
}
#endif

static int simulate_profile(
    const RegCoeffs *coeff,
    std::vector<float> *simulated,
    std::vector<float> *simulated_log)
{
    verify(coeff->x0 > 1);

    double x = -double(coeff->x0);
    for (size_t i = 0; i < coeff->x0; i++)
    {
        simulated->at(i) = (float)exp(regAt(coeff, x));
        x += 1;
    }
    for (size_t i = coeff->x0; i < simulated->size(); i++)
    {
        simulated->at(i) = (float)exp(regAt(coeff, x));
        x += 1;
    }
    simulated_log->clear();
    for (size_t i = 0; i < simulated->size(); i++)
    {
        simulated_log->push_back(log(simulated->at(i)));
    };
    return 0;
}

static void test_singlePeak()
{
    RegCoeffs coeff;
    coeff.b0 = 12;
    coeff.b1 = 0.32;
    coeff.b2 = -0.13;
    coeff.b3 = -0.22;
    coeff.x0 = 5;
    std::vector<float> simulated;
    simulated.resize(11);
    std::vector<float> simulated_log;
    std::vector<float> x_axis = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    simulate_profile(&coeff, &simulated, &simulated_log);

    std::vector<RegressionGauss> validRegs;
    regression_on_continuum(simulated.data(),
                            x_axis.data(),
                            simulated_log.data(),
                            nullptr,
                            simulated.size(),
                            &validRegs);

    assert(validRegs.size() == 1, "incorrect number of regressions found\n", NULL);
    double diff_b0 = abs(coeff.b0 - validRegs.front().coeffs.b0);
    double diff_b1 = abs(coeff.b1 - validRegs.front().coeffs.b1);
    double diff_b2 = abs(coeff.b2 - validRegs.front().coeffs.b2);
    double diff_b3 = abs(coeff.b3 - validRegs.front().coeffs.b3);

    assert(diff_b0 < abs(coeff.b0) * 0.05, "> 5%% Error in b0 estimate\n", NULL);
    assert(diff_b1 < abs(coeff.b1) * 0.05, "> 5%% Error in b1 estimate\n", NULL);
    assert(diff_b2 < abs(coeff.b2) * 0.05, "> 5%% Error in b2 estimate\n", NULL);
    assert(diff_b3 < abs(coeff.b3) * 0.05, "> 5%% Error in b3 estimate\n", NULL);
}

static int simulate_stepwise(
    const RegCoeffs *coeff,
    std::vector<float> *xvec,
    std::vector<float> *simulated,
    float delta_x)
{
    verify(coeff->x0 > 1);

    delta_x /= 20;
    float x = -delta_x * (float)coeff->x0;
    for (size_t i = 0; i < coeff->x0; i++)
    {
        xvec->push_back(x);
        simulated->at(i) = (float)exp(regAt(coeff, x));
        x += delta_x;
    }
    for (size_t i = coeff->x0; i < simulated->size(); i++)
    {
        xvec->push_back(x);
        simulated->at(i) = (float)exp(regAt(coeff, x));
        x += delta_x;
    }
    return 0;
}

static void test_areaPrediction()
{
    RegCoeffs coeff;
    coeff.b0 = 8;
    coeff.b1 = 0.32;
    coeff.b2 = -0.13;
    coeff.b3 = -0.22;
    coeff.x0 = 500;
    coeff.scale = 6;

    std::vector<float> x;
    std::vector<float> y;
    y.resize(1000);

    float x_step = 1;
    simulate_stepwise(&coeff, &x, &y, x_step); // note: the sampling rate is not equal to the stride of x
    double area_e = area_empiric(&x, &y);
    double area_t = peakArea(&coeff, x_step, 0, nullptr);

    assert(abs(area_e - area_t) < 0.01, "Incorrect area calculation when b2 and b3 < 0\n", NULL); // expected area > 10^5

    // also check for positive cases
    coeff.b2 = 0.012;
    x.clear();
    simulate_stepwise(&coeff, &x, &y, x_step);
    // remove values of y before the valley point to have equivalent areas
    double valley = -coeff.b1 / (2 * coeff.b2);
    size_t idx = 0;
    while (x[idx] < valley)
    {
        y[idx] = 0;
        idx += 1;
    }
    area_e = area_empiric(&x, &y);
    area_t = peakArea(&coeff, x_step, 0, nullptr);
    // note: this uses a greater tolerance as the expected difference because when using 1000 steps,
    // since the valley point is not included exactly. Even if the calculation is slightly unstable,
    // this error will always be below the expected variance from instrument noise.
    assert(abs(area_e - area_t) < y[idx] * abs(valley - x[idx]), "Incorrect area calculation when b2 > 0 and b3 < 0\n", NULL);

    coeff.b1 = -0.32; // this is necessary because a positive b1 means b3 is the apex
    coeff.b2 = -0.34;
    coeff.b3 = 0.021;
    x.clear();
    simulate_stepwise(&coeff, &x, &y, x_step);
    valley = -coeff.b1 / (2 * coeff.b3);
    idx = 999;
    while (x[idx] > valley)
    {
        y[idx] = 0;
        idx -= 1;
    }
    area_e = area_empiric(&x, &y);
    area_t = peakArea(&coeff, x_step, 0, nullptr);

    assert(abs(area_e - area_t) < y[idx] * abs(valley - x[idx]), "Incorrect area calculation when b2 < 0 and b3 > 0\n", NULL);

    // run the uncertainty calculation at least once to check for errors
    double uncert = 0;
    double mse = 0.03;
    coeff.b3 = -0.21;
    x.clear();
    simulate_stepwise(&coeff, &x, &y, x_step);
    area_e = area_empiric(&x, &y);
    area_t = peakArea(&coeff, x_step, mse, &uncert);
    assert(area_t > uncert, "Too high uncertainty\n", NULL);
}

int main(int argc, const char *argv[])
{
    // initialise seed for noise generating funtions
    srand(1234); // NOLINT

    if (argc == 2 && std::strcmp(argv[1], "--known-reflection-failure") == 0)
    {
        test_known_reflected_three_peak_grouping();
        return 0;
    }

    if (argc == 2 && std::strcmp(argv[1], "--legacy-pt03") == 0)
        return test_qpeaks_find(&pt_03);
    if (argc == 2 && std::strcmp(argv[1], "--legacy-pt04") == 0)
        return test_qpeaks_find(&pt_04);

    test_apex_conversion_nonunit_spacing();
    test_weighted_mean_ignores_interpolated_values();
    test_peak_at_region_boundaries_and_cumulative_df();
    test_apex_grouping_more_than_128_candidates();
    test_apex_grouping_order_and_reflection();
    test_indexed_grouping_matches_pairwise_graph();
    test_joint_selection_keeps_both_supported_maxima();
    test_known_reflected_three_peak_grouping();
    test_features_follow_eic_retention_times();
    test_centroid_id_above_uint16();
    test_linear_null_model_residuals();
    test_log_fit_acceptance_precedes_height_correction();
    test_exact_log_fit_allows_zero_residual();
    test_f_threshold_respects_significance_level();
    test_joint_search_propagates_constraints_without_dropping_groups();
    test_three_mandatory_groups_inside_two_resolution_intervals();
    test_analytic_peak_conflicts();
    test_conflict_preserves_stable_neighbor();
    test_flat_maximum_is_one_group();
    test_close_maxima_share_one_selected_fit();
    if (argc == 2 && std::strcmp(argv[1], "--regressions-only") == 0)
        return 0;

    test_qpeaks_set();
    test_singlePeak();
    test_areaPrediction();
    control_sim_gauss();

    // survey_EMG();

    return 0;
}
