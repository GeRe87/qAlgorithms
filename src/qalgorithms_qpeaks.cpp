#include "qalgorithms_qpeaks.h"
#include "qalgorithms_datatypes.h"
#include "qalgorithms_logging.h"
#include "qalgorithms_read_file.h"
#include "qalgorithms_utils.h"

#include "libcerf_reduced.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#if defined(_WIN32) && !defined(_USE_MATH_DEFINES)
    #define _USE_MATH_DEFINES // relevant for windows to have math constants
#endif
#include <math.h>
#include <vector>

#include <algorithm> // sorting - consider removing dependence on the STL here

namespace qAlgorithms
{
    // --- internal function declarations --- //

    // copy relevant values from regression struct to centroid struct - not all fields supported
    static inline CentroidPeak regToCen(const RegressionGauss *reg, uint32_t id, uint16_t specNum);

    static size_t groupRegsByApex(const std::vector<RegressionGauss> *validRegressions, int16_t apexGroups[]);

    static int32_t groupApexIsStable(
        const std::vector<RegressionGauss> *validRegressions,
        const int16_t apexGroups[],
        const int16_t groupNum);

    static size_t selectFromGroup(
        const std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const int16_t apexGroups[],
        const int16_t groupNum);

    static void resolveApexGroups(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const size_t length,
        const int16_t apexGroups[],
        const std::vector<bool> &conflictingGroups);

    static int32_t regression_on_continuum(
        const float *intensities,
        const float *x_axis,
        const float *intensities_log,
        const uint16_t *const df,
        const size_t length,
        std::vector<RegressionGauss> *result);

    /// @brief adjust the height of a regression to better fit the exponential data
    /// @param intensities non-logarithmic intensity values the regression was fitted to
    /// @param r range of the regression
    /// @param predicted empty vector that the predicted values for intensity (AFTER correction) are written to
    /// @param coeff coefficients that should be updated
    /// @return used correction factor
    static double correctB0(const float *const intensities,
                            const Span_i32 r,
                            float *predicted,
                            RegCoeffs *coeff);

    /// @brief perform various statistical tests to see if a regression describes a valid peak
    /// @param degreesOfFreedom_cum cumulative degrees of freedom (only relevant for interpolated data)
    /// @param intensities measured intensities
    /// @param intensities_log log of measured intensities - must have same length as intensities
    /// @param mutateReg regression that should be mutated by this function
    /// @return 0 if the regression is valid, otherwise the filter step which kicked it out
    static invalid calcRegressionProperties(
        const float *intensities,
        const float *intensities_log,
        const float *x_axis,
        const float *predict,
        const size_t df_sum,
        const size_t length,
        RegressionGauss *mutateReg);

    static void adjustRegression(RegressionGauss *reg, const float *x_values);

    // take a jacobian matrix as input and return the transpose at scale
    static double matProductReg(const double J[4], const size_t scale);

    static bool insignificantCoefficient(const RegCoeffs *coeffs, const double mse, const size_t df_sum);

    // utility functions for calculating regression values
    static double regAt(const RegCoeffs *coeff, const double x);

    static double fullWidthHalfMax(const RegCoeffs *coeff, const double height, const double delta_x);

    // this one does not include b0
    static double regExp_fac(const RegCoeffs *coeff, const double x);

    static double peakPositionUncert(const RegCoeffs *c, const double mse);
    static double peakHeightUncert(const RegCoeffs *c, const double mse);

    static double peakArea(const RegCoeffs *c, const double delta_x, const double mse, double *uncert);

    static inline double peakPosition(const RegCoeffs *c);

    const double minIntensity_global = 2.7182818284590452; // exp(1) == e

    // -------------------------------------- //

    static void adjustRegression(RegressionGauss *reg, const float *x_values)
    {
        assert(reg->isValid);

        float position = reg->position;
        size_t leftOfApex = size_t(position);
        float delta_x = x_values[leftOfApex + 1] - x_values[leftOfApex];
        // position is determined relative to the point left of the apex
        float apexFraction = delta_x * (position - trunc(position));
        reg->position = x_values[leftOfApex] + apexFraction;
        reg->position_unc *= delta_x;

        reg->fwhm = (float)fullWidthHalfMax(&reg->coeffs, reg->height, delta_x);

        reg->dqs = erfc(reg->area_unc / reg->area);
    }

    int32_t qpeaks_find(
        // SOME_IMPLEMENTATION_OF_LINEAR_ALLOCATOR_HERE
        const float *intensity_base,
        const float *x_values,
        const uint16_t *DF_cum,
        const size_t length,
        std::vector<RegressionGauss> *result)
    {
        // control input for nullpointers, mismatching x and y, and fitting maxscale
        if (intensity_base == nullptr || x_values == nullptr || result == nullptr)
        {
            return -1;
        }
        if (length < MINLENGTH)
        {
            return -2;
        }

        // @todo the assumtion that all values of x are equidistant is conrolled. If gaps are
        // found, the missing values are interpolated assuming an exponential rate of change.
        // this should happen before calling this function (?)

        /*
        The fitting routine assumes that all present peaks have a modified gaussian base function.
        This means that no baseline exists. Baseline substraction, if appropriate, has to be performed
        before calling the qpeaks_find function.
        */
        std::vector<float> intensity_base_log(length + 1);

        for (size_t i = 0; i < length; i++)
        {
            intensity_base_log[i] = log(intensity_base[i]);
        }

        // The minimum intensity is "> e ", because this way we know that all log intensities
        // treated in the program will be > 1. This reduces the risk of high floating point
        // errors when multiplying a number << 1 with one >> 1.

        // Process each complete contiguous region above the intensity threshold,
        // including a region that reaches the final sample.
        size_t rangeStart = 0;
        while (rangeStart < length)
        {
            while (rangeStart < length && !(intensity_base[rangeStart] > minIntensity_global))
                ++rangeStart;
            size_t pointIdx = rangeStart;
            while (pointIdx < length && intensity_base[pointIdx] > minIntensity_global)
                pointIdx += 1;
            const size_t newLen = pointIdx - rangeStart;
            if (newLen >= MINLENGTH)
            {
                if (newLen > UINT16_MAX)
                    throw std::runtime_error("A contiguous profile region exceeds the supported 65535-point limit.");
                const float *intensities = intensity_base + rangeStart;
                const float *intensities_log = intensity_base_log.data() + rangeStart;
                const float *x_axis = x_values + rangeStart;
                std::vector<uint16_t> localDF;
                if (DF_cum != nullptr)
                {
                    const uint16_t previous = rangeStart == 0 ? 0 : DF_cum[rangeStart - 1];
                    localDF.resize(newLen);
                    for (size_t i = 0; i < newLen; ++i)
                        localDF[i] = DF_cum[rangeStart + i] - previous;
                }
                const uint16_t *df = DF_cum == nullptr ? nullptr : localDF.data();

                assert(intensities_log[newLen - 1] > 1);

                const size_t previousResults = result->size();
                const int32_t status = regression_on_continuum(intensities,
                                        x_axis,
                                        intensities_log,
                                        df,
                                        newLen,
                                        result);
                if (status < 0)
                    return status;
                for (size_t i = previousResults; i < result->size(); ++i)
                    result->at(i).span.startIdx += (int32_t)rangeStart;
            }
            rangeStart = pointIdx + 1;
        }

        return (int32_t)result->size();
    }

#pragma region "running regression"

    static Span_i32 calcRegSpan(const RegCoeffs *coeffs);

    static size_t validateRegressions(
        const float *intensities,
        const float *intensities_log,
        const float *x_axis,
        const uint16_t *const degreesOfFreedom_cum,
        const std::vector<RegCoeffs> *coefficients,
        const size_t length,
        std::vector<RegressionGauss> *validRegressions)
    {
        size_t validCount = 0;
        std::vector<float> predict(length, 0);
        for (size_t i = 0; i < coefficients->size(); i++)
        {
            const RegCoeffs *coeffs = coefficients->data() + i;
            // sets the range and checks for validity
            Span_i32 span = calcRegSpan(coeffs);
            assert(span.startIdx >= 0);
            bool span_fail = (span.length < 5) || // this is redundant but kept for clarity
                             (span.startIdx + 2 > (int32_t)coeffs->x0) ||
                             (span.endIdx() - 2 < (int32_t)coeffs->x0);
            if (span_fail)
                continue;

            size_t df_sum = sumOfCumulative(degreesOfFreedom_cum, span.startIdx, span.length);
            if (df_sum < MINLENGTH)
                continue;

            df_sum -= 4; // four coefficients, adjust for components

            RegressionGauss reg;
            reg.coeffs = *coeffs;
            reg.span = span;

            /*
                Adjustment of b0 coefficient:
                When working with log-transformed data, the coefficients are suboptimal for the exponential case.
                Since we must work with a log system to perform a linear regression, there is a bias in the
                results which is somewhat corrected here. While correction could occur before validation, the
                initial tests filter out a lot of bad regressions which reduces processing time. The tests are
                presumed to be better when using the transformed coefficients in terms of applicability of the results.
                This function also modifies the "predict" vector supplied as its argument!
            */
            correctB0(intensities, span, predict.data(), &reg.coeffs);

            invalid failpoint = calcRegressionProperties(
                intensities,
                intensities_log,
                x_axis,
                predict.data(),
                df_sum,
                length,
                &reg);

            validCount += failpoint == ok ? 1 : 0;
            if (failpoint == ok)
            {
                validRegressions->push_back(reg);
            }
        }
        assert(validCount == validRegressions->size());

        return validCount;
    }

#pragma endregion "running regression"

#pragma region "Conflict Elimination"

    static int32_t regression_on_continuum(
        const float *intensities,
        const float *x_axis,
        const float *intensities_log,
        const uint16_t *const df,
        const size_t length,
        std::vector<RegressionGauss> *result)
    {

        // coefficients for single-b0 peaks, spans all regressions over a peak window
        // all entries in coeff are sorted by scale and position in ascending order - this is not checked!
        std::vector<RegCoeffs> coefficients;
        findCoefficients(intensities_log, length, &coefficients);

        std::vector<RegressionGauss> validRegressions; // all independently valid regressions regressions
        validRegressions.reserve(coefficients.size() / 2);
        size_t validCount = validateRegressions(intensities,
                                                intensities_log,
                                                x_axis,
                                                df,
                                                &coefficients,
                                                length,
                                                &validRegressions);
        if (validCount == 0)
            return 0;
        if (validCount == 1)
        {
            RegressionGauss *reg = validRegressions.data();

            // When the data is very large, a singular regression should be viewed with suspicion.
            // The only exception is made for data which is fully covered by the
            // regression (up to +1), since tests on expansion only make sense when this is possible
            bool fullCoverage = length - reg->span.length <= 1;
            if (!fullCoverage)
            {
                int32_t s = reg->span.startIdx;
                int32_t e = reg->span.endIdx();
                float intenstiy_L = s <= 0 ? 0 : intensities[s - 1];
                float intenstiy_R = size_t(e) + 1 >= length ? 0 : intensities[e + 1];
                fullCoverage = reg->height > max(intenstiy_L, intenstiy_R);
            }
            if (fullCoverage)
            {
                reg->isValid = true;
                adjustRegression(reg, x_axis);
                result->push_back(*reg);
            }
            // if there is no full coverage, we must assume that the peak is just noise. There
            // should probably be a more elaborated test than this, especially regarding false
            // detections of ill-formed peaks. False discarding of peaks with a positive coefficient
            // should not be a problem, since if those are formed but no closed form peak, something
            // is probably very wrong.

            return 0;
        }

        // rework of the apex selection, currently without the ability to interpolate points
        std::vector<int16_t> apexGroups(validRegressions.size(), -1);
        size_t apexcount = groupRegsByApex(&validRegressions, apexGroups.data());
        assert(apexcount < INT16_MAX);

        std::vector<bool> conflictingGroups(apexcount, false);
        bool hasConflict = false;
        for (int16_t groupNum = 0; groupNum < (int16_t)apexcount; groupNum++)
        {
            int32_t stableApex = groupApexIsStable(&validRegressions, apexGroups.data(), groupNum);
            if (stableApex == -1)
            {
                // equivalent to only one possible apex being in the data, no complex deconvolution is required
                size_t chosenOne = selectFromGroup(&validRegressions, intensities, apexGroups.data(), groupNum);
                validRegressions[chosenOne].isValid = true;
            }
            else
            {
                conflictingGroups[(size_t)groupNum] = true;
                hasConflict = true;
            }
        }
        if (hasConflict)
            resolveApexGroups(&validRegressions, intensities, length,
                             apexGroups.data(), conflictingGroups);

        for (size_t peak = 0; peak < validRegressions.size(); peak++)
        {
            RegressionGauss *reg = validRegressions.data() + peak;
            if (reg->isValid)
            {
                adjustRegression(reg, x_axis);
                result->push_back(*reg);
            }
        }
        return 0;
    }

    // Refine only conflicting groups. A broad fit whose apex lies in the valley
    // between two peaks can join their candidates transitively, even though the
    // measured maxima and smaller valid fits distinguish the peaks. Such a fit
    // must not become a third peak or suppress a supported neighbouring maximum.
    //
    // Coefficients, fit spans, degrees of freedom and uncertainties are retained.
    // Valleys bound only the common SSE comparison interval; they do not define
    // a new regression window. This selects peak descriptions, not quantitative
    // deconvolution of the overlapping signals.
    static void resolveApexGroups(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const size_t length,
        const int16_t apexGroups[],
        const std::vector<bool> &conflictingGroups)
    {
        struct ObservedMaximum
        {
            size_t first;
            size_t last;
        };
        const size_t count = validRegressions->size();
        std::vector<int32_t> support(count, -1);
        std::vector<ObservedMaximum> maxima;
        std::vector<bool> supportedGroups(conflictingGroups.size(), false);

        for (size_t i = 0; i < count; ++i)
        {
            const RegressionGauss &reg = validRegressions->at(i);
            const size_t group = (size_t)apexGroups[i];
            // Already selected stable representatives remain eligible unchanged.
            // Other candidates from a stable group are outside this refinement.
            if (!conflictingGroups[group] && !reg.isValid)
                continue;

            const size_t left = (size_t)std::floor(reg.position);
            const size_t right = (size_t)std::ceil(reg.position);
            for (size_t p = left; p <= right && p < length; ++p)
            {
                size_t first = p;
                size_t last = p;
                while (first > 0 && intensities[first - 1] == intensities[p])
                    --first;
                while (last + 1 < length && intensities[last + 1] == intensities[p])
                    ++last;
                // Treat a flat maximum as one observed apex. Both outside
                // neighbours must be lower, so a flat baseline is not a peak.
                if (first == 0 || last + 1 >= length ||
                    intensities[first - 1] >= intensities[p] ||
                    intensities[last + 1] >= intensities[p])
                    continue;
                support[i] = (int32_t)first;
                maxima.push_back({first, last});
                supportedGroups[group] = true;
                break;
            }
        }
        for (size_t group = 0; group < conflictingGroups.size(); ++group)
        {
            if (conflictingGroups[group] && !supportedGroups[group])
                throw std::runtime_error("An overlapping apex group has no observed maximum supported by a valid fit.");
        }

        std::sort(maxima.begin(), maxima.end(), [](const ObservedMaximum &a, const ObservedMaximum &b)
                  { return a.first < b.first; });
        maxima.erase(std::unique(maxima.begin(), maxima.end(), [](const ObservedMaximum &a, const ObservedMaximum &b)
                                 { return a.first == b.first; }), maxima.end());
        assert(!maxima.empty());

        const double resolution = 2 * GLOBAL_MINSCALE - 1 + 2 * FLT_EPSILON;
        std::vector<bool> initiallySelected(count, false);
        for (size_t i = 0; i < count; ++i)
            initiallySelected[i] = validRegressions->at(i).isValid;

        // Nearby measured maxima need not represent separable fitted peaks.
        // Combine their alternatives only when the entire continuous apex
        // range remains inside the existing resolution limit and the fit
        // spans overlap. Testing the full range prevents transitive chains
        // from swallowing alternatives farther apart than that limit.
        while (maxima.size() > 1)
        {
            size_t mergeAt = maxima.size();
            double smallestRange = INFINITY;
            bool ambiguous = false;
            for (size_t k = 0; k + 1 < maxima.size(); ++k)
            {
                double lowerApex = INFINITY;
                double upperApex = -INFINITY;
                bool stableLeft = false;
                bool stableRight = false;
                int32_t leftStart = INT32_MAX;
                int32_t leftEnd = 0;
                int32_t rightStart = INT32_MAX;
                int32_t rightEnd = 0;
                for (size_t i = 0; i < count; ++i)
                {
                    const bool left = support[i] == (int32_t)maxima[k].first;
                    const bool right = support[i] == (int32_t)maxima[k + 1].first;
                    if (!left && !right)
                        continue;
                    const RegressionGauss &reg = validRegressions->at(i);
                    lowerApex = std::min(lowerApex, (double)reg.position);
                    upperApex = std::max(upperApex, (double)reg.position);
                    if (left)
                    {
                        stableLeft = stableLeft || reg.isValid;
                        leftStart = std::min(leftStart, reg.span.startIdx);
                        leftEnd = std::max(leftEnd, reg.span.endIdx());
                    }
                    if (right)
                    {
                        stableRight = stableRight || reg.isValid;
                        rightStart = std::min(rightStart, reg.span.startIdx);
                        rightEnd = std::max(rightEnd, reg.span.endIdx());
                    }
                }
                if ((stableLeft && stableRight) || leftEnd < rightStart || rightEnd < leftStart)
                    continue;
                const double apexRange = upperApex - lowerApex;
                if (apexRange > resolution)
                    continue;
                if (apexRange < smallestRange)
                {
                    smallestRange = apexRange;
                    mergeAt = k;
                    ambiguous = false;
                }
                else if (apexRange == smallestRange && k == mergeAt + 1)
                    ambiguous = true;
            }
            // Equally close alternatives sharing a group must not be resolved
            // just by their processing order. The final check retains an
            // explicit error if their selected descriptions remain too close.
            if (mergeAt == maxima.size() || ambiguous)
                break;
            const int32_t previousKey = (int32_t)maxima[mergeAt + 1].first;
            const int32_t combinedKey = (int32_t)maxima[mergeAt].first;
            for (int32_t &key : support)
                if (key == previousKey)
                    key = combinedKey;
            maxima[mergeAt].last = maxima[mergeAt + 1].last;
            maxima.erase(maxima.begin() + (ptrdiff_t)mergeAt + 1);
        }

        // Keep the complete minimum plateau in both adjacent scoring domains.
        // This avoids an arbitrary left/right choice when intensities are tied.
        std::vector<ObservedMaximum> valleys;
        for (size_t k = 1; k < maxima.size(); ++k)
        {
            const ObservedMaximum &left = maxima[k - 1];
            const ObservedMaximum &right = maxima[k];
            ObservedMaximum valley{left.last + 1, left.last + 1};
            for (size_t p = left.last + 2; p < right.first; ++p)
            {
                if (intensities[p] < intensities[valley.first])
                    valley = {p, p};
                else if (intensities[p] == intensities[valley.first])
                    valley.last = p;
            }
            if (!(intensities[valley.first] < std::min(intensities[left.first], intensities[right.first])))
                throw std::runtime_error("Supported neighbouring maxima have no separating observed valley.");
            valleys.push_back(valley);
        }

        for (size_t k = 0; k < maxima.size(); ++k)
        {
            std::vector<size_t> candidates;
            std::vector<size_t> selectedStable;
            for (size_t i = 0; i < count; ++i)
            {
                if (support[i] != (int32_t)maxima[k].first)
                    continue;
                if (validRegressions->at(i).isValid)
                    selectedStable.push_back(i);
                else if (conflictingGroups[(size_t)apexGroups[i]])
                    candidates.push_back(i);
            }
            // Preserve a previously selected stable description exactly. If
            // two stable groups describe the same observed maximum, compare
            // only those representatives and emit that maximum once.
            if (selectedStable.size() == 1)
                continue;
            if (!selectedStable.empty())
                candidates = selectedStable;
            assert(!candidates.empty());

            size_t limitLeft = length;
            size_t limitRight = 0;
            for (size_t i : candidates)
            {
                const Span_i32 span = validRegressions->at(i).span;
                limitLeft = std::min(limitLeft, (size_t)span.startIdx);
                limitRight = std::max(limitRight, (size_t)span.endIdx());
            }
            if (k > 0)
                limitLeft = std::max(limitLeft, valleys[k - 1].first);
            if (k < valleys.size())
                limitRight = std::min(limitRight, valleys[k].last);
            assert(limitLeft <= limitRight);

            double bestSSE = INFINITY;
            size_t best = count;
            for (size_t i : candidates)
            {
                const RegCoeffs &coeffs = validRegressions->at(i).coeffs;
                double sse = 0;
                for (size_t p = limitLeft; p <= limitRight; ++p)
                {
                    const double prediction = std::exp(regAt(&coeffs, (double)p - (double)coeffs.x0));
                    const double difference = (double)intensities[p] - prediction;
                    sse += difference * difference;
                    if (sse > bestSSE)
                        break;
                }
                if (sse < bestSSE)
                {
                    bestSSE = sse;
                    best = i;
                }
            }
            if (best == count)
                throw std::runtime_error("No finite fit error is available for a supported apex group.");
            for (size_t i : selectedStable)
                validRegressions->at(i).isValid = false;
            validRegressions->at(best).isValid = true;
        }

        // Check the actual chosen fits, not rounded observation indices or
        // former group seeds. Independent stable groups are preserved; every
        // newly selected overlapping description must satisfy the old limit.
        for (size_t i = 0; i < count; ++i)
        {
            if (!validRegressions->at(i).isValid)
                continue;
            for (size_t j = i + 1; j < count; ++j)
            {
                if (!validRegressions->at(j).isValid || (initiallySelected[i] && initiallySelected[j]))
                    continue;
                const RegressionGauss &left = validRegressions->at(i);
                const RegressionGauss &right = validRegressions->at(j);
                if (left.span.endIdx() < right.span.startIdx || right.span.endIdx() < left.span.startIdx)
                    continue;
                if (std::abs(left.position - right.position) <= resolution)
                    throw std::runtime_error("Selected overlapping fit apices remain below the existing resolution limit.");
            }
        }
    }

    static size_t groupRegsByApex(const std::vector<RegressionGauss> *validRegressions, int16_t apexGroups[])
    {
        const size_t length = validRegressions->size();
        if (length == 0)
            return 0;
        std::fill(apexGroups, apexGroups + length, int16_t(-1));
        const double resolution = 2 * GLOBAL_MINSCALE - 1 + 2 * FLT_EPSILON;
        size_t groupCount = 0;
        std::vector<size_t> pending;
        for (size_t seed = 0; seed < length; ++seed)
        {
            if (apexGroups[seed] != -1)
                continue;
            if (groupCount >= INT16_MAX)
                throw std::runtime_error("Too many apex groups in a single profile region.");
            const int16_t group = (int16_t)groupCount++;
            apexGroups[seed] = group;
            pending.clear();
            pending.push_back(seed);
            // Connected components of pairwise proximity AND window overlap.
            // Testing against a fixed seed window made membership depend on
            // candidate order; expanded apex bounds could also bridge a gap
            // using proximity and overlap from different candidates.
            for (size_t next = 0; next < pending.size(); ++next)
            {
                const RegressionGauss &member = validRegressions->at(pending[next]);
                for (size_t i = 0; i < length; ++i)
                {
                    if (apexGroups[i] != -1)
                        continue;
                    const RegressionGauss &candidate = validRegressions->at(i);
                    const bool overlap =
                        candidate.span.startIdx < member.span.startIdx + member.span.length &&
                        member.span.startIdx < candidate.span.startIdx + candidate.span.length;
                    if (overlap && std::abs(candidate.position - member.position) < resolution)
                    {
                        apexGroups[i] = group;
                        pending.push_back(i);
                    }
                }
            }
        }
        return groupCount;
    }

    // function that asserts that all regressions within one grouping describe only one apex
    static int32_t groupApexIsStable(
        const std::vector<RegressionGauss> *validRegressions,
        const int16_t apexGroups[],
        const int16_t groupNum)
    {
        const int32_t no_difference = -1;
        const size_t regCount = validRegressions->size();
        assert(regCount > 1);

        double apex_lim_L = INFINITY;
        double apex_lim_R = 0;

        bool differenceCandidate = false;

        for (size_t i = 0; i < regCount; i++)
        {
            if (apexGroups[i] != groupNum)
                continue;
            const RegressionGauss *reg = validRegressions->data() + i;
            apex_lim_L = min(apex_lim_L, reg->position);
            apex_lim_R = max(apex_lim_R, reg->position);

            // To separate two apexes unabiguously, there must be at least three points
            // between them. If this condition is not fulfilled for a group, the individual
            // apexes cannot be separated even if more than one group exists as a ground truth.
            differenceCandidate = apex_lim_R - apex_lim_L > 3 + 2 * FLT_EPSILON;
            if (differenceCandidate)
                break;
        }
        assert(apex_lim_L <= apex_lim_R);

        if (!differenceCandidate)
            return no_difference;

        // @todo a difference is very likely if there is a regression with positive
        // coeff in the data. In such a case, all regressions that are primarily (what does that mean?)
        // outside of that region must be reconsidered. Observed problem: decent regressions in
        // group, but two problematic ones with an apex beyond the first limit

        // also todo: eliminations based on statitics (currently only f test) should be done
        // only after the opinionated grouping step

        // even if more than one apex could exist in the given data, the assignment of different
        // groups only makes sense if the regressions in question do not have total overlap.
        // The reasoning here is that one regression was fit over the partial range of another
        // regression, the discrepancy between the two indicates a suboptimal fit resulting from not
        // observing the entire relevant data. One blindspot of a check just testing the two
        // outermost apexes is that (hypothetically) another pair violating this condition could
        // exist and describe two separate regressions, even if the checked pair is not sufficiently
        // different to demonstrate this. A potential cause for such a situation is present if a
        // small apex to the side of a larger one is clearly defined, but the very large regression
        // window next to it fully covers the fit region. Here, it helps that we know that all
        // regressions are sorted by scale when being processed: If we stop comparing after finding
        // the first potential conflict, by definition we are comparing the two smallest possible
        // cases of describing distinct peaks with each other, so inaccurate full overlap is unlikely.

        // Note on resolving the conflict:
        // If we determine that the group contains at least two more accurate groups, two things have
        // to be observed. Firstly, all present apexes must be assigned to one or the other group -
        // it should be reasonable to weight them the same and assign based on smallest absolute
        // distance. Secondly, determine if a regression range needs to be adjusted in order to
        // prevent selecting an optimal group representative based on profiles containing two apexes.
        // This should be fulfilled by truncating all regions of definition so they at most extend
        // to the point that leaves two observations to the apex.

        size_t bound_reg_L_L = 0;
        size_t bound_reg_L_R = 0;
        size_t bound_reg_R_L = 0;
        size_t bound_reg_R_R = 0;
        // Iterate through all regressions until the two outermost ones are found
        // always iterating everything is inefficient, but this part of the function will not run often
        for (size_t i = 0; i < regCount; i++)
        {
            if (apexGroups[i] != groupNum)
                continue;
            const RegressionGauss *reg = validRegressions->data() + i;
            if (reg->position == apex_lim_L)
            {
                bound_reg_L_L = reg->span.startIdx;
                bound_reg_L_R = reg->span.endIdx();
            }
            if (reg->position == apex_lim_R)
            {
                bound_reg_R_L = reg->span.startIdx;
                bound_reg_R_R = reg->span.endIdx();
            }
        }

        // check that either regression is fully contained within another, meaning its bounds
        // are always closer to the center than those of the competing regression. If either is
        // the case, the two apexes in question cannot be separated and the group apex is stable.
        bool left_in_right = bound_reg_L_L >= bound_reg_R_L && bound_reg_L_R <= bound_reg_R_R;
        bool right_in_left = bound_reg_L_L <= bound_reg_R_L && bound_reg_L_R >= bound_reg_R_R;
        bool stableApex = left_in_right || right_in_left;

        if (stableApex)
            return no_difference;

        // next step: refine the grouping by deciding which of the two conflicting apexes should
        // be moved into its own group on basis of distance to the mean / median apex

        return 0; // this will be the ID of the peak that violates the total group apex estimate the most @todo
    }

    // return the index of the regression that was found to be the best group representative
    static size_t selectFromGroup(
        const std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const int16_t apexGroups[],
        const int16_t groupNum)
    {
        assert(groupNum > -1);
        size_t regCount = validRegressions->size();
        assert(regCount > 1);

        // 1) iterate over the groups to find the relevant region over which to compare
        // the regressions
        int32_t lim_L = INT32_MAX;
        int32_t lim_R = 0;
        for (size_t i = 0; i < regCount; i++)
        {
            if (apexGroups[i] == groupNum)
            {
                const RegressionGauss *reg = validRegressions->data() + i;
                lim_L = min(reg->span.startIdx, lim_L);
                lim_R = max(lim_R, reg->span.endIdx());
            }
        }
        assert(lim_R > 0);

        // calculate SSE in order
        // We use the SSE instead of the MSE since n is identical for this
        // comparison
        double bestSSE = INFINITY;
        size_t bestRegIdx = regCount;

        for (size_t i = 0; i < regCount; i++)
        {
            if (apexGroups[i] != groupNum)
                continue;

            const RegressionGauss *reg = validRegressions->data() + i;

            double currentSSE = 0;
            for (size_t xval = lim_L; xval < (size_t)lim_R + 1; xval++)
            {
                double observed = intensities[xval];
                double x = double(xval) - (double)reg->coeffs.x0;
                double predict = exp(regAt(&reg->coeffs, x));
                double diff = observed - predict;
                currentSSE += diff * diff;

                // early termination of mse calculation should save some time on average, not measured!
                if (currentSSE > bestSSE)
                    break;
            }
            if (currentSSE < bestSSE)
            {
                bestSSE = currentSSE;
                bestRegIdx = i;
            }
        }
        if (bestRegIdx == regCount)
            throw std::runtime_error("No finite fit error is available for an apex group.");
        return bestRegIdx;
    }

#pragma endregion "eliminate conflicting regs"

    void findCoefficients(
        const float *intensity_log,
        const size_t length,
        std::vector<RegCoeffs> *coeffs)
    {
        const size_t maxscale = min(maxscale_global, (length - 1) / 2);
        assert(maxscale >= GLOBAL_MINSCALE);
        assert(maxscale <= QALGORITHMS_MAXSCALE_PRECOMPILED);

        assert(length > 4);

        // for every checked scale, there are length - 2 scale regressions performed
        // totalRegs = sum_i from minscale to maxscale (length - 2 * scale_i)
        // scalediff = maxscale - minscale + 1
        // totalRegs = scalediff * length - 2 * (2 + 3 + 4 + ... + maxscale)
        // the last bit is the (triangle number of maxscale) - (miscale - 1)
        // totalRegs = scalediff * length - 2 * (maxscale * (maxscale + 1) / 2 - minscale + 1)
        // totalRegs = scalediff * length - maxscale * (maxscale + 1) + 2 * minscale - 2)
        const size_t totalRegs = (maxscale - GLOBAL_MINSCALE + 1) * length -
                                 maxscale * (maxscale + 1) +
                                 (size_t)2 * GLOBAL_MINSCALE - 2; // typecast here to silence a clang-tidy warning
        coeffs->resize(totalRegs);

        // the first and last MINSCALE elements of the data do not need to be checked for x0, since they are invalid by definition
        const size_t limit = length - GLOBAL_MINSCALE;
        for (size_t x0 = GLOBAL_MINSCALE; x0 < limit; x0++)
        {
            const float *cen = intensity_log + x0; // this is initially the third real point

            // calculate the convolution with the kernel of the lowest scale - 1 (= 1), i.e. xT * intensity_log[cen - 1 : cen + 1]
            // the product sums are the rows of the design matrix (xT) * intensity_log[i:i+4] (dot product)
            // they are set to scale = 1 so the first value written is at scale = 2
            // b0 is 1, 1, 1,
            double product_sum_b0 = cen[-1] + cen[0] + cen[1];
            // b1 is -1, 0, 1
            double product_sum_b1 = -cen[-1] + cen[1];
            // b2 is 1, 0, 0
            double product_sum_b2 = cen[-1];
            // b3 is 0, 0, 1
            double product_sum_b3 = cen[1];

            size_t maxscale_absolute = 0;
            { // the largest valid scale depends on x0
                size_t maxscale_left = x0;
                size_t maxscale_right = length - x0 - 1;
                size_t maxscale_limits = min(maxscale_left, maxscale_right);
                maxscale_absolute = min(maxscale_limits, maxscale);
            }

            // var for access in inner loop
            size_t offset_prev = 0;
            for (size_t scale = GLOBAL_MINSCALE; scale <= maxscale_absolute; scale++)
            {
                // random access is difficult to vectorise
                double leftVal = cen[-scale];
                double rightVal = cen[scale];

                // expand the kernel to the left and right of the intensity_log.
                // b0 is expanded by the two outer points * 1
                product_sum_b0 += leftVal + rightVal;
                // b1 is expanded by the points * scale, negative to the left
                // product_sum_b1 += -double(scale) * leftVal + double(scale) * rightVal;
                product_sum_b1 += double(scale) * (rightVal - leftVal);
                // b2 and b3 are expanded by scale^2 the outermost point to the left or right
                double scale_sqr = double(scale * scale);
                product_sum_b2 += scale_sqr * leftVal;
                product_sum_b3 += scale_sqr * rightVal;

                const MatInverse inv = qalgo_matInverse[scale];

                const double inv_B_b0 = inv.B * product_sum_b0;
                const double inv_D_b1 = inv.D * product_sum_b1;

                // access is determined by scale and x0.
                // scale 2: idx is x0 - scale
                // scale 3: length - ((scale - 1) * 2) + x0 - scale
                // scale 4: 2 * length - (scale - 1) * 2 - (scale - 2) * 2 + x0 - scale
                //          2 * length - scale * 4 + 6 + x0 - scale
                const size_t offset_front = x0 - scale;
                const size_t offset = offset_front + offset_prev;
                assert(offset < totalRegs);
                RegCoeffs *current = coeffs->data() + offset;
                current->b0 = inv.A * product_sum_b0 + inv.B * (product_sum_b2 + product_sum_b3);
                current->b1 = inv.C * product_sum_b1 + inv.D * (product_sum_b2 - product_sum_b3);
                current->b2 = inv_B_b0 + inv_D_b1 + inv.E * product_sum_b2 + inv.F * product_sum_b3;
                current->b3 = inv_B_b0 - inv_D_b1 + inv.F * product_sum_b2 + inv.E * product_sum_b3;
                current->scale = scale;
                current->x0 = x0;
                // next scale would access front of vector
                offset_prev += length - 2 * scale;
            }
        }
    }

#pragma endregion "running regression"

#pragma region "validate Regression"

    static double correctB0(const float *const intensities,
                            const Span_i32 range,
                            float *predicted,
                            RegCoeffs *coeff) // @todo also correct the error matrix, assess degree of correction applied
    {
        // problem: after the log transform, regression residuals are not directly transferable to
        // the retransformed model. This is corrected by adjusting b0 so that the MSE in the
        // exponential case is minimal. We perform another regression taking e^(x b1 + x^2 b23) as a constant
        // and adjusting e^b0 so that (intensities - e^b0 * constant)^2 is minimal.
        // This is the same as scaling the regression by a constant c.

        // The (XtX)^-1 term collapses into a scalar, 1 / sum(predict^2). Xt * y is also scalar, sum(predict * intensities).
        // The corrective factor is then sum(predict * intensities) / sum(predict^2).

        // predict intensity only within range to prevent unnecessary exp operations.
        // prediction is incomplete, so it has to be multiplied with exp(b0)

        // uncertainties of this prediction are derived from the general equation for every coefficient.
        // the derivative is the following general form where bx' depends on the derived coefficient (1 for b0, x for b1, x^2 for b23)
        // b0' + ln(sum(exp(b0 + b1 x + b23 x^2) * y_i)) * sum(exp(b0 + b1 x + b23 x^2) * y_i * bx')
        //     - ln(sum(exp(b0 + b1 x + b23 x^2)^2)) * sum(exp(b0 + b1 x + b23 x^2)^2 * bx'^2)
        // note that the value of bx' depends on x and cannot be removed from the sum. Precalculation is also not
        // possible due to the value of y being required. This is not a large performance concern since we can combine it
        // with the calculations that are already necessary for obtaining the modified b0.

        // double b0_old = coeff->b0;

        const size_t start = range.startIdx;
        const size_t end = range.endIdx() + 1;
        assert(start < end);
        double x0 = double(coeff->x0);
        for (size_t i = start; i < end; i++)
        {
            double x = double(i) - x0;
            predicted[i] = (float)regExp_fac(coeff, x);
        }

        double sum_predictSq = 0;
        double sum_predictReal = 0;
        double b0_exp = exp(coeff->b0);

        // Regression correction is only calculated from the range in which the regression is relevant initially.
        for (size_t i = start; i < end; i++)
        {
            double pred = predicted[i] * b0_exp;
            sum_predictSq += pred * pred;
            sum_predictReal += pred * intensities[i];
        }
        double correction = sum_predictReal / sum_predictSq;

        // exp(a) * exp(b) == exp(a + b), so b0 + log(correction) is the same as predict * correction
        coeff->b0 += log(correction);
        // assert(abs(coeff->b0 - b0_old) < 0.001);

        // adjust the now incorrect values for predict. Remember that the previous prediciton was incomplete!
        float factor = (float)exp(coeff->b0);
        for (size_t i = start; i < end; i++)
        {
            predicted[i] *= factor;
        }

        return factor;
    }

    static Span_i32 calcRegSpan(const RegCoeffs *coeffs)
    {
        // test regression validity without depending on b0 or the degrees of freedom
        const bool valley_left = coeffs->b2 >= 0;
        const bool valley_right = coeffs->b3 >= 0;
        const bool apexLeft = coeffs->b1 < 0;

        // the peak must have a maximum
        bool noApex = (valley_left && apexLeft) ||
                      (valley_right && (!apexLeft));
        if (noApex)
            return {0, 0};

        Span_i32 span = {0, 0};

        if (valley_left)
        {
            int32_t position_b2 = (int32_t)coeffs->x0 - int32_t(coeffs->b1 / (2 * coeffs->b2));
            int32_t position_sc = (int32_t)coeffs->x0 - coeffs->scale;
            span.startIdx = max(position_b2, position_sc);
        }
        else
        {
            span.startIdx = (int32_t)coeffs->x0 - (int32_t)coeffs->scale;
        }

        if (valley_right)
        {
            int32_t position_b3 = (int32_t)coeffs->x0 + int32_t(-coeffs->b1 / (2 * coeffs->b3));
            int32_t position_sc = (int32_t)coeffs->x0 + coeffs->scale;
            span.set_endIdx(min(position_b3, position_sc));
        }
        else
        {
            span.set_endIdx((int32_t)coeffs->x0 + (int32_t)coeffs->scale);
        }
        return span;
    }

    // static double apexToEdgeRatio(const RegressionGauss *mutateReg, const float *intensities);

    /// @brief calculate the residual sum of squares for the log regression / data
    /// @param mutateReg relevant regression
    /// @param observed log data
    /// @return RSS value
    static double calcRSS_log(const RegressionGauss *mutateReg, const float *observed);

    /// @brief performs two F-tests against the log data. First H0 is the mean, second y = mx + b
    /// @param observed log data (or normal data, depends on the use case)
    /// @param RSS_reg previously calculated residual sum of squares of the complex model. Hard assumpion of four coefficients.
    /// @param range range of the regression.
    /// @return true: Regression is significant; false: Regression is not better than either alternative.
    static bool f_testRegression(const float *intensities_log, double RSS_reg, const Span_i32 range);

    /// @brief Assign the regression area, height, uncertaincies and similar. further, statistical
    /// tests are used to determine if the given regression is a valid peak description
    /// @param intensities The initially measured data
    /// @param intensities_log log transform of the above
    /// @param x_axis what the name implies
    /// @param predict regression evaluated at all points of the integer-transformed x
    /// @param df_sum total degrees of freedom for the considered regression
    /// @param length number of elements in intensities(_log), x_axis and predict
    /// @param mutateReg the regression being evaluated. It is only fully initialised if the
    /// function returns OK
    /// @return success or type of failure of the executed tests
    static invalid calcRegressionProperties( // returns the number of the failed test
        const float *intensities,
        const float *intensities_log,
        const float *x_axis,
        const float *predict,
        const size_t df_sum,
        const size_t length,
        RegressionGauss *mutateReg)
    {
        assert(!mutateReg->isValid);
        const size_t scale = mutateReg->coeffs.scale;
        const RegCoeffs *coeffs = &mutateReg->coeffs;
        const Span_i32 regSpan = mutateReg->span;

        assert(scale > 1);
        assert(coeffs->x0 + scale < length);

        // @todo new error correction here. Previous covariance matrix was mse (log) * (XtX)^-1,
        // multiply with matrix U, where first four terms are partial derivative of equation in
        // correctB0 by original coefficients

        float position = (float)peakPosition(coeffs);
        // to prevent potential errors, reject those regressions which do not have the apex in the
        // regression window with at least one point tp the side.
        if ((int32_t)abs(position) + 2 > coeffs->scale)
        {
            return invalid_apex;
        }

        // @todo this should be a test for similarity, since we cannot assume factor two is a good
        // fit for the worst possible combination of apex and edge, t-test instead?
        // double apexToEdge = apexToEdgeRatio(mutateReg, intensities);
        // if (apexToEdge < 2)
        // {
        //     failstates += 1;
        //     // return invalid_apexToEdge; // invalid apex to edge ratio // b0 independent
        // }

        // @todo differentiate between tests performed in log and exp domain better
        // everything involving the RSS is dependent on b0!
        double RSS_log = calcRSS_log(mutateReg, intensities_log);
        assert(RSS_log > 0);
        // double RSS_exp = calcRSS(predict->data(), intensities, &regSpan);
        // assert(RSS_exp > 0);
        // double mse_exp = RSS_exp / double(df_sum);
        double mse_log = RSS_log / double(df_sum);

        bool badCoeff = insignificantCoefficient(coeffs, mse_log, df_sum);
        if (badCoeff)
        {
            return invalid::invalid_quadratic;
        }

        /*
        We use an F-test instead of a normal height threshold. If the peak shape is a good enough
        fit, even very low intensity signals should be accounted for. The test here has a hard-coded
        alpha of 0.05
        */
        bool f_ok = f_testRegression(intensities_log, RSS_log, regSpan);
        if (!f_ok)
        {
            return f_test_fail; // H0 holds, the two distributions are not noticeably different
        }

        // uncertainty calculation and t-tests against peak properties

        float uncert_position = (float)peakPositionUncert(coeffs, mse_log);
        mutateReg->position_unc = uncert_position;
        float uncert_height = (float)peakHeightUncert(coeffs, mse_log);
        mutateReg->height_unc = uncert_height;

        const float absolutePosition = (float)coeffs->x0 + position;
        if (!std::isfinite(absolutePosition) || absolutePosition < 0 || absolutePosition + 1 >= length)
            return invalid_apex;
        size_t leftOfApex = (size_t)absolutePosition;
        float delta_x = x_axis[leftOfApex + 1] - x_axis[leftOfApex];
        double uncertainty = -1;
        mutateReg->area = (float)peakArea(coeffs, delta_x, mse_log, &uncertainty);
        mutateReg->area_unc = (float)uncertainty;

        // Note: This is not a satistically reasoned limit, but rather the number e regardless of
        // scaling. We introduce this limit here to prevent cases where the log operation introduces
        // negative intensities
        if (mutateReg->area <= minIntensity_global)
            return invalid_height;

        // @todo cramer-von-mises test / modification of KS statistic instead. The property of interest
        // is correctness of the peak shape as a volume, not strictly the fitted curve

        mutateReg->df = df_sum;
        mutateReg->position = (float)coeffs->x0 + position;
        mutateReg->jaccard = (float)calcJaccardIdx(intensities, predict, length);
        mutateReg->height = (float)exp(regAt(coeffs, peakPosition(coeffs)));
        assert(mutateReg->position > 1);

        mutateReg->isValid = false; // note: this field is only required later on for eliminating regressions
        return ok;
    }

    double calcRSS_log(const RegressionGauss *mutateReg, const float *observed)
    {
        double RSS = 0;
        const size_t start = mutateReg->span.startIdx;
        const float *obs = observed + start;
        double x = double(start) - double(mutateReg->coeffs.x0);
        const size_t len = mutateReg->span.length;
        for (size_t i = 0; i < len; i++)
        {
            double pred = regAt(&mutateReg->coeffs, x);
            double difference = obs[i] - pred;
            x += 1.0;
            RSS += difference * difference;
        }
        assert(RSS > 0);
        return RSS;
    }

    static double calcRSS_H0_cf2(const float *observed, const Span_i32 range)
    {
        // this function calculates the RSS for H0: y = b0 + x * b1 (no weights)

        double slope = NAN;
        double intercept = NAN;
        size_t length = range.length;
        assert(length > 0);
        const float *obs = observed + range.startIdx;
        linReg_intx(obs, length, &slope, &intercept);

        double RSS = 0;
        double x = 0;
        for (size_t i = 0; i < length; i++)
        {
            // linReg_intx returns the slope of y = intercept + slope * x.
            // Reversing its sign inflates the null-model RSS and admits
            // regressions which do not improve on a simple linear trend.
            double difference = obs[i] - (intercept + slope * x);
            RSS += difference * difference;
            x += 1;
        }
        assert(RSS >= 0);

        return RSS;
    }

    bool f_testRegression(const float *intensities_log, double RSS_reg_log, const Span_i32 range)
    {
        // during the tests, the RSS for the regression has already been calculated in calcRSS_log
        assert(RSS_reg_log > 0);
        const size_t length = range.length;
        /*
            The following code was initially intended to speed up the f test by performing a simple
            calculation first and then skipping if the result is invalid. Obviously, it is actually
            incorrect to test the same property twice. Further, due to the massive difference in
            degrees of freedom between both cases for regressions of size 5 (that being 3), sometimes
            the constant model was preferable despite the intended linear model being worse than a fit

            bool f_ok = false;
            double RSS_H0_cf1 = calcRSS_H0_cf1(intensities_log, range); // y = b
            f_ok = F_test_regs(RSS_reg_log, RSS_H0_cf1, 4, 1, length, 0.05);
            if (!f_ok)
                return false;
        */

        double RSS_H0_cf2 = calcRSS_H0_cf2(intensities_log, range); // y = b0 + b1 * x
        // A perfect linear fit is valid and cannot be improved by adding
        // parameters. This also keeps zero null-model RSS out of F_value.
        if (RSS_H0_cf2 <= RSS_reg_log)
            return false;
        bool f_ok = F_test_regs(RSS_reg_log, RSS_H0_cf2, 4, 2, length, 0.05);

        return f_ok;
    }

#if 0
    static double apexToEdgeRatio(const RegressionGauss *mutateReg, const float *intensities)
    {
        // is the apex at least twice as large as the outermost point?
        // assumption: outermost point is already near base level
        // addition: this relation should hold for the observed and predicted  values, also compared
        // to each other. Without this filter, there are some cases during feature detection where
        // a given bin contains no peak, but the model is fit to an ascending slope initially and
        // a short, nearly constant positive-coefficient half on the right. As such, the condition
        // of signal-to-noise (approximately) has to hold for the worst possible combination of all
        // candidate values
        const size_t idxApex = size_t(mutateReg->position) + mutateReg->coeffs.x0;

        double maxEdge_O = max(intensities[mutateReg->span.startIdx],
                               intensities[mutateReg->span.endIdx()]);
        double apex_O = intensities[idxApex]; // since this is not the actual apex height, it might be a bad idea to use it

        double x_l = mutateReg->span.startIdx - mutateReg->coeffs.x0;
        double x_r = x_l + mutateReg->span.length;
        double maxEdge_P = max(exp(regAt(&mutateReg->coeffs, x_l)),
                               exp(regAt(&mutateReg->coeffs, x_r)));

        double apex = min(apex_O, exp(regAt(&mutateReg->coeffs, mutateReg->position)));
        double maxEdge = max(maxEdge_O, maxEdge_P);

        return apex / maxEdge;
    }
#endif

    static bool insignificantCoefficient(const RegCoeffs *coeffs, const double mse, const size_t df_sum)
    {
        // If the coefficient for the peak half containing the apex is not significantly
        // different from 0, the produced regression does not reliably describe a system
        // that actually has an apex. On the other hand, the non-apex half can approach
        // zero without impacting the regression correctness. Significance here uses a
        // one-sided t-Test against zero.

        assert(mse > 0);
        double testCoeff = abs(coeffs->b1 < 0 ? coeffs->b2 : coeffs->b3); // absolute difference counts
        const double inv_E = qalgo_matInverse[coeffs->scale].E;
        double divisor = sqrt(inv_E * mse);
        double t_value = testCoeff / divisor; // test against mean 0
        return t_value < t_critical_95(df_sum);
    }

#pragma region "Feature Detection"

    static float weightedMeanAndVariance_EIC(const float *weight,
                                             const float *values,
                                             const size_t length,
                                             float *variance)
    {
        assert(weight[0] > 0);
        assert(values[0] > 0);
        assert(weight[length - 1] > 0);
        assert(values[length - 1] > 0);

        // weighted mean using intensity as weighting factor and left_limit right_limit as range
        size_t realPoints = 0;
        double sum_weighted_x = 0; // sum of values * weight
        double sum_weight = 0;
        for (size_t j = 0; j < length; j++)
        {
            // multiplication with zero is used instead of a continue so this can be vectorised.
            bool measured = values[j] != 0;
            float w = (float)measured * weight[j];
            sum_weighted_x += values[j] * w;
            sum_weight += w;
            realPoints += (size_t)measured; // interpolated points do not count!
        }
        double dpoints = (double)realPoints;

        double weighted_mean = sum_weighted_x / sum_weight;

        // only calculate the mean if the variance is not required
        if (variance == nullptr)
            return (float)weighted_mean;

        double sum_Qxxw = 0.0; // sum of (values - mean)^2 * weight
        for (size_t j = 0; j < length; j++)
        {
            float val = values[j];
            double difference = val - weighted_mean;
            sum_Qxxw += difference * difference * weight[j] * float(val != 0);
        }

        *variance = (float)sqrt(sum_Qxxw / sum_weight / dpoints);
        return (float)weighted_mean;
    };

    static FeaturePeak peakToFeat(const RegressionGauss *peak, const EIC *eic, uint32_t eic_ID)
    {
        const size_t start = peak->span.startIdx;
        const float *area_arr = eic->ints_area.data() + start;
        const float *mz_arr = eic->mz.data() + start;
        const float *dqsc_arr = eic->DQSC.data() + start;
        const float *rt_arr = eic->RT.data() + start;

        float DQSB = -1; // = weightedMeanAndVariance_EIC(area_arr, dqsb_arr, peak->length, nullptr);
        float DQSC = weightedMeanAndVariance_EIC(area_arr, dqsc_arr, peak->span.length, nullptr);
        assert(DQSC > 0);
        assert(DQSC <= 1);

        float mz_uncert = 0;
        float mz = weightedMeanAndVariance_EIC(area_arr, mz_arr, peak->span.length, &mz_uncert);
        assert(mz > 10);
        assert(mz_uncert >= 0);

        return {
            peak->coeffs,
            peak->height,
            peak->area,
            peak->height_unc,
            peak->area_unc,
            peak->dqs,
            DQSB,
            DQSC,
            peak->position,
            mz,
            peak->position_unc,
            mz_uncert,
            eic_ID,
            (uint32_t)peak->span.startIdx,
            (uint32_t)peak->span.length,
            rt_arr[0],
            rt_arr[peak->span.length - 1],
        };
    }

    size_t findFeatures(const std::vector<EIC> *EICs,
                        const std::vector<float> * /*convertRT*/,
                        std::vector<FeaturePeak> *res)
    {
        std::vector<RegressionGauss> peaks;

        // the only relevant change here is that the start point of the x axis (RT) depends on the first scan in an EIC
        // it is already assured that a bin contains continuous data
        size_t eicCount = EICs->size();
        for (size_t eic_ID = 0; eic_ID < eicCount; eic_ID++)
        {
            const EIC *bin = EICs->data() + eic_ID;
            const float *rt = bin->RT.data();
            size_t binLen = bin->df.size();
            int32_t status = 0;
            try
            {
                status = qpeaks_find(bin->ints_area.data(), rt, bin->df.data(), binLen, &peaks);
            }
            catch (const std::exception &error)
            {
                throw std::runtime_error("Feature detection failed in EIC " + std::to_string(eic_ID) + ": " + error.what());
            }
            if (status < 0)
                throw std::runtime_error("Feature detection rejected EIC " + std::to_string(eic_ID) + ".");
            for (size_t peak = 0; peak < peaks.size(); peak++)
            {
                FeaturePeak feat = peakToFeat(&peaks[peak], bin, eic_ID);
                res->push_back(feat);
            }
            peaks.clear();
        }

        // peaks are sorted here so they can be treated as const throughout the rest of the program
        std::sort(res->begin(), res->end(), [](const FeaturePeak &lhs, const FeaturePeak &rhs)
                  { return lhs.retentionTime < rhs.retentionTime; });

        return res->size();
    }

#pragma endregion "Feature Detection"

#pragma region "find centroids"

    static inline CentroidPeak regToCen(const RegressionGauss *reg, uint32_t id, uint16_t specNum)
    {
        return {
            reg->position,
            reg->height,
            reg->area,
            reg->fwhm,
            reg->height_unc,
            reg->area_unc,
            reg->position_unc,
            reg->dqs,
            id,
            specNum,
            reg->coeffs.scale,
            reg->numCompetitors,
        };
    }

    size_t findCentroids(const XML_File *data,
                         const std::vector<uint32_t> *selectedIndices,
                         std::vector<CentroidPeak> *centroids)
    {
        assert(!data->defective);
        assert(centroids->empty());
        if (selectedIndices->size() > UINT16_MAX)
            throw std::runtime_error("More than 65535 MS1 profile scans per polarity are not supported in this working version.");

        std::vector<float> spectrum_mz;
        spectrum_mz.reserve(1000);
        std::vector<float> spectrum_int;
        spectrum_int.reserve(1000);
        std::vector<RegressionGauss> ret;

        for (size_t specNum = 0; specNum < selectedIndices->size(); specNum++)
        {
            int32_t ok = get_spectrum(data,
                                      &spectrum_mz,
                                      &spectrum_int,
                                      selectedIndices->at(specNum));
            if (ok != 0)
                throw std::runtime_error("Could not decode an mzML spectrum; this polarity's feature list was not completed.");
            if (spectrum_int.size() < MINLENGTH)
            {
                spectrum_mz.clear();
                spectrum_int.clear();
                continue;
            }

            int32_t peaksFound = 0;
            try
            {
                peaksFound = qpeaks_find(spectrum_int.data(), spectrum_mz.data(), nullptr, spectrum_int.size(), &ret);
            }
            catch (const std::exception &error)
            {
                throw std::runtime_error("Centroid detection failed at spectrum index " +
                                         std::to_string(selectedIndices->at(specNum)) + ": " + error.what());
            }
            if (peaksFound < 0)
                throw std::runtime_error("Centroid detection rejected a spectrum.");

            for (size_t p = 0; p < (size_t)peaksFound; p++)
            {
                const RegressionGauss *peak = ret.data() + p;
                size_t id = centroids->size();
                if (id >= UINT32_MAX)
                    throw std::runtime_error("The centroid count exceeds the supported range.");
                centroids->push_back(regToCen(peak, id, specNum));
            }

            spectrum_mz.clear();
            spectrum_int.clear();
            ret.clear();
        }
        return centroids->size();
    }

    void centroids_to_mzml(const std::filesystem::path *pathSource)
    {
        (void)fprintf(stderr, "Warning: Due to the way processing is handled internally, it is not possible\n"
                              "to use qAlgorithms for centroiding to mzML and to produce feature lists in one run.\n");

        std::filesystem::path pathTarget = *pathSource;
        pathTarget.replace_filename(pathSource->stem().string() + "_qcentroid.mzML");

        // only work on the copy
        XML_File source_file(pathSource->c_str(), mzML);
        const size_t numSpecs = source_file.number_spectra;

        std::vector<float> spectrum_mz;
        std::vector<float> spectrum_int;
        std::vector<RegressionGauss> result;
        std::vector<double> res_spectrum_mz;
        std::vector<double> res_spectrum_int;
        std::vector<char> char_spectrum_mz;
        std::vector<char> char_spectrum_int;

        for (size_t specNum = 0; specNum < numSpecs; specNum++)
        {
            // 1) obtain next spectrum that is MS level 1 and not centroided
            const pugi::xml_node *spec = source_file.linknodes->data() + specNum;

            bool isProfile = spectrum_is_profile(&source_file, specNum);
            if (!isProfile)
                continue;

            uint32_t ms_level = spectrum_ms_level(&source_file, specNum);
            if (ms_level != 1)
                continue; // centroiding only makes sense for ms level 1

            get_spectrum(&source_file, &spectrum_mz, &spectrum_int, specNum);

            // 3) extract spectum data, process with qpeaks_find
            int32_t numPeaks = qpeaks_find(spectrum_int.data(),
                                           spectrum_mz.data(),
                                           nullptr,
                                           spectrum_mz.size(),
                                           &result);
            assert(numPeaks > 0);
            const size_t numPeaks_u = (size_t)numPeaks;

            // 4) write centroid m/z and intensity into two float arrays
            res_spectrum_int.resize(numPeaks_u);
            res_spectrum_mz.resize(numPeaks_u);
            for (size_t peak = 0; peak < numPeaks_u; peak++)
            {
                res_spectrum_int[peak] = result[peak].area;
                res_spectrum_mz[peak] = result[peak].position;
            }

            // 5) compress arrays to byte vector
            compress_and_encode((uint8_t *)res_spectrum_int.data(),
                                numPeaks * sizeof(float),
                                &char_spectrum_int);
            compress_and_encode((uint8_t *)res_spectrum_mz.data(),
                                numPeaks * sizeof(float),
                                &char_spectrum_int);

            // 6) set spectrum mode to centroided. This means replacing the attribute name and value
            spec->find_child_by_attribute("cvParam", "name", "profile spectrum").set_value("centroid spectrum");
            spec->find_child_by_attribute("cvParam", "accession", "MS:1000128").set_value("MS:1000127");

            // 7) overwrite intensity and m/z for spectrum, also set length of array
            pugi::xml_named_node_iterator dataArray = spec->child("binaryDataArrayList").children("binaryDataArray").begin();
            dataArray->attribute("encodedLength").set_value(char_spectrum_mz.size());
            dataArray->child("binary").set_value(char_spectrum_mz.data());
            dataArray->next_sibling();
            dataArray->attribute("encodedLength").set_value(char_spectrum_int.size());
            dataArray->child("binary").set_value(char_spectrum_int.data());

            // reset storage vectors to avoid in-loop allocation
            spectrum_mz.clear();
            spectrum_int.clear();
            result.clear();
            res_spectrum_mz.clear();
            res_spectrum_int.clear();
            char_spectrum_mz.clear();
            char_spectrum_int.clear();
        }
        source_file.mzml_base_document.save_file(pathTarget.c_str());
        source_file.free_linknodes();
    }

#pragma region "Helper math"

    inline double regAt(const RegCoeffs *coeff, const double x)
    {
        const double b23 = x < 0 ? coeff->b2 : coeff->b3;
        return coeff->b0 + (coeff->b1 + x * b23) * x;
    }

    double regExp_fac(const RegCoeffs *coeff, const double x)
    {
        double b23 = x < 0 ? coeff->b2 : coeff->b3;
        return exp((coeff->b1 + x * b23) * x);
    }

    double fullWidthHalfMax(const RegCoeffs *coeff, const double height, const double delta_x)
    {
        assert(height > 2);
        // solve height / 2 = exp(b0 + x b1 + x^2 b2)
        double y = log(height / 2);
        // use quadratic formula for a x^2 + b x + c = 0
        double a_l = coeff->b2;
        double a_r = coeff->b3;
        double b = coeff->b1;
        double c = coeff->b0 - y;

        double x_l = NAN;
        double x_r = NAN;
        double dummy = NAN;
        solveQuadratic(a_l, b, c, &x_l, &dummy);
        solveQuadratic(a_r, b, c, &dummy, &x_r);

        // special case: valley is such that one half of the peak is not defined at half maximum

        if (x_l == INFINITY)
        {
            assert(a_l > 0);
            x_l = -b / (2 * a_l);
        }
        if (x_r == INFINITY)
        {
            assert(a_r > 0);
            x_r = -b / (2 * a_r);
        }

        assert(x_l < x_r);

        return (x_r - x_l) * delta_x;
    }

    double matProductReg(const double J[4], const size_t scale)
    {
        // Calculate the Matrix Product of J * Xinv * J^T for uncertainty calculation,
        // Where J is the jacobian matrix of the regression property A for which the
        // uncertainty should be calculated (J = [d A / d b0, d A / d b1, d A / d b2, d A / d b3])

        ///                   | A  0  B  B |    | J0 |
        ///                   | 0  C  D -D |    | J1 |
        /// | J0 J1 J2 J3 | x | B  D  E  F | x  | J2 |
        ///                   | B -D  F  E |    | J3 |
        ///
        ///                                                                    | J0 |
        /// |   J0 A + J1 0;   J0 0 + J1 C;   J0 B + J1 D;   J0 B - J1 D |     | J1 |
        /// | + J2 B + J3 B; + J2 D - J3 D; + J2 E + J3 F; + J2 F + J3 E |  x  | J2 |
        ///                                                                    | J3 |
        /// The product is a scalar
        ///
        ///   J0^2 A + J0 J2 B + J0 J3 B
        /// + J1^2 C + J1 J2 D - J1 J3 D
        /// + J0 J2 B + J1 J2 D + J2^2 E + J2 J3 F
        /// + J0 J3 B - J1 J3 D + J2 J3 F + J3^2 E
        ///
        /// simplify
        ///
        ///   J0^2 A + J1^2 C + (J2^2 + J3^2) E
        /// + J0 B * (J2 + J3) + J1 D * (J2 - J3)
        /// + J0 B * (J2 + J3) + J1 D * (J2 - J3)
        /// + 2 * J2 J3 F
        ///
        ///   J0^2 A + J1^2 C + (J2^2 + J3^2) E
        /// + 2 * J0 B * (J2 + J3) + 2 * J1 D * (J2 - J3) + 2 * J2 J3 F
        ///
        /// = J0^2 A + J1^2 C + (J2^2 + J3^2) E + 2 * (J0 B * (J2 + J3) + J1 D * (J2 - J3) + J2 J3 F )

        assert(scale > 1);
        const MatInverse inv = qalgo_matInverse[scale];
        double vecMatrxTranspose = J[0] * J[0] * inv.A +
                                   J[1] * J[1] * inv.C +
                                   (J[2] * J[2] + J[3] * J[3]) * inv.E +
                                   2 * (J[0] * inv.B * (J[2] + J[3]) + J[1] * inv.D * (J[2] - J[3]) + J[2] * J[3] * inv.F);
        return vecMatrxTranspose;
    }

    double peakArea(const RegCoeffs *c, const double delta_x, const double mse, double *uncert)
    {
        // base function: integral of e^(b0 + b1 x + b2 x^2) dx =
        // [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * erfi( (b1 + 2 b2 x) / (2 sqrt(b2)) ) ] / (2 sqrt(b2))   // source: wolfram alpha
        // erfi(x) = i * erf(i * x) * -1
        // under the condition b2 < 0: sqrt(b2) = i * sqrt(-b2), where sqrt(-b2) is a real number
        // [(real part) * -1 * i * erf( i * (real) / (i * 2 * sqrt(-b2)) ) ] / (i * 2 * sqrt(-b2))
        // i within and outside of the error function cancel each other out: i / i = 1
        // erf(0) = 0 ; erf(-inf) = -1 ; erf(inf) = 1
        // F(-inf) = [ (...) * -1 * erf( (b1 + 2 b2 * -inf) / (> 0) ) ] / (...)
        // when evaluating F towards infinity, b2 or b3 always have a negative sign
        // b2 * -inf = +inf ; b2 * +inf = -inf
        // F(-inf) = [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * -1 * +1 ] / (2 sqrt(-b2))
        // F(+inf) = [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * -1 * -1 ] / (2 sqrt(-b2))
        // F(0)    = [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * -1 * erf( b1 / (2 sqrt(-b2)) ) ] / (2 sqrt(-b2))
        // if b2 or b3 are positive, erfi has to be used. The positive part of the function replaces F(+-inf)
        // for reasons of numerical stability, the constant scaling factor    b0_exp * sqrt(pi) / 2    is factored out.
        // the final equation is composed as follows:
        // A = b0_exp * sqrt(pi) / 2 * (F(0)_b2 - F(-inf)_b2 + F(inf)_b3 - F(0)_b3)
        // for every peak half, F can be written as exp_b23 * erf_b23 / sqrt_b23, where only erf_b23 changes. As such:
        // A = b0_exp * sqrt(pi) / 2 * (exp_b2 / sqrt_b2 * (erf_b2_0 - erf_b2_-inf) + exp_b3 / sqrt_b3 * (erf_b3_inf - erf_b3_0))
        // apply erf(-x) = -erf(x); erfc(-x) = erf(x) + 1:
        // left:  erf_b2_0 - erf_b2_-inf =    -erf(b1 / (2 sqrt(-b2)) + 1 = erfc( b1 / sqrt(-b2))
        // right: erf_b3_inf - erf_b3_0  = 1 + erf(b1 / (2 sqrt(-b3))     = erfc(-b1 / sqrt(-b3))
        // A = b0_exp * (sqrt(pi) / 2) * (
        //      e^(-b1^2/(4 b2)) * erfc( b1 / (2 sqrt(-b2)) ) / sqrt(-b2) +
        //      e^(-b1^2/(4 b3)) * erfc(-b1 / (2 sqrt(-b3)) ) / sqrt(-b3)
        // )
        // erfcx(x) = e^(x^2) * erfc(x); (b1 / (2 sqrt(-b2)))^2 = -b1^2 / (4 b2); (-b1 / (2 sqrt(-b2)))^2 = -b1^2 / (4 b2)
        // A = b0_exp * (sqrt(pi) / 2) * (
        //      erfcx( b1 / (2 sqrt(-b2)) ) / sqrt(-b2) +
        //      erfcx(-b1 / (2 sqrt(-b3)) ) / sqrt(-b3)
        // )
        // in implementation, sqrt(pi) / 2 is multiplied with each half individually so that the case
        // of a positive coefficient is handled cleanly.

        // if there is a valley point, the antiderivative evaluated at the valley is always zero. This
        // is because erfi(b1 + 2 b2 x) resolves to erfi(b1 + 2 b2 * (-b1 / (2 b2))) == 0
        // (the valley point is always at -b1 / (2 b2), refer to apex position calculation).
        // The area of the half is then:
        // A_L = e^(b0) * sqrt(pi) / 2 * e^(-b1^2/(4 b2)) * erfi(b1 / (2 sqrt(b2)))  / sqrt(b2)
        // This suffers from the same algorithmic instability as the positive case.
        // Here, we can use Dawson's integral D(x) = 1/2 sqrt(pi) * e^(-x^2) * erfi(x)
        // A_L = e^(b0) * D(b1 / (2 sqrt(b2))) / sqrt(b2)
        // The same transformation applies to b3. Depending on left / right valley position,
        // we calculate the partial area as A_L - 0 or 0 - A_R

        // further consideration has to be made for the case of a peak with the non-apex coefficient
        // being effectively zero. This is determined to be the case if the peak area of that half
        // is infinity - as can be seen above, the error function of a very large number is divided
        // by a very small number for b2 or b3 close to zero. No compensation in the code here since
        // there never were any cases of this occuring, add in case it happens

        const double b0 = c->b0;
        const double b1 = c->b1;
        const double b2 = c->b2;
        const double b3 = c->b3;
        const double sqrt_pi_2 = 1.7724538509055158819 / 2; // sqrt(M_PI);

        bool b2_neg = b2 < 0;
        bool b3_neg = b3 < 0;
        assert(b2_neg || b3_neg);
        // the apex cannot be left if b2 is positive, same for b3
        // (NOT b2_neg) NAND (b1 < 0)
        // NOT ((NOT b2_neg) AND (b1 < 0))
        // b2_neg OR (NOT b1 < 0) // b1 != 0
        // b2_neg OR b1 > 0
        assert(b2_neg || b1 > 0);
        assert(b3_neg || b1 < 0);
        assert(delta_x > 0);

        double sqrt_b2 = sqrt(abs(b2));
        double sqrt_b3 = sqrt(abs(b3));
        double b0_exp = exp(b0);

        double area_L = b2_neg ? libcerf::erfcx(b1 / sqrt_b2 / 2) / sqrt_b2 * sqrt_pi_2
                               : libcerf::dawson(b1 / sqrt_b2 / 2) / sqrt_b2;
        assert(area_L > 0 && area_L < INFINITY);

        double area_R = b3_neg ? libcerf::erfcx(-b1 / sqrt_b3 / 2) / sqrt_b3 * sqrt_pi_2
                               : -libcerf::dawson(b1 / sqrt_b3 / 2) / sqrt_b3;
        assert(area_R > 0 && area_R < INFINITY);

        double area_F = (area_L + area_R) * b0_exp;

        if (uncert == nullptr) // not calculating the uncertainty means mse can be ignored
        {
            return area_F * delta_x;
        }
        assert(mse > 0);

        {
            // We cannot directly calculate the uncertainty from the exponential form (true peak area).
            // This is because the mse required for the uncertainty is only defined as relating to the
            // coefficients in the log-transformed space.
            // Instead, since the coefficients apply to the logarithmic transform, we calculate the
            // uncertainty of the logarithm of the area as defined above.
            //     A  = e^b0 *  (area_L + area_R) * deleta_x
            // log(A) = b0 + log(area_L + area_R) + deleta_x
            // The uncertainty is calculated using the Jacobian J, which contains the per-coefficient
            // partial derivatives. delta_x is just scaling and never relevant to uncertainty.
            // It is obvious that the derivative by b0 is 1.
            double J[4];
            J[0] = 1;

            // b1:
            // The derivative of log(x) is x' / x. Apply the chain rule and we obtain:
            // log(A)'_b1 = (area_L' + area_R') / (area_L + area_R)
            // For negative b2, the area is calculated as:
            // area_L  = erfcx(b1 / (2 sqrt(b2))) * sqrt(pi) / (2 sqrt(b2))
            // area_L' = -H_(-2)(b1 / (2 sqrt(b2))) * 1/b2
            // Where H_(-2) is the hermite polynomial evaluated at -2. The following transformation applies:
            // H_(-2)(x) = 1/2 - sqrt(pi)/2 * x * e^(x^2) * erfc(x) = 1/2 - sqrt(pi)/2 * x * erfcx(x)
            // (see https://math.stackexchange.com/questions/2238662/do-hermite-polynomials-exist-for-negative-integers)
            // replace H in the above equation:
            // area_L' = 1/(2 b2) - sqrt(pi)/2 * b1 / (2 sqrt(b2) * erfcx(b1 / (2 sqrt(b2)) * 1/b2
            // notice that we can extract the definition of the area from the derivative and obtain:
            // area_L' = (area_L * b1 - 1) / (2 b2)
            // For b2 > 0, area_L = D(b1 / 2 sqrt(b2)) / sqrt(b2), where D is dawson's integral.
            // D(a x)' = a - 2 a^2 x D(a x)
            // area_L' = 1 / (2 sqrt(b2)) * 1 / sqrt(b2) -
            //           2 ( 1 / (2 sqrt(b2)) )^2 *
            //           b1 * D(b1 / 2 sqrt(b2)) / sqrt(b2)
            // area_L' = 1 / (2 b2) -
            //           1 / (2 b2) *
            //           b1 * D(b1 / 2 sqrt(b2)) / sqrt(b2)
            // factor out the area itself and  1 / (2 b2):
            // area_L' = 1 / (2 b2) * (1 - b1 * area_L)
            // We see that the first derivative of area_L is a scaling of itself that is fully
            // independent of the sign of b2.

            // If we work with b3, we use erfcx(-x). This is equivalent to substituting b1 with -b1 and multiplying with -1
            // Consequently, we can back-substitute b1 = -b1 into the final result:
            // area_R' = 1 / (2 b3) * (1 + b1 * area_R) * -1
            // For positive values of b3, divide the full term by  -1 / sqrt(b3)  instead of +1:
            // area_R' = 1 / (2 sqrt(b3)) *  -1 / sqrt(b3) -
            //           2 ( 1 / (2 sqrt(b3)) )^2 * b1 *
            //           -D(b1 / 2 sqrt(b3)) / sqrt(b3)
            // area_R' = -1 / (2 b3) - 1 / (2 b3) * b1 * area_R
            // area_R' = -1 / (2 b3) * (1 + b1 * area_R)
            // Here, too, the two terms converge to the same scaling of the initial area

            double area_L_db1 = (1 - b1 * area_L) / (2 * b2);
            double area_R_db1 = -(1 + b1 * area_R) / (2 * b3);
            assert(area_L_db1 != NAN);
            assert(area_R_db1 != NAN);
            J[1] = (area_L_db1 + area_R_db1) / (area_L + area_R);

            // for the partial derivative by b2 and b3 for negative coefficents, there is a dependence on the
            // hermite polynomials H_(-2)(x) and H_(-1)(x):
            // area_L  = erfcx(a / (2 sqrt(x))) / sqrt(x) * sqrt(pi)/2 ; get first derivative with wolfram alpha
            // area_L' = ( a H_(-2)(a/(2 sqrt(x))) ) / (2 x^2) - ( H_(-1)(a/(2 sqrt(x))) ) / (2 x^(3/2))
            // area_L' = b1 * H_(-2)( b1 / (2 sqrt(b2)) ) / (2 b2^2) -
            //           H_(-1)( b1 / (2 sqrt(b2)) ) / (2 b2 * sqrt(b2))
            // multiply the second term with (2 sqrt(pi)) / (2 sqrt(pi)) to convert H_(-1)(x) to erfcx(x) * sqrt(pi)/2
            // area_L' = b1 * H_(-2)( b1 / (2 sqrt(b2)) ) / (2 b2^2) -
            //           erfcx( b1 / (2 sqrt(b2)) ) / (2 b2 * sqrt(b2)) * sqrt(pi)/2
            // simplyfy by applying the definition of the area to the second term
            // area_L' = b1 / (2 b2^2) * H_(-2)( b1 / (2 sqrt(b2)) ) -
            //           area_L / (2 b2)
            // apply the identity of H_(-2)(x) with x = b1 / (2 sqrt(b2)), see above
            // area_L' = b1 / (2 b2^2) * ( 1/2 - sqrt(pi)/2 * b1 / (2 sqrt(b2)) * erfcx(b1 / (2 sqrt(b2))) ) -
            //           area_L / (2 b2)
            // area_L' = b1 / (4 b2^2) -
            //           b1 / (2 b2^2) * sqrt(pi)/2 * b1 / (2 sqrt(b2)) * erfcx(b1 / (2 sqrt(b2))) -
            //           area_L / (2 b2)
            // extract the area from the second term
            // area_L' = b1 / (4 b2^2) -
            //           b1 / (2 b2^2) * b1 / 2 * area_L -
            //           area_L / (2 b2)
            // area_L' = b1 / (4 b2^2) -
            //           b1^2 / (4 b2^2) * area_L -
            //           area_L / (2 b2) * (2 b2) / (2 b2)
            // area_L' = b1 / (4 b2^2) -
            //           area_L * b1^2 / (4 b2^2) -
            //           area_L * (2 b2) / (4 b2^2)
            // area_L' = ( b1 - area_L * b1^2 - area_L * (2 b2) ) / (4 b2^2)
            // area_L' = ( b1 - area_L * (b1^2 + 2 b2) ) / (4 b2^2)

            // For positive b2 (using wolfram alpha):
            // area_L = D(a / (2 sqrt(x))) / sqrt(x); a = b2; x = b2
            // area_L' = ( (a^2 - 2 x) D(a/(2 sqrt(x))) )/(4 x^(5/2)) - a/(4 x^2)
            // area_L' = ((b1^2 - 2 b2) * area_L) / (4 b2^2) - b1 / (4 b2^2)
            // area_L' = ( area_L * (b1^2 - 2 b2) - b1 ) / (4 b2^2)

            // unlike for b1, the derivative changes for positive or negative b2.
            // This is not just a sign flip, notice that here we use a difference instead
            // of a sum. This discrepancy originates in the fact that we initially took
            // the absolute value of b2 for the negative case, but not the positive one.
            // @todo is it correct to not write "sqrt(-b2)" for the uncertainty?

            // for negative b3, set b1 = -b1:
            // area_R' = ( - b1 - area_R * (b1^2 + 2 b3) ) / (4 b3^2)
            // area_R' =  -( b1 + area_R * (b1^2 + 2 b3) ) / (4 b3^2)
            // for positive b3, multiply the whole term with -1:
            // area_R' = - ( area_R * (b1^2 - 2 b3) - b1 ) / (4 b3^2)
            // area_R' =   ( b1 - area_R * (b1^2 - 2 b3) ) / (4 b3^2)

            // @todo results are not fully symmetric - is there an error here?

            if (b2_neg)
            {
                J[2] = (b1 - area_L * (b1 * b1 + 2 * b2)) / (4 * b2 * b2);
            }
            else
            {
                J[2] = (area_L * (b1 * b1 - 2 * b2) - b1) / (4 * b2 * b2);
            }

            if (b3_neg)
            {
                J[3] = -(b1 + area_R * (b1 * b1 + 2 * b3)) / (4 * b3 * b3);
            }
            else
            {
                J[3] = (b1 - area_R * (b1 * b1 - 2 * b3)) / (4 * b3 * b3);
            }

            double u = matProductReg(J, c->scale);
            assert(u > 0);
            *uncert = sqrt(u * mse); // this is the uncertainty of the logarithmic area.
        }

        // scale by delta_x, the assumed difference between two neighbouring points in x.
        // this is done explicitly in this function since it is assumed when using the precalculated inverse.
        *uncert *= delta_x;
        return area_F * delta_x;
    }

    double peakPositionUncert(const RegCoeffs *c, const double mse)
    {
        // peak position: -b1 / (2 b23)
        // Jacobian (derivatives by b0 to b3) is 0 for b0 and b2 or b3, depending on which peak half contains the apex
        double J[4] = {0, 0, 0, 0};
        bool apex_left = c->b1 < 0;
        double b23 = apex_left ? c->b2 : c->b3;
        J[1] = -1 / (2 * b23); // dx d b1
        int32_t idx = apex_left ? 2 : 3;
        J[idx] = c->b1 / (2 * b23 * b23); // -b1 * -1

        double u = matProductReg(J, c->scale);
        assert(u > 0);
        double uncert = sqrt(u * mse);
        return uncert;
    }

    inline double peakPosition(const RegCoeffs *c)
    {
        double b23 = c->b1 < 0 ? c->b2 : c->b3;
        return -c->b1 / (2 * b23);
    }

    double peakHeightUncert(const RegCoeffs *c, const double mse)
    {
        // height: exp( b0 + b1 * -b1 / (2 b23) + b23 * (-b1 / (2 b23))^2 )
        // = exp( b0 + (-b1^2 / (2 b23)) + (b23 * b1^2 / (4 b23^2)) )
        // = exp( b0 + (-b1^2 / (2 b23)) * (2 / 2) + b1^2 / (4 b23) )
        // = exp( b0 - b1^2 / (4 b23) )
        // We only deal with the log form here: b0 - b1^2 / (4 b23)

        bool apex_left = c->b1 < 0;
        double b23 = apex_left ? c->b2 : c->b3;
        double b1_sq = c->b1 * c->b1;
        // double height = exp(c->b0 - b1_sq / (4 * b23));

        // uncertainty: calculate uncertainty for the logarithmic case, then transform back
        // note: this leads to two values for uncertainty
        double deriv_b0 = -b1_sq / (4 * b23);
        double deriv_b1 = c->b0 - c->b1 / (2 * b23); // b0 + apex position
        double deriv_b23 = c->b0 + b1_sq / (4 * b23 * b23);
        // assign jacobian
        double J[4] = {deriv_b0, deriv_b1, 0, 0};
        J[apex_left ? 2 : 3] = deriv_b23;

        // double J_exp[4] = {height, -c->b1 * height / (2 * b23), 0, 0};
        // J_exp[apex_left ? 2 : 3] = b1_sq * height / (4 * b23 * b23);
        // double u_exp = matProductReg(J_exp, c->scale);
        // double test_u = sqrt(u_exp * mse_exp); // both u_exp and the MSE are too large. The combination via square root is also wrong (?)

        double u = matProductReg(J, c->scale);
        assert(u > 0);
        double uncert = sqrt(u * mse);
        // assert(uncert == test_u); // calculating the error for the exponential case directly leads to wrong results
        return uncert;
    }

} // namespace qAlgorithms
