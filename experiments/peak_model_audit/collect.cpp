#include "qalgorithms_peak_search.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

namespace
{
    constexpr double fixedStep = .35;

    void number(std::ostream &out, double value)
    {
        if (std::isfinite(value))
        {
            out << value;
        }
        else
        {
            out << "null";
        }
    }

    template <class Values>
    void array(std::ostream &out, const Values &values)
    {
        out << '[';
        bool first = true;
        for (double value : values)
        {
            if (!first)
            {
                out << ',';
            }
            first = false;
            number(out, value);
        }
        out << ']';
    }

    std::vector<size_t> subsample(const std::vector<size_t> &indices, size_t limit)
    {
        if (indices.size() <= limit)
        {
            return indices;
        }
        std::vector<size_t> result;
        for (size_t k = 0; k < limit; ++k)
        {
            result.push_back(indices[(2 * k + 1) * indices.size() / (2 * limit)]);
        }
        return result;
    }

    void writeWindow(std::ostream &out, size_t scan, const char *cohort, size_t center, size_t scale,
                     const std::vector<float> &mz, const std::vector<float> &intensity, double globalGap)
    {
        std::vector<float> raw, logs;
        std::vector<double> coordinates, gaps;
        const char *reason = "eligible";
        for (size_t i = center - scale; i <= center + scale; ++i)
        {
            raw.push_back(intensity[i]);
            coordinates.push_back(mz[i]);
            if (!(intensity[i] > 0 && std::isfinite(intensity[i])))
            {
                reason = "nonpositive_or_nonfinite";
            }
            if (i > center - scale)
            {
                gaps.push_back(std::sqrt(static_cast<double>(mz[i])) - std::sqrt(static_cast<double>(mz[i - 1])));
            }
        }
        const double mean = std::accumulate(gaps.begin(), gaps.end(), 0.0) / gaps.size();
        double variance = 0;
        bool regular = globalGap > 0 && std::isfinite(globalGap);
        for (double gap : gaps)
        {
            variance += (gap - mean) * (gap - mean);
            regular = regular && std::isfinite(gap) && gap >= .55 * globalGap && gap <= 1.45 * globalGap;
        }
        const double cv = mean > 0 ? std::sqrt(variance / gaps.size()) / mean : std::numeric_limits<double>::quiet_NaN();
        regular = regular && std::isfinite(cv) && cv <= .20;
        if (std::string(reason) == "eligible" && !regular)
        {
            reason = "irregular_or_gapped_coordinates";
        }
        out << "{\"type\":\"window\",\"scan\":" << scan << ",\"cohort\":\"" << cohort
            << "\",\"center\":" << center << ",\"scale\":" << scale << ",\"reason\":\"" << reason << "\",\"spacing_cv\":";
        number(out, cv);
        out << ",\"raw\":";
        array(out, raw);
        out << ",\"mz\":";
        array(out, coordinates);
        if (std::string(reason) == "eligible")
        {
            for (float value : raw)
            {
                logs.push_back(static_cast<float>(std::log(static_cast<double>(value))));
            }
            HyperbolicModel::Geometry cache(scale, fixedStep);
            HyperbolicModel::Fit selected;
            HyperbolicModel::runningFit({logs.data()}, logs.size(), cache, [&](const HyperbolicModel::Fit &fit)
            {
                if (fit.center == scale && fit.scale == scale)
                {
                    selected = fit;
                }
            });
            const auto inference = inferHyperbolicWindow({raw.data()}, raw.size(), scale, scale, fixedStep);
            const auto properties = evaluateHyperbolicWindow(selected, 0, cache, raw.data(), raw.size());
            out << ",\"coefficients\":";
            const double coefficients[4] = {selected.locals[0].intercept, selected.sharedSlope,
                selected.locals[0].curvature, selected.locals[0].skew};
            array(out, coefficients);
            out << ",\"logs\":";
            array(out, logs);
            out << ",\"inference_status\":" << static_cast<int>(inference.status)
                << ",\"decision\":" << static_cast<int>(inference.groupDecision)
                << ",\"integrable\":" << properties.integrable
                << ",\"apex_supported\":" << properties.apexInObservedSupport
                << ",\"refit_valid\":" << properties.amplitude.valid << ",\"area_u\":";
            number(out, properties.areaU.value);
            out << ",\"area_se_u\":";
            number(out, inference.properties.empty() ? std::numeric_limits<double>::quiet_NaN() : inference.properties[0].areaUncertainty.value);
        }
        out << "}\n";
    }
}

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            throw std::invalid_argument("Usage: real_peak_audit INPUT.mzML OUTPUT.jsonl");
        }
        const std::filesystem::path input(argv[1]);
        XML_File file(input.c_str(), mzML);
        if (file.defective)
        {
            throw std::runtime_error("Invalid source mzML");
        }
        const auto scans = filter_spectra(&file, true, Polarities::positive, true);
        std::vector<size_t> scanIndices(scans.begin(), scans.end());
        const auto selected = subsample(scanIndices, 80);
        if (selected.empty())
        {
            throw std::runtime_error("No positive profile MS1 scans");
        }
        std::ofstream out(argv[2]);
        out << std::setprecision(17);
        std::vector<float> mz, intensity;
        for (size_t scan : selected)
        {
            if (get_spectrum(&file, &mz, &intensity, scan) != 0 || mz.size() != intensity.size())
            {
                throw std::runtime_error("Invalid selected spectrum");
            }
            if (mz.size() < 17)
            {
                out << "{\"type\":\"scan\",\"scan\":" << scan << ",\"status\":\"too_short\"}\n";
                continue;
            }
            std::vector<size_t> centers, maxima;
            std::vector<double> gaps;
            for (size_t i = 0; i < mz.size(); ++i)
            {
                if (!(mz[i] > 0 && std::isfinite(mz[i])))
                {
                    throw std::runtime_error("Invalid m/z");
                }
                if (i > 0)
                {
                    const double gap = std::sqrt(static_cast<double>(mz[i])) - std::sqrt(static_cast<double>(mz[i - 1]));
                    if (!(gap > 0))
                    {
                        throw std::runtime_error("Non-increasing coordinates");
                    }
                    gaps.push_back(gap);
                }
                if (i >= 8 && i + 8 < mz.size())
                {
                    centers.push_back(i);
                    if (intensity[i] >= 100 && intensity[i] > intensity[i - 1] && intensity[i] >= intensity[i + 1])
                    {
                        maxima.push_back(i);
                    }
                }
            }
            std::nth_element(gaps.begin(), gaps.begin() + gaps.size() / 3, gaps.end());
            const double globalGap = gaps[gaps.size() / 3];
            centers = subsample(centers, 64);
            maxima = subsample(maxima, 64);
            std::vector<PeakSearchWindow> windows;
            for (size_t center : centers)
            {
                for (size_t scale : {4u, 8u})
                {
                    windows.push_back({center, scale});
                }
            }
            const HyperbolicSearchPlan plan(mz.size(), fixedStep, windows);
            const char *planStatus = "evaluated_exploratory";
            size_t candidateCount = 0;
            try
            {
                candidateCount = plan.evaluate({intensity.data()}).candidateIndices.size();
            }
            catch (const std::invalid_argument &)
            {
                const bool invalid = std::any_of(intensity.begin(), intensity.end(), [](float value)
                {
                    return !(value > 0 && std::isfinite(value));
                });
                if (!invalid)
                {
                    throw;
                }
                planStatus = "rejected_nonpositive_or_nonfinite";
            }
            const size_t zeros = std::count(intensity.begin(), intensity.end(), 0.0f);
            const size_t negatives = std::count_if(intensity.begin(), intensity.end(), [](float value) { return value < 0; });
            const size_t nonfinite = std::count_if(intensity.begin(), intensity.end(), [](float value) { return !std::isfinite(value); });
            out << "{\"type\":\"scan\",\"scan\":" << scan << ",\"status\":\"" << planStatus
                << "\",\"points\":" << mz.size() << ",\"plan_windows\":" << windows.size()
                << ",\"candidates\":" << candidateCount << ",\"zeros\":" << zeros
                << ",\"negatives\":" << negatives << ",\"nonfinite\":" << nonfinite << "}\n";
            for (size_t center : centers)
            {
                for (size_t scale : {4u, 8u})
                {
                    writeWindow(out, scan, "fixed_index", center, scale, mz, intensity, globalGap);
                }
            }
            for (size_t center : maxima)
            {
                for (size_t scale : {4u, 8u})
                {
                    writeWindow(out, scan, "observed_maximum", center, scale, mz, intensity, globalGap);
                }
            }
        }
        out.close();
        if (!out)
        {
            throw std::runtime_error("Failed to write complete audit output");
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
