#include "qalgorithms_input_output.h"
#include "qalgorithms_qbin.h"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_read_file.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace fs = std::filesystem;
using namespace qAlgorithms;

static void exportProfile(const XML_File &file, const fs::path &out, Polarities polarity,
                          PeakFitOptions fit)
{
    const auto scans = filter_spectra(&file, true, polarity, true);
    if (scans.empty())
        return;
    fs::create_directories(out);
    std::vector<float> rt;
    get_spectra_RT(&file, &scans, &rt);
    if (rt.size() != scans.size())
        throw std::runtime_error("Could not read retention times.");

    std::vector<CentroidPeak> centroids;
    std::vector<Span_i32> spans;
    findCentroids(&file, &scans, &centroids, fit, &spans);
    const auto bins = centroids.size() >= 5 ? performQbinning(&centroids, &rt) : std::vector<EIC>{};
    std::vector<FeaturePeak> features;
    findFeatures(&bins, &rt, &features, nullptr, fit);
    if (!printFeatureList(&features, out, "profile", &bins, true, false))
        throw std::runtime_error("Could not write feature statistics.");

    std::ofstream scanCsv(out / "scans.csv"), centroidCsv(out / "centroids.csv");
    std::ofstream binCsv(out / "bins.csv"), memberCsv(out / "feature_centroids.csv");
    std::ofstream backgroundCsv(out / "background.csv");
    std::ofstream relativeCsv(out / "relative_noise.csv");
    scanCsv << std::setprecision(12) << "scan,rt_s,points\n";
    centroidCsv << std::setprecision(12)
                << "id,scan,start,end,mz,area,height,fwhm,mz_se,area_se,height_se,mz_step,bin,right_left_halfwidth\n";
    binCsv << std::setprecision(12) << "bin,centroid,mz,rt_s,mz_prediction_interval\n";
    memberCsv << "feature,centroid\n";
    backgroundCsv << std::setprecision(12) << "scan,mz,intensity\n";
    relativeCsv << std::setprecision(12) << "log_residual_proxy\n";

    std::vector<int64_t> binIds(centroids.size(), -1);
    for (size_t i = 0; i < bins.size(); ++i)
    {
        const auto &bin = bins[i];
        for (size_t j = 0; j < bin.cenID.size(); ++j)
        {
            if (j && bin.df[j] == bin.df[j - 1])
                continue;
            binIds.at(bin.cenID[j]) = static_cast<int64_t>(i);
            binCsv << i << ',' << bin.cenID[j] << ',' << bin.mz[j] << ',' << bin.RT[j]
                   << ',' << bin.predInterval[j] << '\n';
        }
    }
    for (size_t i = 0; i < features.size(); ++i)
    {
        const auto &feature = features[i];
        if (!feature.unresolvedMaximaRT.empty())
            continue;
        const auto &bin = bins.at(feature.eic_ID);
        for (size_t j = feature.idxPeakStart; j < size_t(feature.idxPeakStart) + feature.length; ++j)
            if (j == 0 || bin.df.at(j) != bin.df.at(j - 1))
                memberCsv << i + 1 << ',' << bin.cenID.at(j) << '\n';
    }

    size_t first = 0;
    std::vector<double> mzValues, intensities;
    for (size_t i = 0; i < scans.size(); ++i)
    {
        if (get_spectrum(&file, &mzValues, &intensities, scans[i]) != 0)
            throw std::runtime_error("Could not reread an MS1 profile scan.");
        scanCsv << scans[i] << ',' << rt[i] << ',' << mzValues.size() << '\n';
        size_t last = first;
        std::vector<bool> covered(mzValues.size(), false);
        while (last < centroids.size() && centroids[last].number_MS1 == i)
        {
            const auto &peak = centroids[last];
            const auto span = spans.at(last);
            if (span.startIdx < 0 || span.length < 2 || size_t(span.endIdx()) >= mzValues.size())
                throw std::runtime_error("Invalid profile span returned by centroiding.");
            const double step = (mzValues[span.endIdx()] - mzValues[span.startIdx]) / (span.length - 1);
            std::fill(covered.begin() + span.startIdx, covered.begin() + span.endIdx() + 1, true);
            const auto apex = std::max_element(intensities.begin() + span.startIdx,
                                               intensities.begin() + span.endIdx() + 1);
            const double half = *apex / 2;
            auto left = apex, right = apex;
            while (left != intensities.begin() + span.startIdx && *left > half)
                --left;
            while (right != intensities.begin() + span.endIdx() && *right > half)
                ++right;
            centroidCsv << last << ',' << scans[i] << ',' << span.startIdx << ',' << span.endIdx() + 1
                        << ',' << peak.mz << ',' << peak.area << ',' << peak.height << ',' << peak.width
                        << ',' << peak.mzUncertainty << ',' << peak.areaUncertainty << ',' << peak.heightUncertainty
                        << ',' << step << ',' << binIds[last] << ',';
            if (*left <= half && *right <= half && left != apex && right != apex)
            {
                const double leftWidth = mzValues[apex - intensities.begin()] - mzValues[left - intensities.begin()];
                if (leftWidth > 0)
                    centroidCsv << (mzValues[right - intensities.begin()] - mzValues[apex - intensities.begin()]) / leftWidth;
            }
            centroidCsv << '\n';
            const auto center = apex - intensities.begin();
            if (center > span.startIdx && center + 2 <= span.endIdx() &&
                intensities[center - 1] > half && intensities[center + 1] > half &&
                intensities[center + 2] > half && half > 0)
                relativeCsv << (std::log(intensities[center + 2]) - 3 * std::log(intensities[center + 1]) +
                    3 * std::log(intensities[center]) - std::log(intensities[center - 1])) / std::sqrt(20.) << '\n';
            ++last;
        }
        first = last;
        const size_t stride = std::max(size_t(1), mzValues.size() / 64);
        for (size_t j = i % stride; j < mzValues.size(); j += stride)
            if (!covered[j] && std::isfinite(intensities[j]) && intensities[j] >= 0)
                backgroundCsv << scans[i] << ',' << mzValues[j] << ',' << intensities[j] << '\n';
    }
    for (auto *stream : {&scanCsv, &centroidCsv, &binCsv, &memberCsv, &backgroundCsv, &relativeCsv})
    {
        stream->close();
        if (!*stream)
            throw std::runtime_error("Could not finish writing the calibration tables.");
    }
}

static int run(const fs::path &input, const fs::path &out, std::string_view model, bool benchmarkMode)
{
    PeakFitOptions fit;
    fit.benchmarkMode = benchmarkMode;
    if (model == "hyperbolic")
        fit.model = PeakModelKind::Hyperbolic;
    else if (model != "quadratic")
        throw std::invalid_argument("Choose quadratic or hyperbolic for the calibration.");
    XML_File file(input.c_str(), mzML);
    if (file.defective)
        throw std::runtime_error("Could not open the mzML file.");
    exportProfile(file, out / "positive", Polarities::positive, fit);
    exportProfile(file, out / "negative", Polarities::negative, fit);
    if (!fs::exists(out / "positive" / "scans.csv") && !fs::exists(out / "negative" / "scans.csv"))
        throw std::runtime_error("No supported MS1 profile scans were found.");
    return 0;
}

#ifdef _WIN32
int wmain(int argc, wchar_t **argv)
#else
int main(int argc, char **argv)
#endif
{
    try
    {
        if (argc != 4 && (argc != 5 || fs::path(argv[4]) != "--benchmark-mode"))
            throw std::invalid_argument("Usage: qmzml-profile INPUT.mzML OUTPUT_DIRECTORY KERNEL [--benchmark-mode]");
        return run(argv[1], argv[2], fs::path(argv[3]).string(), argc == 5);
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
