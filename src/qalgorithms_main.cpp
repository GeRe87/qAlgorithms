// Command-line workflow: mzML profile spectra -> centroids -> bins -> feature CSV.
#include "qalgorithms_datatypes.h"
#include "qalgorithms_input_output.h"
#include "qalgorithms_qbin.h"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_peak_module.h"
#include "qalgorithms_peak_search.h"
#include "qalgorithms_read_file.h"

#include <cstdio>
#include <cctype>
#include <charconv>
#include <ctime>
#include <exception>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using namespace qAlgorithms;

static std::string outputKey(std::string name)
{
#ifdef _WIN32
    for (char &c : name)
        c = (char)std::tolower((unsigned char)c);
#endif
    return name;
}

// Research output only: one grouped linear fit per MS1 profile spectrum.
// Does not modify the existing centroid/bin/feature processing path.
static void processFoldingSpectra(
    const XML_File &inputFile, const std::vector<uint32_t> &selectedIndices,
    const std::vector<float> &retentionTimes,
    const std::filesystem::path &outputDirectory, const std::string &filename,
    const UserInputSettings &args)
{
    const auto output = outputDirectory / (filename + "_folding.csv");
    if (std::filesystem::exists(output) && args.noOverwrite)
    {
        if (!args.silent)
            printf("    Folding output exists, leaving unchanged: %s\n",
                   output.string().c_str());
        return;
    }

    const auto temp = std::filesystem::path(output.string() + ".partial");
    // Never overwrite another process's unfinished work or any existing
    // user-provided .partial file.
    if (std::filesystem::exists(temp))
        throw std::runtime_error("Folding partial output already exists: " + temp.string());

    std::ofstream stream(temp);
    if (!stream.is_open())
        throw std::runtime_error("Could not open temporary folding CSV: " + temp.string());

    const FoldingSamplingAxis axis = args.foldingAxis == FoldingCliAxis::tof
        ? FoldingSamplingAxis::tofSqrtMz : FoldingSamplingAxis::orbitrapInvSqrtMz;
    const FoldingCurvatureMode mode = args.foldingSymmetric
        ? FoldingCurvatureMode::symmetric : FoldingCurvatureMode::asymmetric;
    const char *const axisName = axis == FoldingSamplingAxis::tofSqrtMz
        ? "tof" : "orbitrap";
    const char *const modeName = mode == FoldingCurvatureMode::symmetric
        ? "symmetric" : "asymmetric";
    const char *const polarityName = inputFile.polarityMode == Polarities::mixed
        ? (spectrum_polarity(&inputFile, selectedIndices.front()) == Polarities::positive
                ? "positive" : "negative")
        : (inputFile.polarityMode == Polarities::positive ? "positive" : "negative");

    size_t fittedScans = 0, fittedPeaks = 0, rejectedCandidates = 0;
    std::vector<float> mz, intensities;
    stream << "spectrumIndex,retentionTimeSeconds,polarity,axis,model,"
              "samplingStep,profilesInFit,anchorIndex,anchorMz,anchorIntensity,"
              "samplesUsed,validApex,apexMz,apexUncertaintyMz,"
              "b0,b1,b2,b3,b0StandardError,b1StandardError,"
              "b2StandardError,b3StandardError,logRSS,pooledLogMSE,degreesOfFreedom\n";
    stream << std::setprecision(14);

    try
    {
        for (size_t i = 0; i < selectedIndices.size(); ++i)
        {
            if (get_spectrum(&inputFile, &mz, &intensities, selectedIndices[i]) != 0)
                throw std::runtime_error("Could not decode MS1 profile spectrum " +
                                         std::to_string(selectedIndices[i]));
            if (mz.size() != intensities.size())
                throw std::runtime_error("m/z and intensity arrays have different lengths");

            FoldingSpectrumFit folding;
            try
            {
                folding = foldSpectrumProfiles(mz.data(), intensities.data(),
                                               mz.size(), axis, mode,
                                               args.foldingMinIntensity);
            }
            catch (const std::exception &error)
            {
                throw std::runtime_error("Folding spectrum index " +
                    std::to_string(selectedIndices[i]) + ": " + error.what());
            }
            rejectedCandidates += folding.observedMaxima - folding.selectedCandidates;
            if (!folding.fitted) continue;
            ++fittedScans;
            const auto &fit = folding.fit;
            const double commonSE2 = std::sqrt(std::max(0.0,
                fit.curvatureCovarianceFactor[0] * fit.mse));
            const double commonSE3 = std::sqrt(std::max(0.0,
                fit.curvatureCovarianceFactor[3] * fit.mse));
            for (const FoldingSpectrumPeak &peak : folding.peaks)
            {
                const auto &b = peak.coefficients;
                stream << selectedIndices[i] << ',' << retentionTimes[i] << ','
                       << polarityName << ',' << axisName << ',' << modeName << ','
                       << folding.samplingStep << ',' << folding.peaks.size() << ','
                       << peak.anchorIndex << ',' << peak.anchorMz << ','
                       << peak.anchorIntensity << ',' << peak.pointCount << ','
                       << (peak.validApex ? 1 : 0) << ',';
                if (peak.validApex)
                    stream << peak.apexMz << ',' << peak.apexUncertaintyMz;
                else
                    stream << ',';
                stream << ',' << b.b0 << ',' << b.b1 << ','
                       << fit.b2 << ',' << fit.b3 << ','
                       << std::sqrt(std::max(0.0, b.varianceB0Factor * fit.mse)) << ','
                       << std::sqrt(std::max(0.0, b.varianceB1Factor * fit.mse)) << ','
                       << commonSE2 << ',' << commonSE3 << ','
                       << b.residualSumSquares << ',' << fit.mse << ','
                       << fit.degreesOfFreedom << '\n';
                ++fittedPeaks;
            }
        }

        stream.close();
        if (!stream)
            throw std::runtime_error("Could not complete folding output CSV");

        // Write the complete output before replacing an existing CSV. This
        // avoids publishing a partially analyzed file on decoding failure.
        if (std::filesystem::exists(output))
            std::filesystem::remove(output);
        std::filesystem::rename(temp, output);
    }
    catch (...)
    {
        stream.close();
        std::error_code ec;
        std::filesystem::remove(temp, ec);
        throw;
    }

    if (!args.silent)
        printf("    Experimental folding: %zu peaks in %zu of %zu MS1 scans "
               "(%zu candidates rejected); output: %s\n",
               fittedPeaks, fittedScans, selectedIndices.size(),
               rejectedCandidates, output.string().c_str());
    if (fittedScans == 0)
        (void)fprintf(stderr, "Warning: no eligible grouped profiles in %s. "
                              "Folding CSV contains only its header.\n", filename.c_str());
}

static void processFile(const std::filesystem::path &pathSource, const UserInputSettings &args,
                        std::set<std::string> &outputNames)
{
    if (args.printCentroidsMZML)
    {
        centroids_to_mzml(&pathSource);
        return;
    }

    XML_File inputFile(pathSource.c_str(), mzML);
    if (inputFile.defective)
        throw std::runtime_error("The mzML file could not be read.");

    std::vector<Polarities> polarities{inputFile.polarityMode};
    if (inputFile.polarityMode == Polarities::mixed)
        polarities = {Polarities::positive, Polarities::negative};

    size_t processedPolarities = 0;
    for (const Polarities polarity : polarities)
    {
        std::vector<uint32_t> selectedIndices = filter_spectra(&inputFile, true, polarity, true);
        if (selectedIndices.empty())
            continue;
        if (args.peakFolding && args.foldingMaxScans != 0 &&
            selectedIndices.size() > args.foldingMaxScans)
            selectedIndices.resize(args.foldingMaxScans);

        std::string filename = pathSource.stem().string();
        if (inputFile.polarityMode == Polarities::mixed)
            filename += polarity == Polarities::positive ? "_positive" : "_negative";
        if (!args.outputPath.empty() && !outputNames.insert(outputKey(filename)).second)
            throw std::runtime_error("Two inputs would write the same output filename. Process them into separate directories.");
        if (!args.silent)
            printf("    Processing %zu %s MS1 profile spectra\n", selectedIndices.size(),
                   polarity == Polarities::positive ? "positive" : "negative");

        std::vector<float> retentionTimes;
        get_spectra_RT(&inputFile, &selectedIndices, &retentionTimes);
        if (retentionTimes.size() != selectedIndices.size())
            throw std::runtime_error("Could not read retention times for the selected spectra.");

        if (args.peakFolding)
        {
            processFoldingSpectra(inputFile, selectedIndices, retentionTimes,
                                  args.outputPath, filename, args);
            ++processedPolarities;
            continue;
        }

        if (args.printProfileSection)
            printProfileSections(&inputFile, &args, &selectedIndices, filename);
        if (args.term == TerminateAfter::profilesec)
        {
            ++processedPolarities;
            continue;
        }

        std::vector<CentroidPeak> centroids;
        findCentroids(&inputFile, &selectedIndices, &centroids);
        if (!args.silent)
        {
            printf("    Centroiding complete: %zu centroids\n", centroids.size());
            (void)fflush(stdout);
        }
        if (args.printCentroids)
            printCentroids(&centroids, &retentionTimes, args.outputPath, filename, args.silent, args.noOverwrite);
        if (args.term == TerminateAfter::centroids)
        {
            ++processedPolarities;
            continue;
        }

        std::vector<EIC> bins;
        if (centroids.size() >= 5)
            bins = performQbinning(&centroids, &retentionTimes);
        if (!args.silent)
        {
            printf("    Binning complete: %zu extracted ion chromatograms\n", bins.size());
            (void)fflush(stdout);
        }
        if (args.printBins)
            printBins(&centroids, &bins, args.outputPath, filename, args.silent, args.noOverwrite);
        if (args.term == TerminateAfter::binning)
        {
            ++processedPolarities;
            continue;
        }

        std::vector<FeaturePeak> features;
        findFeatures(&bins, &retentionTimes, &features);
        size_t unresolvedCount = 0;
        for (const auto &feature : features) unresolvedCount += !feature.unresolvedMaximaRT.empty();
        if (unresolvedCount)
            (void)fprintf(stderr, "Warning: %s contains %zu explicitly marked unresolved signal descriptions; these are not resolved individual features.\n", filename.c_str(), unresolvedCount);
        if (features.empty())
            (void)fprintf(stderr, "Warning: no features detected in %s; the feature CSV contains its header only.\n",
                          filename.c_str());
        if (args.printFeatures &&
            !printFeatureList(&features, args.outputPath, filename, &bins, args.silent, args.noOverwrite))
            throw std::runtime_error("Could not write the feature CSV.");
        if (!args.silent)
            printf("    produced %zu centroids, %zu bins and %zu features\n",
                   centroids.size(), bins.size(), features.size());
        ++processedPolarities;
    }
    if (processedPolarities == 0)
        throw std::runtime_error("No MS1 profile spectra with a supported polarity were found. Centroid-only data is not supported.");
}

// Parse the explicit coordinates and thresholds used by native-model diagnostics.
static double parsePeakModuleNumber(const std::string &value)
{
    size_t end = 0;
    const double number = std::stod(value, &end);
    if (end != value.size() || !std::isfinite(number))
    {
        throw std::invalid_argument("Expected a finite numeric module argument");
    }
    return number;
}

static size_t parsePeakModuleCount(std::string_view value)
{
    size_t number = 0;
    const auto end = value.data() + value.size();
    const auto parsed = std::from_chars(value.data(), end, number);
    if (parsed.ec != std::errc{} || parsed.ptr != end)
    {
        throw std::invalid_argument("Expected a nonnegative integer module argument");
    }
    return number;
}

static std::vector<float> readPeakModuleSamples(int start, int argc, const char *argv[])
{
    std::vector<float> raw;
    for (int i = start; i < argc; ++i)
    {
        const double value = parsePeakModuleNumber(argv[i]);
        if (std::abs(value) > std::numeric_limits<float>::max())
        {
            throw std::invalid_argument("Intensity is outside the supported range");
        }
        raw.push_back(static_cast<float>(value));
    }
    return raw;
}

static const char *peakStatusName(peak_model::PeakValidationStatus value)
{
    switch (value)
    {
    case peak_model::PeakValidationStatus::passed: return "passed";
    case peak_model::PeakValidationStatus::failed: return "failed";
    case peak_model::PeakValidationStatus::notEvaluated: return "notEvaluated";
    case peak_model::PeakValidationStatus::unavailable: return "unavailable";
    }
    return "unavailable";
}

// Diagnostic entry point for a prespecified window and an external model.
static int runPeakModelModule(int argc, const char *argv[])
{
    if (argc < 11)
    {
        throw std::invalid_argument(
            "Usage: qAlgorithms --peak-model-module PATH STEP MIN_AREA MIN_COVERAGE I0 I1 I2 I3 I4 [...]");
    }
    const auto raw = readPeakModuleSamples(6, argc, argv);
    const peak_model::FullPeakCandidateOptions options{parsePeakModuleNumber(argv[4]), parsePeakModuleNumber(argv[5]), .05, {}};
    const peak_model::NativePeakModel model{std::filesystem::path(argv[2])};
    const auto result = model.validateCandidate(raw.data(), raw.size(), parsePeakModuleNumber(argv[3]), options);
    printf("model=%s\nstatus=%s\nform=%s\napex=%s\nformConfidence=%s\narea=%s\ncoverage=%s\nbackground=%s\n",
        model.name().c_str(), peakStatusName(result.status), peakStatusName(result.stages.geometry.form),
        peakStatusName(result.stages.geometry.apexSupport), peakStatusName(result.stages.formConfidence),
        peakStatusName(result.stages.areaEvidence.status), peakStatusName(result.stages.windowCoverage.status), peakStatusName(result.background.status));
    return 0;
}

// Keep descriptor bytes from colliding with the window record delimiters.
static std::string nativePeakModelLabel(const std::string &name)
{
    constexpr char hex[] = "0123456789ABCDEF";
    std::string label = "native:";
    for (unsigned char value : name)
    {
        if ((value >= 'a' && value <= 'z') || (value >= 'A' && value <= 'Z') || (value >= '0' && value <= '9')
            || value == '-' || value == '_' || value == '.')
        {
            label += static_cast<char>(value);
        }
        else
        {
            label += '%';
            label += hex[value >> 4];
            label += hex[value & 15];
        }
    }
    return label;
}

// Run a uniform declared family, optionally including the quadratic baseline.
static int runPeakModelModuleSearch(int argc, const char *argv[])
{
    using namespace peak_model;
    if (argc < 15)
    {
        throw std::invalid_argument(
            "Usage: qAlgorithms --peak-model-module-search PATH STEP MIN_AREA MIN_COVERAGE ALPHA ROTATIONS "
            "MIN_SCALE MAX_SCALE [--quadratic MIN_AREA MIN_COVERAGE] I0 I1 I2 I3 I4 [...]");
    }
    FullPeakCandidateOptions options{parsePeakModuleNumber(argv[4]), parsePeakModuleNumber(argv[5]), .05, {}};
    options.background.rotations = parsePeakModuleCount(argv[7]);
    const bool includeQuadratic = std::string_view(argv[10]) == "--quadratic";
    const int sampleStart = includeQuadratic ? 13 : 10;
    if (argc < sampleStart + 5)
    {
        throw std::invalid_argument("Native model search requires at least five intensities");
    }
    const auto raw = readPeakModuleSamples(sampleStart, argc, argv);
    const auto windows = makeUniformPeakSearchWindows(raw.size(), parsePeakModuleCount(argv[8]), parsePeakModuleCount(argv[9]));
    std::vector<CandidateModelWindows> builtins;
    if (includeQuadratic)
    {
        auto quadratic = options;
        quadratic.minimumArea = parsePeakModuleNumber(argv[11]);
        quadratic.minimumCoverage = parsePeakModuleNumber(argv[12]);
        builtins.push_back({PeakModelKind::PiecewiseQuadratic, 1, quadratic, windows});
    }
    auto module = std::make_shared<const NativePeakModel>(std::filesystem::path(argv[2]));
    const CandidateSearchPlan plan(raw.size(), std::move(builtins), parsePeakModuleNumber(argv[6]),
        {{std::move(module), parsePeakModuleNumber(argv[3]), options, windows}});
    const auto result = findPeaks(plan, raw.data(), raw.size());
    printf("familySize=%zu\nwindowAlpha=%.17g\ncandidates=%zu\n",
        result.familySize, result.backgroundAlphaPerWindow, result.candidateIndices.size());
    for (const auto &entry : result.windows)
    {
        printf("window=%s,%zu,%zu,%s,%s,%.17g\n",
            entry.model ? "builtin:piecewise-quadratic" : nativePeakModelLabel(entry.nativeModelName).c_str(),
            entry.window.center, entry.window.scale,
            peakStatusName(entry.validation.status), peakStatusName(entry.validation.background.status), entry.adjustedBackgroundP);
    }
    return 0;
}

int main(int argc, const char *argv[])
{
    try
    {
        if (argc > 1 && std::string_view(argv[1]) == "--peak-model-module")
        {
            return runPeakModelModule(argc, argv);
        }
        if (argc > 1 && std::string_view(argv[1]) == "--peak-model-module-search")
        {
            return runPeakModelModuleSearch(argc, argv);
        }
        UserInputSettings args = passCliArgs(argc, argv);
        if (!inputsAreSensible(args))
            return 1;
        const std::vector<std::filesystem::path> tasks = controlInput(&args);

        // A shared output directory must not silently merge distinct samples.
        std::set<std::string> outputStems;
        if (!args.outputPath.empty())
        {
            for (const auto &path : tasks)
                if (!outputStems.insert(outputKey(path.stem().string())).second)
                    throw std::runtime_error("Input files have the same filename stem. Process them into separate output directories.");
        }

        const clock_t start = clock();
        size_t failures = 0;
        std::set<std::string> outputNames;
        for (size_t i = 0; i < tasks.size(); ++i)
        {
            if (!args.silent)
                printf("\nReading file %zu of %zu: %s\n", i + 1, tasks.size(), tasks[i].string().c_str());
            try
            {
                processFile(tasks[i], args, outputNames);
            }
            catch (const std::exception &error)
            {
                (void)fprintf(stderr, "Error processing %s: %s\n", tasks[i].string().c_str(), error.what());
                ++failures;
                if (!args.skipError)
                    return 1;
            }
        }
        if (!args.silent)
            printf("Completed %zu of %zu files successfully (%zu failed) in %.3f s.\n",
                   tasks.size() - failures, tasks.size(), failures,
                   double(clock() - start) / CLOCKS_PER_SEC);
        return failures == 0 ? 0 : 1;
    }
    catch (const std::exception &error)
    {
        (void)fprintf(stderr, "Error: %s\n", error.what());
        return 1;
    }
}
