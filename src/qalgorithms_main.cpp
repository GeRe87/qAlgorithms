// Command-line workflow: mzML profile spectra -> centroids -> bins -> feature CSV.
#include "qalgorithms_datatypes.h"
#include "qalgorithms_input_output.h"
#include "qalgorithms_qbin.h"
#include "qalgorithms_qpeaks.h"
#include "qalgorithms_read_file.h"

#include <cstdio>
#include <cctype>
#include <ctime>
#include <exception>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
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
        const std::vector<uint32_t> selectedIndices = filter_spectra(&inputFile, true, polarity, true);
        if (selectedIndices.empty())
            continue;

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

int main(int argc, const char *argv[])
{
    try
    {
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
