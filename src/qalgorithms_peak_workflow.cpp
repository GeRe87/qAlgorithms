#include "qalgorithms_peak_workflow.h"
#include "qalgorithms_peak_module.h"
#include "qalgorithms_peak_search.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

// mzML adapter for the declared qPeaks search. This file checks the
// measurement grid and writes scan-level results with joint covariance.
namespace qAlgorithms
{
    using namespace peak_model;

    namespace
    {
        std::optional<AffineSampleAxis> affineMzGrid(const std::vector<float> &mz)
        {
            const double origin = mz.front();
            const double step = (static_cast<double>(mz.back()) - origin) / static_cast<double>(mz.size() - 1);
            if (!(origin > 0 && step > 0 && std::isfinite(step)))
            {
                return std::nullopt;
            }
            for (size_t i = 0; i < mz.size(); ++i)
            {
                if (!std::isfinite(mz[i]) || (i && mz[i] <= mz[i - 1]))
                {
                    return std::nullopt;
                }
                // Allow storage rounding of an affine grid decoded as float.
                const double ulp = static_cast<double>(mz[i]) - std::nextafter(mz[i], 0.0f);
                if (std::abs(std::fma(static_cast<double>(i), step, origin) - mz[i]) > ulp)
                {
                    return std::nullopt;
                }
            }
            return AffineSampleAxis{origin, step, AxisUnit::MassToCharge};
        }

        const char *selectionStatus(ApexSelectionStatus status)
        {
            switch (status)
            {
            case ApexSelectionStatus::selected: return "measurementUnavailable";
            case ApexSelectionStatus::unsupported: return "unsupportedMaximum";
            case ApexSelectionStatus::insufficientCommonSupport: return "insufficientCommonSupport";
            case ApexSelectionStatus::ambiguous: return "ambiguous";
            case ApexSelectionStatus::numericalFailure: return "numericalFailure";
            }
            return "unavailable";
        }

        void writeModel(std::ostream &stream, const std::string &model)
        {
            stream << '"';
            for (char c : model)
            {
                if (c == '"')
                {
                    stream << '"';
                }
                stream << c;
            }
            stream << '"';
        }
    }

    void processCalibratedSpectra(const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<float> &retentionTimes, const std::string &filename, const UserInputSettings &args)
    {
        const auto &settings = args.calibratedPeaks.value();
        const auto output = std::filesystem::path(args.outputPath) / (filename + "_calibrated_peaks.csv");
        if (args.noOverwrite && std::filesystem::exists(output))
        {
            return;
        }
        if (indices.size() != retentionTimes.size())
        {
            throw std::invalid_argument("Calibrated peak workflow requires one retention time per scan");
        }
        std::shared_ptr<const NativePeakModel> module;
        if (settings.model.rfind("module:", 0) == 0)
        {
            module = std::make_shared<const NativePeakModel>(std::filesystem::path(settings.model.substr(7)));
            const size_t parameters = module->parameterCount();
            if (settings.minScale < parameters / 2 + parameters % 2)
            {
                throw std::invalid_argument("Calibrated native peak minimum scale leaves no residual degrees of freedom");
            }
        }
        const auto kind = settings.model == "quadratic" ? PeakModelKind::PiecewiseQuadratic : PeakModelKind::Hyperbolic;
        const std::string modelName = module ? "native:" + module->name() : "builtin:" + settings.model;
        auto staging = output;
        staging += ".partial";
        if (!std::filesystem::create_directory(staging))
        {
            throw std::runtime_error("Calibrated peak partial output already exists: " + staging.string());
        }
        const auto partial = staging / "peaks.csv";
        std::ofstream csv;
        try
        {
            csv.open(partial);
            if (!csv)
            {
                throw std::runtime_error("Could not open calibrated peak CSV: " + partial.string());
            }
            csv << "spectrumIndex,retentionTimeSeconds,model,record,status,familySize,backgroundAlphaPerWindow,"
                   "candidateCount,maximumFirst,maximumLast,windowCenter,windowHalfWidth,mz,areaIntensityMz,"
                   "varianceMz,covarianceMzArea,varianceArea\n";
            csv << std::setprecision(std::numeric_limits<double>::max_digits10);
            std::vector<float> mz, raw;
            for (size_t i = 0; i < indices.size(); ++i)
            {
                auto prefix = [&](const char *record, const char *status, size_t size, double alpha, size_t candidates)
                {
                    csv << indices[i] << ',' << retentionTimes[i] << ',';
                    writeModel(csv, modelName);
                    csv << ',' << record << ',' << status << ',' << size << ',' << alpha << ',' << candidates;
                };
                if (get_spectrum(&file, &mz, &raw, indices[i]) != 0 || mz.size() != raw.size())
                {
                    throw std::runtime_error("Could not decode calibrated peak spectrum " + std::to_string(indices[i]));
                }
                if (raw.size() < 5 || settings.minScale > (raw.size() - 1) / 2)
                {
                    prefix("scan", "noWindows", 0, 0, 0);
                    csv << ",,,,,,,,,\n";
                    continue;
                }
                const auto axis = affineMzGrid(mz);
                if (!axis)
                {
                    prefix("scan", "unsupportedAxis", 0, 0, 0);
                    csv << ",,,,,,,,,\n";
                    continue;
                }
                const size_t maximum = std::min(settings.maxScale, (raw.size() - 1) / 2);
                auto windows = makeUniformPeakSearchWindows(raw.size(), settings.minScale, maximum);
                // Convert the physical minimum using the same Jacobian as
                // the measured area. Window enumeration uses geometry only.
                const double minimumArea = static_cast<double>(
                    static_cast<long double>(settings.minimumAreaMz) * settings.step / axis->unitsPerSample);
                FullPeakCandidateOptions options{minimumArea, settings.minimumCoverage, .05, {}};
                options.background.rotations = settings.rotations;
                std::vector<CandidateModelWindows> builtins;
                std::vector<NativeCandidateModelWindows> native;
                if (module)
                {
                    native.push_back({module, settings.step, options, std::move(windows)});
                }
                else
                {
                    builtins.push_back({kind, settings.step, options, std::move(windows)});
                }
                const auto result = [&]()
                {
                    try
                    {
                        const CandidateSearchPlan plan(raw.size(), std::move(builtins), settings.alpha, std::move(native));
                        return findPeaks(plan, raw.data(), raw.size(), &*axis);
                    }
                    catch (const std::exception &error)
                    {
                        throw std::runtime_error("Calibrated peak spectrum " + std::to_string(indices[i]) + ": " + error.what());
                    }
                }();
                const size_t candidates = result.candidateIndices.size();
                prefix("scan", "evaluated", result.familySize, result.backgroundAlphaPerWindow, candidates);
                csv << ",,,,,,,,,\n";
                for (const auto &peak : result.measurements->peaks)
                {
                    const auto &selection = peak.selection;
                    prefix("peak", peak.properties ? "measured" : selectionStatus(selection.status),
                        result.familySize, result.backgroundAlphaPerWindow, candidates);
                    csv << ',' << selection.first << ',' << selection.last;
                    if (selection.status == ApexSelectionStatus::selected)
                    {
                        const auto &window = result.windows[selection.windowIndex].window;
                        csv << ',' << window.center << ',' << window.scale;
                    }
                    else
                    {
                        csv << ",,";
                    }
                    if (peak.properties)
                    {
                        const auto &p = *peak.properties;
                        csv << ',' << p.position << ',' << p.area << ',' << p.covariance[0]
                            << ',' << p.covariance[1] << ',' << p.covariance[3] << '\n';
                    }
                    else
                    {
                        csv << ",,,,,\n";
                    }
                }
                for (size_t index : result.measurements->unassignedWindowIndices)
                {
                    const auto &window = result.windows[index].window;
                    prefix("window", "unassigned", result.familySize, result.backgroundAlphaPerWindow, candidates);
                    csv << ",,," << window.center << ',' << window.scale << ",,,,,\n";
                }
            }
            csv.close();
            if (!csv)
            {
                throw std::runtime_error("Could not complete calibrated peak CSV: " + output.string());
            }
            if (!args.noOverwrite || !std::filesystem::exists(output))
            {
                std::filesystem::rename(partial, output);
            }
            else
            {
                std::filesystem::remove(partial);
            }
            std::filesystem::remove(staging);
        }
        catch (...)
        {
            csv.close();
            std::error_code ec;
            std::filesystem::remove(partial, ec);
            std::filesystem::remove(staging, ec);
            throw;
        }
    }
}
