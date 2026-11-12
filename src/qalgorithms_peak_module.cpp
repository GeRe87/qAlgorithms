#include "qalgorithms_peak_module.h"
#include "qalgorithms_peak_model_abi.h"
#include "detail/qalgorithms_peak_background.h"
#include "detail/qalgorithms_peak_cascade.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// Native module ownership and the bridge from the C ABI to the qPeaks
// validation stages. The loader does not interpret model coefficients.
namespace qAlgorithms
{
    namespace peak_model
    {
        struct NativePeakModel::Module
        {
#ifdef _WIN32
            HMODULE handle = nullptr;
#else
            void *handle = nullptr;
#endif
            const QAlgorithmsPeakModelV1 *model = nullptr;

            ~Module()
            {
                if (handle != nullptr)
                {
#ifdef _WIN32
                    FreeLibrary(handle);
#else
                    dlclose(handle);
#endif
                }
            }
        };

        namespace
        {
            PeakValidationStatus moduleStatus(int32_t status)
            {
                switch (status)
                {
                case QALGORITHMS_PEAK_PASSED: return PeakValidationStatus::passed;
                case QALGORITHMS_PEAK_FAILED: return PeakValidationStatus::failed;
                case QALGORITHMS_PEAK_UNAVAILABLE: return PeakValidationStatus::unavailable;
                default: throw std::runtime_error("Peak model module returned an unknown status");
                }
            }

            using Fit = std::unique_ptr<void, void (*)(void *)>;

            Fit fitWindow(const QAlgorithmsPeakModelV1 &model, const std::vector<double> &raw,
                double step, std::vector<double> &logs, PeakValidationStatus *status)
            {
                Fit fit(nullptr, model.destroy_fit);
                *status = PeakValidationStatus::unavailable;
                for (size_t i = 0; i < raw.size(); ++i)
                {
                    if (!(raw[i] > 0 && std::isfinite(raw[i])))
                    {
                        return fit;
                    }
                    logs[i] = std::log(raw[i]);
                }
                void *state = nullptr;
                const auto fitStatus = model.fit_window(logs.data(), raw.size(), step, &state);
                fit.reset(state);
                *status = moduleStatus(fitStatus);
                if (*status == PeakValidationStatus::passed && !fit)
                {
                    throw std::runtime_error("Peak model module reported success without a fit");
                }
                return fit;
            }

            PeakGeometryValidation checkModuleGeometry(const QAlgorithmsPeakModelV1 &model, const Fit &fit,
                size_t scale, double step, size_t minimum)
            {
                PeakGeometryValidation result;
                result.form = moduleStatus(model.has_peak_shape(fit.get()));
                if (result.form != PeakValidationStatus::passed)
                {
                    return result;
                }
                double apex = std::numeric_limits<double>::quiet_NaN();
                result.apexSupport = moduleStatus(model.apex(fit.get(), &apex));
                if (result.apexSupport == PeakValidationStatus::passed)
                {
                    return checkPeakGeometry(true, apex, scale, step, minimum);
                }
                return result;
            }

            bool predictLog(const QAlgorithmsPeakModelV1 &model, const Fit &fit,
                const std::vector<double> &grid, double step, std::vector<double> &prediction)
            {
                for (size_t i = 0; i < grid.size(); ++i)
                {
                    prediction[i] = std::numeric_limits<double>::quiet_NaN();
                    if (moduleStatus(model.predict_log(fit.get(), grid[i] * step, &prediction[i])) != PeakValidationStatus::passed
                        || !std::isfinite(prediction[i]))
                    {
                        return false;
                    }
                }
                return true;
            }

            void normalizeTemplate(std::vector<double> &prediction)
            {
                const double maximum = *std::max_element(prediction.begin(), prediction.end());
                for (double &value : prediction)
                {
                    value = std::exp(value - maximum);
                }
            }
        }

        NativePeakModel::NativePeakModel(const std::filesystem::path &path) : module_(std::make_unique<Module>())
        {
            const auto absolute = std::filesystem::absolute(path);
#ifdef _WIN32
            module_->handle = LoadLibraryExW(absolute.c_str(), nullptr,
                LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
            const auto symbol = module_->handle == nullptr ? nullptr
                : GetProcAddress(module_->handle, QALGORITHMS_PEAK_MODEL_ENTRY_V1);
            QAlgorithmsPeakModelEntryV1 entry = nullptr;
            static_assert(sizeof(entry) == sizeof(symbol), "Windows module function pointer size");
            std::memcpy(&entry, &symbol, sizeof(entry));
#else
            module_->handle = dlopen(absolute.c_str(), RTLD_NOW | RTLD_LOCAL);
            const auto entry = module_->handle == nullptr ? nullptr
                : reinterpret_cast<QAlgorithmsPeakModelEntryV1>(
                    dlsym(module_->handle, QALGORITHMS_PEAK_MODEL_ENTRY_V1));
#endif
            if (entry == nullptr)
            {
                throw std::runtime_error("Cannot load the peak model module or find its v1 entry point: " + absolute.string());
            }
            const auto *model = entry();
            if (model == nullptr || model->abi_version != QALGORITHMS_PEAK_MODEL_ABI_V1
                || model->struct_size < sizeof(QAlgorithmsPeakModelV1) || model->name == nullptr || *model->name == '\0'
                || model->parameter_count == 0 || !model->fit_window || !model->destroy_fit || !model->has_peak_shape
                || !model->apex || !model->predict_log || !model->form_confidence || !model->area_evidence
                || !model->window_coverage)
            {
                throw std::runtime_error("Invalid v1 peak model module descriptor: " + absolute.string());
            }
            module_->model = model;
        }

        NativePeakModel::~NativePeakModel() = default;
        NativePeakModel::NativePeakModel(NativePeakModel &&) noexcept = default;
        NativePeakModel &NativePeakModel::operator=(NativePeakModel &&) noexcept = default;

        std::string NativePeakModel::name() const
        {
            if (!module_)
            {
                throw std::logic_error("Peak model module has been moved");
            }
            return module_->model->name;
        }

        FullPeakCandidateValidation NativePeakModel::validateCandidate(const float *intensities, size_t length,
            double step, const FullPeakCandidateOptions &options) const
        {
            if (!module_)
            {
                throw std::logic_error("Peak model module has been moved");
            }
            const auto &model = *module_->model;
            const size_t scale = length / 2, minimum = options.background.minimumSamplesPerSide;
            const double edge = static_cast<double>(scale) * step;
            if (intensities == nullptr || length < 5 || length % 2 == 0 || length <= model.parameter_count
                || !(step > 0 && std::isfinite(edge)) || options.background.rotations == 0
                || !(options.background.alpha > 0 && options.background.alpha < 1)
                || !(options.logAlpha > 0 && options.logAlpha < .5)
                || !(options.minimumArea > 0 && std::isfinite(options.minimumArea))
                || !(options.minimumCoverage > 0 && options.minimumCoverage < 1))
            {
                throw std::invalid_argument("Invalid native peak candidate request");
            }
            const std::vector<double> raw(intensities, intensities + length);
            if (!std::all_of(raw.begin(), raw.end(), [](double value) { return std::isfinite(value); }))
            {
                throw std::invalid_argument("Native peak validation requires finite intensities");
            }
            std::vector<double> grid(length), logs(length), prediction(length);
            std::iota(grid.begin(), grid.end(), -static_cast<double>(scale));
            FullPeakCandidateValidation result;
            auto fit = fitWindow(model, raw, step, logs, &result.status);
            result.stages.geometry.form = result.status;
            if (result.status != PeakValidationStatus::passed)
            {
                return result;
            }
            const auto geometry = checkModuleGeometry(model, fit, scale, step, minimum);
            const size_t df = length - model.parameter_count;
            QAlgorithmsPeakNoiseV1 noise{0, df, options.logAlpha};
            result.stages = runPeakCascade(geometry, [&]()
            {
                if (!predictLog(model, fit, grid, step, prediction))
                {
                    return PeakValidationStatus::unavailable;
                }
                double rss = 0;
                for (size_t i = 0; i < length; ++i)
                {
                    const double residual = logs[i] - prediction[i];
                    rss += residual * residual;
                }
                noise.log_mse = rss / static_cast<double>(df);
                result.logNoise = {noise.log_mse, df, noise.alpha};
                if (!(noise.log_mse > 0 && std::isfinite(noise.log_mse)))
                {
                    return PeakValidationStatus::unavailable;
                }
                return moduleStatus(model.form_confidence(fit.get(), &noise));
            }, options.minimumArea, [&](double threshold)
            {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                QAlgorithmsPeakAreaV1 area{0, nan, nan, nan};
                const auto status = moduleStatus(model.area_evidence(fit.get(), &noise, threshold, &area));
                if (status == PeakValidationStatus::passed && (area.shape_certified != 1
                    || !(std::isfinite(area.log_lower_bound) && area.log_lower_bound > std::log(threshold))
                    || !(area.left >= -edge && area.left < area.right && area.right <= edge)))
                {
                    throw std::runtime_error("Peak model module returned an invalid passing area bound");
                }
                return PeakAreaEvidence{status, area.shape_certified == 1, area.log_lower_bound, area.left, area.right};
            }, options.minimumCoverage, [&](double threshold)
            {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                QAlgorithmsPeakCoverageV1 coverage{nan, nan, nan};
                const auto status = moduleStatus(model.window_coverage(fit.get(), &noise, -edge, edge, threshold, &coverage));
                if (status == PeakValidationStatus::passed && (!(coverage.fraction >= threshold && coverage.fraction <= 1)
                    || !(coverage.lower_confidence_bound >= threshold && coverage.lower_confidence_bound <= coverage.fraction)
                    || !(coverage.standard_error >= 0 && std::isfinite(coverage.standard_error))))
                {
                    throw std::runtime_error("Peak model module returned an invalid passing coverage bound");
                }
                return PeakWindowCoverage{status, coverage.fraction, coverage.standard_error, coverage.lower_confidence_bound};
            });
            result.status = peakCascadeStatus(result.stages);
            if (result.status != PeakValidationStatus::passed)
            {
                return result;
            }
            normalizeTemplate(prediction);
            fit.reset();
            result.background = runPeakBackground(raw, grid, prediction, true,
                [&](const auto &draw, auto &shape)
                {
                    PeakValidationStatus status;
                    const auto replay = fitWindow(model, draw, step, logs, &status);
                    if (status != PeakValidationStatus::passed
                        || checkModuleGeometry(model, replay, scale, step, minimum).apexSupport != PeakValidationStatus::passed
                        || !predictLog(model, replay, grid, step, shape))
                    {
                        return false;
                    }
                    normalizeTemplate(shape);
                    return true;
                }, options.background);
            result.status = result.background.status;
            return result;
        }
    }
}
