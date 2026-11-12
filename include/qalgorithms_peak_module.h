#ifndef QALGORITHMS_PEAK_MODULE_H
#define QALGORITHMS_PEAK_MODULE_H

#include "qalgorithms_peak_validation.h"

#include <filesystem>
#include <memory>
#include <string>

namespace qAlgorithms
{
    namespace peak_model
    {
        // Native model loaded by explicit path. The ABI is defined separately
        // in qalgorithms_peak_model_abi.h and does not expose C++ types.
        // Concurrent calls require a module that supports them.
        class NativePeakModel
        {
        public:
            explicit NativePeakModel(const std::filesystem::path &path);
            ~NativePeakModel();
            NativePeakModel(const NativePeakModel &) noexcept;
            NativePeakModel &operator=(const NativePeakModel &) noexcept;
            NativePeakModel(NativePeakModel &&) noexcept;
            NativePeakModel &operator=(NativePeakModel &&) noexcept;

            std::string name() const;
            size_t parameterCount() const;
            FullPeakCandidateValidation validateCandidate(const float *intensities, size_t length,
                double step, const FullPeakCandidateOptions &options, CandidateFitOutput *output = nullptr) const;

        private:
            struct Module;
            std::shared_ptr<const Module> module_;
        };
    }
}

#endif
