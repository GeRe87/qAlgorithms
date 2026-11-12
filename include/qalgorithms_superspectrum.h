#ifndef QALGORITHMS_SUPERSPECTRUM_H
#define QALGORITHMS_SUPERSPECTRUM_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace qAlgorithms
{
    struct XML_File;

    struct SuperspectrumSettings
    {
        double gridStep;
        double maxGap;
        size_t maxHalfWidth;
    };

    void processSuperspectra(const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<float> &retentionTimes, const std::filesystem::path &output,
        const SuperspectrumSettings *fixed, bool noOverwrite, bool silent);
}

#endif
