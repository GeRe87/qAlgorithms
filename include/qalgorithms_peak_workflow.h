#ifndef QALGORITHMS_PEAK_WORKFLOW_H
#define QALGORITHMS_PEAK_WORKFLOW_H

#include "qalgorithms_input_output.h"
#include "qalgorithms_read_file.h"

namespace qAlgorithms
{
    // Write calibrated scan-level peaks for one polarity of an mzML input.
    void processCalibratedSpectra(const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<float> &retentionTimes, const std::string &filename, const UserInputSettings &args);
}

#endif
