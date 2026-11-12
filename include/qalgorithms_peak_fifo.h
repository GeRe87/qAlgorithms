#ifndef QALGORITHMS_PEAK_FIFO_H
#define QALGORITHMS_PEAK_FIFO_H

#include "qalgorithms_datatypes.h"
#include "qalgorithms_read_file.h"

#include <array>
#include <cstddef>
#include <optional>

namespace qAlgorithms
{
    // One fixed m/z segment. Scan IDs count all selected MS1 scans, including gaps.
    class PeakSegmentFifo
    {
    public:
        PeakSegmentFifo(double lower, double upper, size_t windowScans);
        void append(size_t scan, std::vector<double> mz, std::vector<double> intensity);
        void advance(size_t scan);
        std::optional<SharedMassFit> fitWindow() const;
        std::optional<SharedMassFit> fitAll() const;

    private:
        struct Scan
        {
            size_t index;
            std::vector<double> mz, intensity;
            // Centered cross-products of x, x² and log intensity.
            std::array<long double, 6> sums{};
        };

        struct Totals
        {
            std::array<long double, 6> sums{};
            size_t scans = 0, points = 0;
        };

        static void update(Totals &totals, const Scan &scan, bool add);
        std::optional<SharedMassFit> fit(const Totals &totals, size_t firstScan) const;

        double lower, upper, origin, scale;
        size_t windowScans, firstActive = 0;
        std::optional<size_t> lastScan;
        std::vector<Scan> scans;
        Totals window, expired;
    };

    std::vector<std::optional<SharedMassFit>> refineFeatureMasses(
        const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<CentroidPeak> &centroids, const std::vector<EIC> &bins,
        const std::vector<FeaturePeak> &features);
}

#endif
