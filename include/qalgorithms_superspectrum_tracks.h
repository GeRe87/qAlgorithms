#ifndef QALGORITHMS_SUPERSPECTRUM_TRACKS_H
#define QALGORITHMS_SUPERSPECTRUM_TRACKS_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <map>
#include <optional>
#include <vector>

namespace qAlgorithms
{
    struct SuperspectrumScanHeight
    {
        uint32_t scan;
        float retentionTime;
        double height;
        size_t distanceToWindowCenter;
    };

    struct SuperspectrumObservation
    {
        double mz, lower, upper, step;
        std::vector<SuperspectrumScanHeight> heights;
    };

    struct SuperspectrumTrackHeight
    {
        SuperspectrumScanHeight measurement;
        uint32_t firstScan;
    };

    using SuperspectrumTrackVisitor = std::function<void(size_t, double, std::optional<double>,
        const std::map<uint32_t, SuperspectrumTrackHeight> &)>;

    // Links overlapping windows. Writes each finished track as one row per scan.
    class SuperspectrumTracks
    {
    public:
        SuperspectrumTracks(std::ostream &output, size_t windowScans, SuperspectrumTrackVisitor visit = {});
        std::vector<size_t> append(const std::vector<SuperspectrumObservation> &observations, uint32_t firstScan);
        void finish();

    private:
        struct Track
        {
            size_t id = 0, windows = 0, lastWindow = 0;
            double mean = 0, squares = 0;
            SuperspectrumObservation last;
            std::map<uint32_t, SuperspectrumTrackHeight> heights;
        };

        void write(const Track &track);
        std::ostream &output;
        size_t windowScans;
        SuperspectrumTrackVisitor visit;
        std::vector<Track> active;
        size_t nextId = 1, window = 0;
    };
}

#endif
