#include "qalgorithms_peak_fifo.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Recover profile traces for the measured scans in each resolved feature.
namespace qAlgorithms
{
    std::vector<std::optional<SharedMassFit>> refineFeatureMasses(
        const XML_File &file, const std::vector<uint32_t> &indices,
        const std::vector<CentroidPeak> &centroids, const std::vector<EIC> &bins,
        const std::vector<FeaturePeak> &features)
    {
        constexpr size_t windowScans = 5;
        struct Segment
        {
            double lower = 0, upper = 0;
            size_t lastScan = 0;
            std::optional<PeakSegmentFifo> fifo;
        };
        std::vector<std::optional<SharedMassFit>> fits(features.size());
        std::vector<Segment> segments(features.size());
        std::vector<std::vector<size_t>> byScan(indices.size());
        for (size_t i = 0; i < features.size(); ++i)
        {
            const auto &feature = features[i];
            if (!feature.unresolvedMaximaRT.empty())
                continue;

            const auto &bin = bins.at(feature.eic_ID);
            double lower = 0, upper = std::numeric_limits<double>::infinity();
            std::vector<size_t> scans;
            for (size_t j = 0; j < feature.length; ++j)
            {
                const size_t point = static_cast<size_t>(feature.idxPeakStart) + j;
                if (point && bin.df.at(point) == bin.df.at(point - 1))
                    continue;
                const auto &centroid = centroids.at(bin.cenID.at(point));
                const size_t scan = bin.scanNumbers.at(point);
                if (centroid.number_MS1 != scan || (!scans.empty() && scan <= scans.back()) ||
                    !(centroid.width > 0 && std::isfinite(centroid.width)))
                {
                    scans.clear();
                    break;
                }
                const double mz = centroid.mz;
                lower = std::max(lower, mz - centroid.width);
                upper = std::min(upper, mz + centroid.width);
                scans.push_back(scan);
            }
            if (scans.size() < 2 || !(lower > 0 && upper > lower))
                continue;

            // Fix the segment before accumulation so an evicted scan can
            // subtract exactly the same cross-products it added.
            segments[i] = {lower, upper, scans.back(), std::nullopt};
            for (size_t scan : scans)
                byScan.at(scan).push_back(i);
        }

        std::vector<double> mz, raw;
        for (size_t scan = 0; scan < indices.size(); ++scan)
        {
            if (byScan[scan].empty())
                continue;
            if (get_spectrum(&file, &mz, &raw, indices[scan]) != 0 || mz.size() != raw.size())
                throw std::runtime_error("Could not decode FIFO profile spectrum " + std::to_string(indices[scan]));
            const bool ordered = std::is_sorted(mz.begin(), mz.end());
            for (size_t i : byScan[scan])
            {
                auto &segment = segments[i];
                if (!segment.fifo)
                    segment.fifo.emplace(segment.lower, segment.upper, windowScans);
                auto &fifo = *segment.fifo;
                fifo.advance(scan);
                if (ordered)
                {
                    const auto first = std::lower_bound(mz.begin(), mz.end(), segment.lower);
                    const auto last = std::upper_bound(first, mz.end(), segment.upper);
                    const auto rawFirst = raw.begin() + (first - mz.begin());
                    const auto rawLast = rawFirst + (last - first);
                    if (last - first >= 3 && std::adjacent_find(first, last) == last &&
                        std::all_of(rawFirst, rawLast, [](double value) { return value > 0; }))
                    {
                        fifo.append(scan, {first, last}, {rawFirst, rawLast});
                    }
                }
                if (scan == segment.lastScan)
                {
                    fits[i] = fifo.fitAll();
                    segment.fifo.reset();
                }
            }
        }
        return fits;
    }
}
