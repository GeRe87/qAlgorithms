#include "qalgorithms_superspectrum_tracks.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <tuple>
#include <utility>

namespace qAlgorithms
{
    SuperspectrumTracks::SuperspectrumTracks(std::ostream &output, size_t windowScans,
        SuperspectrumTrackVisitor visit)
        : output(output), windowScans(windowScans), visit(std::move(visit))
    {
        output << "track,windows,mzMean,mzStdDev,scan,retentionTimeSeconds,height,heightWindowFirstScan\n";
    }

    void SuperspectrumTracks::write(const Track &track)
    {
        std::optional<double> spread;
        if (track.windows > 1)
            spread = std::sqrt(std::max(0.0, track.squares) / (track.windows - 1));
        for (const auto &[scan, height] : track.heights)
        {
            output << track.id << ',' << track.windows << ',' << track.mean << ',';
            if (spread)
                output << *spread;
            output << ',' << scan << ',' << height.measurement.retentionTime << ','
                   << height.measurement.height << ',' << height.firstScan << '\n';
        }
        if (visit)
            visit(track.id, track.mean, spread, track.heights);
    }

    std::vector<size_t> SuperspectrumTracks::append(
        const std::vector<SuperspectrumObservation> &observations, uint32_t firstScan)
    {
        ++window;
        struct Match
        {
            double distance;
            size_t gap, previous, current;
        };
        std::vector<size_t> ordered(observations.size());
        for (size_t i = 0; i < ordered.size(); ++i)
            ordered[i] = i;
        std::sort(ordered.begin(), ordered.end(), [&](size_t a, size_t b)
        {
            return observations[a].mz < observations[b].mz;
        });
        std::vector<Match> matches;
        for (size_t i = 0; i < active.size(); ++i)
        {
            const auto &prior = active[i].last;
            auto first = std::lower_bound(ordered.begin(), ordered.end(), prior.lower,
                [&](size_t index, double lower) { return observations[index].mz < lower; });
            for (; first != ordered.end() && observations[*first].mz <= prior.upper; ++first)
            {
                const auto &current = observations[*first];
                if (prior.mz < current.lower || prior.mz > current.upper)
                    continue;
                const double distance = std::abs(prior.mz - current.mz);
                const double limit = std::max(2 * std::max(prior.step, current.step),
                    std::min(prior.upper - prior.lower, current.upper - current.lower) / 4);
                if (distance > limit)
                    continue;
                const bool shared = std::any_of(prior.heights.begin(), prior.heights.end(), [&](const auto &a)
                {
                    return std::any_of(current.heights.begin(), current.heights.end(),
                        [&](const auto &b) { return a.scan == b.scan; });
                });
                if (shared)
                    matches.push_back({distance, window - active[i].lastWindow, i, *first});
            }
        }
        std::sort(matches.begin(), matches.end(), [&](const Match &a, const Match &b)
        {
            return std::tie(a.distance, a.gap, active[a.previous].id, a.current) <
                   std::tie(b.distance, b.gap, active[b.previous].id, b.current);
        });
        std::vector<bool> used(active.size());
        std::vector<Track> next(observations.size());
        for (const auto &match : matches)
        {
            if (!used[match.previous] && next[match.current].id == 0)
            {
                used[match.previous] = true;
                next[match.current] = std::move(active[match.previous]);
            }
        }
        std::vector<size_t> ids;
        ids.reserve(observations.size());
        for (size_t i = 0; i < observations.size(); ++i)
        {
            auto &track = next[i];
            const auto &observation = observations[i];
            if (!track.id)
                track.id = nextId++;
            ++track.windows;
            track.lastWindow = window;
            const double delta = observation.mz - track.mean;
            track.mean += delta / track.windows;
            track.squares += delta * (observation.mz - track.mean);
            track.last = observation;
            for (const auto &height : observation.heights)
            {
                const auto inserted = track.heights.try_emplace(height.scan, SuperspectrumTrackHeight{height, firstScan});
                if (!inserted.second && height.distanceToWindowCenter <
                    inserted.first->second.measurement.distanceToWindowCenter)
                    inserted.first->second = {height, firstScan};
            }
            ids.push_back(track.id);
        }
        for (size_t i = 0; i < active.size(); ++i)
            if (!used[i])
            {
                if (window - active[i].lastWindow + 1 < windowScans)
                    next.push_back(std::move(active[i]));
                else
                    write(active[i]);
            }
        active = std::move(next);
        return ids;
    }

    void SuperspectrumTracks::finish()
    {
        for (const auto &track : active)
            write(track);
        active.clear();
    }
}
