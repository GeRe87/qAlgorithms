#include "qalgorithms_peak_fifo.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Accumulate a common local log-quadratic shape with one intercept per scan.
namespace qAlgorithms
{
    PeakSegmentFifo::PeakSegmentFifo(double lower, double upper, size_t windowScans)
        : lower(lower), upper(upper), origin(lower + (upper - lower) / 2),
          scale((upper - lower) / 2), windowScans(windowScans)
    {
        if (!(lower > 0 && scale > 0 && std::isfinite(upper) && windowScans >= 2))
            throw std::invalid_argument("Invalid FIFO peak segment or window size");
    }

    void PeakSegmentFifo::update(Totals &totals, const Scan &scan, bool add)
    {
        for (size_t i = 0; i < totals.sums.size(); ++i)
            totals.sums[i] += add ? scan.sums[i] : -scan.sums[i];
        if (add)
        {
            ++totals.scans;
            totals.points += scan.mz.size();
        }
        else
        {
            --totals.scans;
            totals.points -= scan.mz.size();
            if (totals.scans == 0)
                totals = {};
        }
    }

    void PeakSegmentFifo::advance(size_t scan)
    {
        if (lastScan && scan < *lastScan)
            throw std::invalid_argument("FIFO scans must be advanced in acquisition order");
        lastScan = scan;
        while (firstActive < scans.size() && scan - scans[firstActive].index >= windowScans)
        {
            update(window, scans[firstActive], false);
            update(expired, scans[firstActive], true);
            ++firstActive;
        }
    }

    void PeakSegmentFifo::append(size_t index, std::vector<double> mz, std::vector<double> intensity)
    {
        if ((lastScan && index < *lastScan) || (!scans.empty() && index <= scans.back().index))
            throw std::invalid_argument("Each FIFO scan may contribute at most once in acquisition order");
        if (mz.size() != intensity.size() || mz.size() < 3)
            throw std::invalid_argument("A FIFO trace needs at least three paired measurements");

        Scan scan{index, std::move(mz), std::move(intensity)};
        std::array<long double, 3> mean{};
        for (size_t i = 0; i < scan.mz.size(); ++i)
        {
            if (!(scan.mz[i] >= lower && scan.mz[i] <= upper &&
                  (!i || scan.mz[i] > scan.mz[i - 1]) &&
                  scan.intensity[i] > 0 && std::isfinite(scan.intensity[i])))
                throw std::invalid_argument("FIFO trace must have ordered in-segment masses and finite positive intensities");
            const long double x = (static_cast<long double>(scan.mz[i]) - origin) / scale;
            const std::array<long double, 3> point{x, x * x, std::log(static_cast<long double>(scan.intensity[i]))};
            const std::array<long double, 3> delta{point[0] - mean[0], point[1] - mean[1], point[2] - mean[2]};
            const long double weight = static_cast<long double>(i) / (i + 1);
            // Online centering eliminates the scan-specific log amplitude.
            scan.sums[0] += weight * delta[0] * delta[0];
            scan.sums[1] += weight * delta[0] * delta[1];
            scan.sums[2] += weight * delta[1] * delta[1];
            scan.sums[3] += weight * delta[0] * delta[2];
            scan.sums[4] += weight * delta[1] * delta[2];
            scan.sums[5] += weight * delta[2] * delta[2];
            for (size_t j = 0; j < mean.size(); ++j)
                mean[j] += delta[j] / (i + 1);
        }
        // Keep original samples and scan ownership until the segment closes.
        scans.push_back(std::move(scan));
        update(window, scans.back(), true);
        advance(index);
    }

    std::optional<SharedMassFit> PeakSegmentFifo::fit(const Totals &totals, size_t firstScan) const
    {
        if (totals.scans < 2 || totals.points <= totals.scans + 2)
            return std::nullopt;
        const auto &[xx, xq, qq, xy, qy, yy] = totals.sums;
        const long double determinant = xx * qq - xq * xq;
        const long double epsilon = 64 * std::numeric_limits<long double>::epsilon();
        if (!(xx > 0 && qq > 0 && determinant > epsilon * xx * qq))
            return std::nullopt;
        const long double slope = (qq * xy - xq * qy) / determinant;
        const long double curvature = (xx * qy - xq * xy) / determinant;
        if (!(curvature < 0))
            return std::nullopt;
        const long double apex = -slope / (2 * curvature);
        const double mz = static_cast<double>(origin + scale * apex);
        for (size_t i = firstScan; i < scans.size(); ++i)
        {
            if (!(mz > scans[i].mz.front() && mz < scans[i].mz.back()))
                return std::nullopt;
        }

        const long double explained = slope * xy + curvature * qy;
        const long double residual = yy - explained;
        if (residual < -epsilon * (yy + std::abs(explained)))
            return std::nullopt;
        const long double variance = std::max(0.0L, residual) / (totals.points - totals.scans - 2);
        // Delta method for the apex -slope / (2 * curvature), converted to m/z.
        // Evaluate the inverse-Gram quadratic form as a sum of squares.
        const long double gradient = -scale / (2 * curvature);
        const long double apexVariance = variance * gradient * gradient / xx *
            (1 + xx * xx / determinant * (2 * apex - xq / xx) * (2 * apex - xq / xx));
        const double error = static_cast<double>(std::sqrt(apexVariance));
        if (!std::isfinite(error))
            return std::nullopt;
        return SharedMassFit{mz, error, totals.scans, totals.points};
    }

    std::optional<SharedMassFit> PeakSegmentFifo::fitWindow() const
    {
        return fit(window, firstActive);
    }

    std::optional<SharedMassFit> PeakSegmentFifo::fitAll() const
    {
        auto all = expired;
        for (size_t i = 0; i < all.sums.size(); ++i)
            all.sums[i] += window.sums[i];
        all.scans += window.scans;
        all.points += window.points;
        return fit(all, 0);
    }
}
