#include "qalgorithms_peak_fifo.h"

#include <cmath>
#include <stdexcept>

namespace
{
    constexpr double center = 100;

    void require(bool condition)
    {
        if (!condition)
            throw std::runtime_error("FIFO joint mass regression failed");
    }

    // Both nonuniform grids are symmetric, so the slope column is
    // orthogonal to the curvature and to all scan intercepts.
    void append(qAlgorithms::PeakSegmentFifo &fifo, size_t scan, double &xx, double &residual)
    {
        const double inner = scan % 2 ? .25 : .5;
        const std::vector<double> mz{99, 100 - inner, 100, 100 + inner, 101};
        std::vector<double> intensity;
        for (double mass : mz)
        {
            const double x = mass - center;
            const double noise = std::abs(x) > .9 ? .012 * x : -.012 * x / (inner * inner);
            intensity.push_back(std::exp(2 + .13 * scan - .8 * x * x + noise));
            xx += x * x;
            residual += noise * noise;
        }
        fifo.append(scan, mz, intensity);
    }

    void testOffCenter()
    {
        constexpr double apex = 100.27;
        qAlgorithms::PeakSegmentFifo fifo(99, 101, 5);
        double xx = 0, xq = 0, qq = 0, residual = 0;
        for (size_t scan = 0; scan < 6; ++scan)
        {
            const std::vector<double> mz{99.1, 99.55 + .03 * scan, 100.2, 100.95 - .01 * scan};
            double xMean = 0, qMean = 0;
            for (double mass : mz)
            {
                xMean += (mass - apex) / mz.size();
                qMean += (mass - apex) * (mass - apex) / mz.size();
            }
            std::vector<double> intensity;
            for (size_t i = 0; i < mz.size(); ++i)
            {
                const double x = mz[i] - apex;
                double noise = .005;
                // These weights annihilate all polynomials of degree <= 2.
                for (size_t j = 0; j < mz.size(); ++j)
                    if (i != j)
                        noise /= mz[i] - mz[j];
                intensity.push_back(std::exp(2 + scan * .2 - .8 * x * x + noise));
                residual += noise * noise;
                xx += (x - xMean) * (x - xMean);
                xq += (x - xMean) * (x * x - qMean);
                qq += (x * x - qMean) * (x * x - qMean);
            }
            fifo.append(scan, mz, intensity);
        }
        const auto fit = fifo.fitAll();
        const double expected = std::sqrt((residual / 16) / (xx - xq * xq / qq)) / 1.6;
        require(fit && std::abs(fit->mz - apex) < 1e-10);
        require(std::abs(fit->standardError - expected) < 1e-10);
    }
}

int main()
{
    qAlgorithms::PeakSegmentFifo fifo(99, 101, 5);
    double xx = 0, residual = 0;
    const std::vector<size_t> scans{0, 1, 3, 4, 6, 7, 10};
    for (size_t scan : scans)
    {
        append(fifo, scan, xx, residual);
        qAlgorithms::PeakSegmentFifo batch(99, 101, 5);
        double batchXx = 0, batchResidual = 0;
        for (size_t index : scans)
            if (index <= scan && scan - index < 5)
                append(batch, index, batchXx, batchResidual);
        const auto running = fifo.fitWindow();
        const auto direct = batch.fitAll();
        require(running.has_value() == direct.has_value());
        if (running)
        {
            require(running->scans == direct->scans && running->points == direct->points);
            require(std::abs(running->mz - direct->mz) < 1e-10);
            require(std::abs(running->standardError - direct->standardError) < 1e-10);
        }
    }

    const auto all = fifo.fitAll();
    require(all && all->scans == scans.size() && all->points == scans.size() * 5);
    require(std::abs(all->mz - center) < 1e-7);
    const double expected = std::sqrt(residual / (all->points - all->scans - 2) / xx) / 1.6;
    require(std::abs(all->standardError - expected) < 1e-7);
    fifo.advance(15);
    require(!fifo.fitWindow() && fifo.fitAll()->points == all->points);

    bool duplicateRejected = false;
    try
    {
        append(fifo, 10, xx, residual);
    }
    catch (const std::invalid_argument &)
    {
        duplicateRejected = true;
    }
    require(duplicateRejected && fifo.fitAll()->points == all->points);

    qAlgorithms::PeakSegmentFifo unsupported(99, 101, 5);
    unsupported.append(0, {99, 99.5f, 100}, {1, 2, 3});
    unsupported.append(1, {100, 100.5f, 101}, {3, 2, 1});
    require(!unsupported.fitAll());
    testOffCenter();
}
