#include "qalgorithms_qpeaks.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace qAlgorithms;

namespace
{
    using LD = long double;
    using Matrix = std::vector<std::vector<LD>>;

    bool near(const LD actual, const LD expected,
              const LD relative = 2e-7L, const LD absolute = 3e-9L)
    {
        return std::abs(actual - expected) <= absolute + relative * std::abs(expected);
    }

    Matrix inverse(Matrix a)
    {
        const size_t p = a.size();
        Matrix out(p, std::vector<LD>(p, 0.0L));
        for (size_t j = 0; j < p; ++j) out[j][j] = 1.0L;
        for (size_t col = 0; col < p; ++col)
        {
            size_t pivot = col;
            for (size_t r = col + 1; r < p; ++r)
                if (std::abs(a[r][col]) > std::abs(a[pivot][col])) pivot = r;
            if (!(std::abs(a[pivot][col]) > 1e-22L))
                throw std::runtime_error("Independent inverse is singular");
            std::swap(a[col], a[pivot]);
            std::swap(out[col], out[pivot]);
            const LD divisor = a[col][col];
            for (size_t j = 0; j < p; ++j)
            {
                a[col][j] /= divisor;
                out[col][j] /= divisor;
            }
            for (size_t r = 0; r < p; ++r)
            {
                if (r == col) continue;
                const LD factor = a[r][col];
                for (size_t j = 0; j < p; ++j)
                {
                    a[r][j] -= factor * a[col][j];
                    out[r][j] -= factor * out[col][j];
                }
            }
        }
        return out;
    }

    struct Reference
    {
        std::vector<LD> coefficients;
        Matrix covarianceFactor;
        LD rss = 0.0L;
        size_t df = 0;
    };

    Reference fullOLS(const std::vector<FoldingProfileView> &profiles,
                      const FoldingCurvatureMode mode)
    {
        const size_t m = profiles.size();
        const size_t p = mode == FoldingCurvatureMode::symmetric ? 1 : 2;
        const size_t columns = 2 * m + p;
        size_t samples = 0;
        Matrix gram(columns, std::vector<LD>(columns, 0.0L));
        std::vector<LD> rhs(columns, 0.0L);

        for (size_t j = 0; j < m; ++j)
            for (size_t i = 0; i < profiles[j].length; ++i)
            {
                const LD x = profiles[j].u[i];
                const LD y = profiles[j].logIntensity[i];
                std::vector<LD> row(columns, 0.0L);
                row[2 * j] = 1.0L;
                row[2 * j + 1] = x;
                if (p == 1) row[2 * m] = x * x;
                else if (x < 0) row[2 * m] = x * x;
                else row[2 * m + 1] = x * x;
                for (size_t k = 0; k < columns; ++k)
                {
                    rhs[k] += row[k] * y;
                    for (size_t l = 0; l < columns; ++l)
                        gram[k][l] += row[k] * row[l];
                }
                ++samples;
            }

        const Matrix cov = inverse(gram);
        std::vector<LD> coeff(columns, 0.0L);
        for (size_t k = 0; k < columns; ++k)
            for (size_t l = 0; l < columns; ++l)
                coeff[k] += cov[k][l] * rhs[l];

        LD rss = 0.0L;
        for (size_t j = 0; j < m; ++j)
            for (size_t i = 0; i < profiles[j].length; ++i)
            {
                const LD x = profiles[j].u[i];
                const LD y = profiles[j].logIntensity[i];
                const LD curvature = p == 1 || x < 0
                    ? coeff[2 * m] : coeff[2 * m + 1];
                const LD predicted = coeff[2 * j] + coeff[2 * j + 1] * x
                                   + curvature * x * x;
                const LD residual = y - predicted;
                rss += residual * residual;
            }
        return {coeff, cov, rss, samples - columns};
    }

    struct DataSet
    {
        std::vector<std::vector<double>> coordinates;
        std::vector<std::vector<float>> intensities;

        std::vector<FoldingProfileView> views() const
        {
            std::vector<FoldingProfileView> result;
            for (size_t j = 0; j < coordinates.size(); ++j)
                result.push_back({coordinates[j].data(), intensities[j].data(),
                                  coordinates[j].size()});
            return result;
        }
    };

    DataSet makeData(const size_t traceCount, const FoldingCurvatureMode mode)
    {
        DataSet data;
        for (size_t j = 0; j < traceCount; ++j)
        {
            const size_t n = 7 + j % 5;
            std::vector<double> x(n);
            std::vector<float> y(n);
            for (size_t i = 0; i < n; ++i)
            {
                const double xi = -3.25 + 6.65 * static_cast<double>(i) /
                                    static_cast<double>(n - 1)
                                  + 0.06 * static_cast<double>(j % 4)
                                  + 0.015 * std::sin(static_cast<double>(i + 3 * j));
                const double curvature = mode == FoldingCurvatureMode::symmetric
                    ? -0.175 : (xi < 0 ? -0.24 : -0.13);
                const double offset = 6.0 + 0.12 * static_cast<double>(j);
                const double slope = 0.07 * static_cast<double>(static_cast<int>(j % 4) - 2);
                const double noise = 0.009 * std::sin(static_cast<double>(2 * i + 5 * j));
                x[i] = xi;
                y[i] = static_cast<float>(offset + slope * xi
                                       + curvature * xi * xi + noise);
            }
            data.coordinates.push_back(std::move(x));
            data.intensities.push_back(std::move(y));
        }
        return data;
    }

    bool checkAgainstIndependentInverse(const size_t m,
                                        const FoldingCurvatureMode mode)
    {
        const DataSet data = makeData(m, mode);
        const auto views = data.views();
        const FoldingLinearFit fit = fitPeakFoldingLinear(views, mode);
        const Reference ref = fullOLS(views, mode);
        const size_t p = mode == FoldingCurvatureMode::symmetric ? 1 : 2;
        const size_t global = 2 * m;

        if (fit.mode != mode || fit.peaks.size() != m ||
            fit.degreesOfFreedom != ref.df ||
            !near(fit.residualSumSquares, ref.rss) ||
            !near(fit.mse, ref.rss / static_cast<LD>(ref.df)))
            return false;
        if (!near(fit.b2, ref.coefficients[global]) ||
            !near(fit.b3, ref.coefficients[global + p - 1]))
            return false;

        const size_t c0 = global, c1 = global + p - 1;
        if (!near(fit.curvatureCovarianceFactor[0], ref.covarianceFactor[c0][c0]) ||
            !near(fit.curvatureCovarianceFactor[1], ref.covarianceFactor[c0][c1]) ||
            !near(fit.curvatureCovarianceFactor[2], ref.covarianceFactor[c1][c0]) ||
            !near(fit.curvatureCovarianceFactor[3], ref.covarianceFactor[c1][c1]))
            return false;

        for (size_t j = 0; j < m; ++j)
        {
            const size_t a = 2 * j, b = 2 * j + 1;
            const auto &o = fit.peaks[j];
            if (!near(o.b0, ref.coefficients[a]) ||
                !near(o.b1, ref.coefficients[b]) ||
                !near(o.varianceB0Factor, ref.covarianceFactor[a][a]) ||
                !near(o.varianceB1Factor, ref.covarianceFactor[b][b]) ||
                !near(o.covarianceB0B1Factor, ref.covarianceFactor[a][b]) ||
                !near(o.covarianceB1B2Factor, ref.covarianceFactor[b][c0]) ||
                !near(o.covarianceB1B3Factor, ref.covarianceFactor[b][c1]))
                return false;

            if (!o.hasInteriorApex)
            {
                std::cerr << "Synthetic peak should have an interior apex\n";
                return false;
            }
            const bool left = o.b1 < 0;
            const LD curvature = left ? ref.coefficients[c0] : ref.coefficients[c1];
            const LD xApex = -ref.coefficients[b] / (2 * curvature);
            const LD db = -1 / (2 * curvature);
            const LD dc = ref.coefficients[b] / (2 * curvature * curvature);
            const size_t c = left ? c0 : c1;
            const LD apexVar = db * db * ref.covarianceFactor[b][b]
                + dc * dc * ref.covarianceFactor[c][c]
                + 2 * db * dc * ref.covarianceFactor[b][c];
            if (!near(o.apexU, xApex) ||
                !near(o.apexVarianceFactor, apexVar) ||
                !near(o.apexStandardError, std::sqrt(fit.mse * static_cast<double>(apexVar))))
                return false;
        }
        return true;
    }

    bool testPhysicalPositionAndPooling()
    {
        // A translated symmetric quadratic is exactly represented by each
        // trace's own b0 and b1. No nonlinear alignment is required.
        DataSet data;
        for (size_t j = 0; j < 6; ++j)
        {
            std::vector<double> x;
            std::vector<float> y;
            const double delta = 0.5 * static_cast<double>(static_cast<int>(j) - 3) / 3.0;
            for (int i = -4; i <= 4; ++i)
            {
                const double xi = static_cast<double>(i) + 0.07 * static_cast<double>(j % 3);
                x.push_back(xi);
                y.push_back(static_cast<float>(5.0 + 0.3 * static_cast<double>(j)
                    - 0.18 * (xi - delta) * (xi - delta)));
            }
            data.coordinates.push_back(std::move(x));
            data.intensities.push_back(std::move(y));
        }
        const auto views = data.views();
        const auto result = fitPeakFoldingLinear(views, FoldingCurvatureMode::symmetric);
        for (size_t j = 0; j < views.size(); ++j)
        {
            const double delta = 0.5 * static_cast<double>(static_cast<int>(j) - 3) / 3.0;
            if (!result.peaks[j].hasInteriorApex ||
                !near(result.peaks[j].apexU, delta, 2e-5L, 2e-5L))
                return false;
        }

        // Pooling common curvature must not increase the model-based slope
        // variance factor relative to a single trace (correct shared model).
        for (size_t j = 0; j < views.size(); ++j)
        {
            const auto single = fitPeakFoldingLinear({views[j]}, FoldingCurvatureMode::symmetric);
            if (result.peaks[j].varianceB1Factor >
                single.peaks[0].varianceB1Factor + 1e-10)
                return false;
        }
        return true;
    }

    bool testNonPeakCurvature()
    {
        // One side is convex. Even when the other side has a stationary
        // point, it is NOT an isolated mass peak with a reliable apex.
        const double u[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
        float y[7] = {};
        for (size_t i = 0; i < 7; ++i)
        {
            const double x = u[i];
            const double curvature = x < 0 ? -0.18 : 0.14;
            y[i] = static_cast<float>(5.0 - 0.14 * x + curvature * x * x);
        }
        const auto fit = fitPeakFoldingLinear({{u, y, 7}},
                                             FoldingCurvatureMode::asymmetric);
        return fit.b2 < 0.0 && fit.b3 > 0.0 && !fit.peaks[0].hasInteriorApex;
    }

    bool testInvalidAndSingular()
    {
        const auto data = makeData(1, FoldingCurvatureMode::asymmetric);
        auto v = data.views();
        const auto expectInvalid = [](const auto &call)
        {
            try { call(); }
            catch (const std::invalid_argument &) { return true; }
            return false;
        };
        const auto expectSingular = [](const auto &call)
        {
            try { call(); }
            catch (const std::domain_error &) { return true; }
            return false;
        };

        if (!expectInvalid([&](){ fitPeakFoldingLinear({}, FoldingCurvatureMode::symmetric); }))
            return false;
        if (!expectInvalid([&](){
            fitPeakFoldingLinear(v, static_cast<FoldingCurvatureMode>(100)); }))
            return false;
        auto invalid = v;
        invalid[0].u = nullptr;
        if (!expectInvalid([&](){ fitPeakFoldingLinear(invalid); })) return false;
        invalid = v;
        invalid[0].logIntensity = nullptr;
        if (!expectInvalid([&](){ fitPeakFoldingLinear(invalid); })) return false;
        invalid = v;
        invalid[0].length = 2;
        if (!expectInvalid([&](){ fitPeakFoldingLinear(invalid); })) return false;

        std::vector<double> constant(6, 1.0);
        std::vector<float> signal(6, 3.0f);
        if (!expectSingular([&](){
            fitPeakFoldingLinear({{constant.data(), signal.data(), 6}}); }))
            return false;
        const double onlyPositive[6] = {0.3, 0.5, 0.7, 0.9, 1.1, 1.3};
        if (!expectSingular([&](){
            fitPeakFoldingLinear({{onlyPositive, signal.data(), 6}},
                                 FoldingCurvatureMode::asymmetric); }))
            return false;
        const double x3[3] = {-1.0, 0.0, 1.0};
        const float y3[3] = {2.0f, 3.0f, 2.0f};
        if (!expectInvalid([&](){
            fitPeakFoldingLinear({{x3, y3, 3}}, FoldingCurvatureMode::symmetric); }))
            return false;
        std::vector<double> nanCoordinates = data.coordinates[0];
        nanCoordinates[2] = std::numeric_limits<double>::quiet_NaN();
        if (!expectInvalid([&](){
            fitPeakFoldingLinear({{nanCoordinates.data(), data.intensities[0].data(),
                                   nanCoordinates.size()}}); }))
            return false;
        std::vector<float> infSignal = data.intensities[0];
        infSignal[1] = std::numeric_limits<float>::infinity();
        if (!expectInvalid([&](){
            fitPeakFoldingLinear({{data.coordinates[0].data(), infSignal.data(),
                                   infSignal.size()}}); }))
            return false;
        return true;
    }
}

int main()
{
    for (const auto mode : {FoldingCurvatureMode::symmetric,
                            FoldingCurvatureMode::asymmetric})
        for (const size_t m : {size_t(1), size_t(2), size_t(4), size_t(10)})
        {
            if (!checkAgainstIndependentInverse(m, mode))
            {
                std::cerr << "Folding failed full-OLS comparison, m=" << m
                          << " mode=" << (mode == FoldingCurvatureMode::asymmetric
                                                ? "asymmetric" : "symmetric") << '\n';
                return 1;
            }
        }

    if (!testPhysicalPositionAndPooling() || !testNonPeakCurvature() ||
        !testInvalidAndSingular())
    {
        std::cerr << "Folding shifted-apex, pooling or edge-case test failed\n";
        return 1;
    }
    std::cout << "Linear peak folding, covariance and apex uncertainty: PASS\n";
    return 0;
}
