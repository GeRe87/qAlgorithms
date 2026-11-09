#include "qalgorithms_multitrace.h"
#include "qalgorithms_qpeaks.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace qAlgorithms;

namespace
{
    using Matrix = std::vector<std::vector<long double>>;

    bool near(const long double actual, const long double expected,
              const long double relative = 3e-6L, const long double absolute = 1e-15L)
    {
        return std::abs(actual - expected) <= absolute + relative * std::abs(expected);
    }

    // Independent Gauss-Jordan inverse of the actual (m+3)-column design.
    Matrix invert(Matrix matrix)
    {
        const size_t n = matrix.size();
        Matrix inverse(n, std::vector<long double>(n, 0.0L));
        for (size_t i = 0; i < n; ++i) inverse[i][i] = 1.0L;
        for (size_t i = 0; i < n; ++i)
        {
            size_t pivot = i;
            for (size_t j = i + 1; j < n; ++j)
                if (std::abs(matrix[j][i]) > std::abs(matrix[pivot][i])) pivot = j;
            if (std::abs(matrix[pivot][i]) < 1e-25L)
                throw std::runtime_error("Singular independent design");
            std::swap(matrix[i], matrix[pivot]);
            std::swap(inverse[i], inverse[pivot]);
            const long double divisor = matrix[i][i];
            for (size_t j = 0; j < n; ++j)
            {
                matrix[i][j] /= divisor;
                inverse[i][j] /= divisor;
            }
            for (size_t r = 0; r < n; ++r)
            {
                if (r == i) continue;
                const long double factor = matrix[r][i];
                for (size_t j = 0; j < n; ++j)
                {
                    matrix[r][j] -= factor * matrix[i][j];
                    inverse[r][j] -= factor * inverse[i][j];
                }
            }
        }
        return inverse;
    }

    Matrix designGram(const size_t scale, const size_t traces)
    {
        const size_t p = traces + 3;
        Matrix gram(p, std::vector<long double>(p, 0.0L));
        const int h = static_cast<int>(scale);
        for (size_t j = 0; j < traces; ++j)
            for (int x = -h; x <= h; ++x)
            {
                std::vector<long double> row(p, 0.0L);
                row[j] = 1.0L;
                row[traces] = static_cast<long double>(x);
                row[traces + 1] = x < 0 ? static_cast<long double>(x) * x : 0.0L;
                row[traces + 2] = x > 0 ? static_cast<long double>(x) * x : 0.0L;
                for (size_t a = 0; a < p; ++a)
                    for (size_t b = 0; b < p; ++b)
                        gram[a][b] += row[a] * row[b];
            }
        return gram;
    }

    bool checkInverseStructure(const size_t h, const size_t m)
    {
        const Matrix actual = invert(designGram(h, m));
        const SharedShapeMatInverse coeffs = sharedShapeInverseForScale(h, m);
        Matrix expected(m + 3, std::vector<long double>(m + 3, 0.0L));
        for (size_t j = 0; j < m; ++j)
        {
            for (size_t k = 0; k < m; ++k)
                expected[j][k] = j == k ? coeffs.interceptDiagonal
                                          : coeffs.interceptOffDiagonal;
            expected[j][m + 1] = expected[m + 1][j] = coeffs.interceptCurvature;
            expected[j][m + 2] = expected[m + 2][j] = coeffs.interceptCurvature;
        }
        expected[m][m] = coeffs.slopeVariance;
        expected[m][m + 1] = expected[m + 1][m] = coeffs.slopeCurvature;
        expected[m][m + 2] = expected[m + 2][m] = -coeffs.slopeCurvature;
        expected[m + 1][m + 1] = expected[m + 2][m + 2] = coeffs.curvatureVariance;
        expected[m + 1][m + 2] = expected[m + 2][m + 1] = coeffs.curvatureCovariance;
        for (size_t i = 0; i < m + 3; ++i)
            for (size_t j = 0; j < m + 3; ++j)
                if (!near(actual[i][j], expected[i][j]))
                {
                    std::cerr << "Inverse mismatch at h=" << h << " m=" << m
                              << " row=" << i << " col=" << j << " got="
                              << static_cast<double>(actual[i][j]) << " expected="
                              << static_cast<double>(expected[i][j]) << '\n';
                    return false;
                }
        return true;
    }

    bool checkRegression(const size_t h, const size_t m)
    {
        const size_t length = 2 * h + 1;
        std::vector<std::vector<float>> data(m, std::vector<float>(length));
        std::vector<const float *> pointers(m);
        for (size_t j = 0; j < m; ++j)
        {
            pointers[j] = data[j].data();
            for (size_t k = 0; k < length; ++k)
            {
                const double x = static_cast<double>(k) - static_cast<double>(h);
                const double curvature = x < 0 ? -0.002 : -0.003;
                const double perturb = 0.014 * std::sin(static_cast<double>(k + 13*j) * 0.7);
                data[j][k] = static_cast<float>(8.0 + static_cast<double>(j) * 0.3
                                     + 0.04 * x + curvature * x * x + perturb);
            }
        }
        // Independent normal-equation solution using the same float observations.
        const Matrix inverse = invert(designGram(h, m));
        std::vector<long double> rhs(m + 3, 0.0L);
        for (size_t j = 0; j < m; ++j)
            for (size_t k = 0; k < length; ++k)
            {
                const long double y = data[j][k];
                const long double x = static_cast<long double>(k) - static_cast<long double>(h);
                rhs[j] += y;
                rhs[m] += x * y;
                if (x < 0) rhs[m + 1] += x * x * y;
                else if (x > 0) rhs[m + 2] += x * x * y;
            }
        std::vector<long double> reference(m + 3, 0.0L);
        for (size_t i = 0; i < m + 3; ++i)
            for (size_t j = 0; j < m + 3; ++j)
                reference[i] += inverse[i][j] * rhs[j];

        SharedShapeRegCoeffs observed;
        size_t matches = 0;
        size_t visited = 0;
        forEachSharedShapeCoefficients(pointers, length,
            [&](const SharedShapeRegCoeffs &fit)
            {
                ++visited;
                if (fit.x0 == h && fit.scale == h)
                {
                    observed = fit;
                    ++matches;
                }
            });
        if (matches != 1 || observed.b0.size() != m || visited == 0)
            return false;
        for (size_t j = 0; j < m; ++j)
            if (!near(observed.b0[j], reference[j], 4e-5L, 3e-5L))
                return false;
        return near(observed.b1, reference[m], 4e-5L, 3e-5L) &&
               near(observed.b2, reference[m + 1], 4e-5L, 3e-5L) &&
               near(observed.b3, reference[m + 2], 4e-5L, 3e-5L);
    }

    bool checkCollectedOrder()
    {
        const size_t length = 9;
        std::vector<float> a{2, 3, 4, 7, 11, 8, 4, 2, 1};
        std::vector<float> b{6, 7, 8, 11, 15, 12, 8, 6, 5};
        std::vector<SharedShapeRegCoeffs> all;
        findSharedShapeCoefficients({a.data(), b.data()}, length, &all);
        size_t pos = 0;
        for (size_t scale = 2; scale <= 4; ++scale)
            for (size_t centre = scale; centre + scale < length; ++centre)
            {
                if (pos >= all.size() || all[pos].scale != scale || all[pos].x0 != centre ||
                    all[pos].b0.size() != 2 || !near(all[pos].b0[1] - all[pos].b0[0], 4))
                    return false;
                ++pos;
            }
        return pos == all.size();
    }

    bool invalidInputsThrow()
    {
        const std::vector<float> one(5, 1.0f);
        const SharedShapeFitVisitor visit = [](const SharedShapeRegCoeffs &) {};
        try { forEachSharedShapeCoefficients({}, 5, visit); return false; }
        catch (const std::invalid_argument &) {}
        try { forEachSharedShapeCoefficients({one.data()}, 4, visit); return false; }
        catch (const std::invalid_argument &) {}
        try { forEachSharedShapeCoefficients({nullptr}, 5, visit); return false; }
        catch (const std::invalid_argument &) {}
        try { forEachSharedShapeCoefficients({one.data()}, 5, {}); return false; }
        catch (const std::invalid_argument &) {}
        try { findSharedShapeCoefficients({one.data()}, 5, nullptr); return false; }
        catch (const std::invalid_argument &) {}
        try { sharedShapeInverseForScale(1, 1); return false; }
        catch (const std::invalid_argument &) {}
        try { sharedShapeInverseForScale(216, 1); return false; }
        catch (const std::invalid_argument &) {}
        try { sharedShapeInverseForScale(2, 0); return false; }
        catch (const std::invalid_argument &) {}
        return true;
    }
} // namespace

int main()
{
    for (const size_t h : {size_t(2),size_t(3),size_t(8),size_t(25),size_t(100),size_t(215)})
        for (const size_t m : {size_t(1),size_t(2),size_t(3),size_t(10)})
        {
            if (!checkInverseStructure(h, m) || !checkRegression(h, m))
            {
                std::cerr << "Multi-trace regression failed at scale=" << h
                          << " traces=" << m << '\n';
                return 1;
            }
        }
    if (!checkCollectedOrder() || !invalidInputsThrow())
    {
        std::cerr << "Multi-trace ordering or invalid-input tests failed\n";
        return 1;
    }
    std::cout << "Multi-trace Gram inverses and regressions: PASS\n";
    return 0;
}
