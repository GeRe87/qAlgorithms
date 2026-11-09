#include "qalgorithms_multitrace.h"
#include "qalgorithms_qpeaks.h" // existing six-value (XtX)^-1 cache

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace qAlgorithms
{
    SharedShapeMatInverse sharedShapeInverseForScale(const size_t scale,
                                                      const size_t traceCount)
    {
        if (scale < GLOBAL_MINSCALE || scale > maxscale_global || traceCount == 0)
            throw std::invalid_argument("Invalid scale or number of shared-shape traces");
        const MatInverse &inv = qalgo_matInverse[scale];
        const double invM = 1.0 / static_cast<double>(traceCount);
        const double invN = 1.0 / static_cast<double>(2 * scale + 1);
        const double interceptCross = (inv.A - invN) * invM;
        return {
            invN + interceptCross, interceptCross,
            inv.B * invM, inv.C * invM, inv.D * invM,
            inv.E * invM, inv.F * invM
        };
    }

    void forEachSharedShapeCoefficients(
        const std::vector<const float *> &logTraces,
        const size_t length,
        const SharedShapeFitVisitor &visit)
    {
        if (logTraces.empty() || length < MINLENGTH || !visit)
            throw std::invalid_argument("Shared-shape regression requires traces, at least five samples and a visitor");
        if (std::any_of(logTraces.begin(), logTraces.end(),
                        [](const float *trace) { return trace == nullptr; }))
            throw std::invalid_argument("Shared-shape regression received a null trace");

        const size_t traceCount = logTraces.size();
        const double invTraceCount = 1.0 / static_cast<double>(traceCount);
        const size_t maxScale = std::min(maxscale_global, (length - 1) / 2);

        std::vector<double> sum0(traceCount, 0.0);
        SharedShapeRegCoeffs fit;
        fit.b0.resize(traceCount);

        // Reuse the rolling X^T*y sums for each centre, just as in
        // findCoefficients. Only the intercept sum is trace-specific.
        for (size_t centre = GLOBAL_MINSCALE; centre < length - GLOBAL_MINSCALE; ++centre)
        {
            double total0 = 0.0, total1 = 0.0, total2 = 0.0, total3 = 0.0;
            for (size_t j = 0; j < traceCount; ++j)
            {
                const float *data = logTraces[j] + centre;
                const double left = static_cast<double>(data[-1]);
                const double middle = static_cast<double>(data[0]);
                const double right = static_cast<double>(data[1]);
                sum0[j] = left + middle + right;
                total0 += sum0[j];
                total1 += right - left;
                total2 += left;
                total3 += right;
            }

            const size_t maxScaleAtCentre =
                std::min(maxScale, std::min(centre, length - centre - 1));
            for (size_t scale = GLOBAL_MINSCALE; scale <= maxScaleAtCentre; ++scale)
            {
                double leftTotal = 0.0, rightTotal = 0.0;
                for (size_t j = 0; j < traceCount; ++j)
                {
                    const double left = static_cast<double>(logTraces[j][centre - scale]);
                    const double right = static_cast<double>(logTraces[j][centre + scale]);
                    sum0[j] += left + right;
                    leftTotal += left;
                    rightTotal += right;
                }

                const double x = static_cast<double>(scale);
                const double x2 = x * x;
                total0 += leftTotal + rightTotal;
                total1 += x * (rightTotal - leftTotal);
                total2 += x2 * leftTotal;
                total3 += x2 * rightTotal;

                // The common shape equals the single-trace LS fit to the
                // per-point trace mean. Therefore the same A..F cache is exact.
                const double u0 = total0 * invTraceCount;
                const double u1 = total1 * invTraceCount;
                const double u2 = total2 * invTraceCount;
                const double u3 = total3 * invTraceCount;
                const MatInverse &inv = qalgo_matInverse[scale];
                const double meanB0 = inv.A * u0 + inv.B * (u2 + u3);

                fit.b1 = inv.C * u1 + inv.D * (u2 - u3);
                fit.b2 = inv.B * u0 + inv.D * u1 + inv.E * u2 + inv.F * u3;
                fit.b3 = inv.B * u0 - inv.D * u1 + inv.F * u2 + inv.E * u3;
                fit.x0 = centre;
                fit.scale = static_cast<uint16_t>(scale);

                const double invN = 1.0 / static_cast<double>(2 * scale + 1);
                for (size_t j = 0; j < traceCount; ++j)
                    fit.b0[j] = meanB0 + (sum0[j] - u0) * invN;

                visit(fit);
            }
        }
    }

    void findSharedShapeCoefficients(
        const std::vector<const float *> &logTraces,
        const size_t length,
        std::vector<SharedShapeRegCoeffs> *coeffs)
    {
        if (coeffs == nullptr)
            throw std::invalid_argument("Shared-shape output must not be null");
        if (logTraces.empty() || length < MINLENGTH)
            throw std::invalid_argument("Shared-shape regression requires traces and at least five samples");

        const size_t maxScale = std::min(maxscale_global, (length - 1) / 2);
        std::vector<size_t> offsets(maxScale + 1, 0);
        size_t count = 0;
        for (size_t scale = GLOBAL_MINSCALE; scale <= maxScale; ++scale)
        {
            offsets[scale] = count;
            count += length - 2 * scale;
        }
        coeffs->clear();
        coeffs->resize(count);
        forEachSharedShapeCoefficients(logTraces, length,
            [coeffs, &offsets](const SharedShapeRegCoeffs &fit)
            {
                const size_t offset = offsets[fit.scale] + fit.x0 - fit.scale;
                (*coeffs)[offset] = fit;
            });
    }
} // namespace qAlgorithms
