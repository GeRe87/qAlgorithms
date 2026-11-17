#ifndef QALGORITHMS_DETAIL_SINGLE_TRACE_HYPERBOLIC_H
#define QALGORITHMS_DETAIL_SINGLE_TRACE_HYPERBOLIC_H

#include "qalgorithms_qpeaks.h"

namespace qAlgorithms::detail
{
    struct SingleTraceHyperbolic
    {
        struct Coefficients
        {
            RegCoeffs values;
            double logIntercept;
        };

        struct Result : RegressionGauss
        {
            Coefficients coeffs;
        };

        static constexpr size_t parameterCount = 4;
        static constexpr bool omitUnsupportedApexGroups = true;
        static Span_i32 similaritySpan(Span_i32 span, size_t) { return span; }
        static size_t center(const Coefficients &c) { return c.values.x0; }
        static size_t halfWidth(const Coefficients &c) { return c.values.scale; }
        static bool supportsMeasuredCount(size_t measured, Span_i32 span) { return measured == size_t(span.length); }
        static void runningFit(const float *logs, size_t length, std::vector<Coefficients> *output);
        static Span_i32 regressionSpan(const Coefficients *c);
        static double logValue(const Coefficients *c, double offset);
        static double apex(const Coefficients *c);
        static bool apexHasSupport(float position, const Coefficients *c);
        static bool insignificantCurvature(const Coefficients *c, double mse, size_t df);
        static double refitAmplitude(const float *raw, Span_i32 span, float *prediction, Coefficients *c);
        static bool measure(Result &result, const float *raw, const float *logs, size_t length, double mse, double step);
        static double width(const Coefficients *c, double height, double step);
    };
}

#endif
