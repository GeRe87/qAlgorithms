#ifndef QALGORITHMS_MULTITRACE_H
#define QALGORITHMS_MULTITRACE_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace qAlgorithms
{
    // The intercepts belong to separate, equally sampled log-intensity traces.
    // The slope and both piecewise quadratic coefficients are shared:
    // y_j(x) = b0[j] + b1*x + (x < 0 ? b2 : b3)*x*x.
    struct SharedShapeRegCoeffs
    {
        std::vector<double> b0;
        double b1 = 0, b2 = 0, b3 = 0;
        size_t x0 = 0;
        uint16_t scale = 0;
    };

    // The seven unique nonzero entries of the inverse normal matrix for
    // any m >= 2. Index order: b0[0..m-1], b1, b2, b3. For m=1, the
    // diagonal/off-diagonal intercept distinction degenerates harmlessly.
    struct SharedShapeMatInverse
    {
        double interceptDiagonal = 0;    // each (b0_j, b0_j)
        double interceptOffDiagonal = 0; // each (b0_j, b0_k), j != k
        double interceptCurvature = 0;   // each (b0_j, b2) and (b0_j, b3)
        double slopeVariance = 0;        // (b1, b1)
        double slopeCurvature = 0;       // (b1, b2); (b1, b3) is negative
        double curvatureVariance = 0;    // (b2, b2) and (b3, b3)
        double curvatureCovariance = 0;  // (b2, b3)
    };

    // Derive all seven multi-trace inverse values from the EXISTING six-value
    // single-trace cache. Neither a new cache nor a larger matrix is needed.
    // Throws std::invalid_argument for scale outside [2,215] or traceCount=0.
    SharedShapeMatInverse sharedShapeInverseForScale(size_t scale, size_t traceCount);

    using SharedShapeFitVisitor = std::function<void(const SharedShapeRegCoeffs &)>;

    // Streaming, O(number of traces) scratch-space running regression.
    // Each trace must point to 'length' log-intensity samples at the SAME
    // equidistant x positions. Every window uses -scale ... +scale with scale
    // in [2, min(215, (length-1)/2)], including its centre point.
    //
    // Fits are visited in centre-major, then increasing-scale order. The b0
    // storage in the callback view is reused: copy it if it must be retained.
    // Throws std::invalid_argument on empty input, nullptrs or length < 5.
    void forEachSharedShapeCoefficients(
        const std::vector<const float *> &logTraces,
        size_t length,
        const SharedShapeFitVisitor &visit);

    // Collecting equivalent of the streaming API. Output order is scale-major,
    // matching findCoefficients; memory is proportional to the result size.
    void findSharedShapeCoefficients(
        const std::vector<const float *> &logTraces,
        size_t length,
        std::vector<SharedShapeRegCoeffs> *coeffs);
} // namespace qAlgorithms

#endif // QALGORITHMS_MULTITRACE_H
