// qalgorithms_qpeaks.h
#ifndef QALGORITHMS_QPEAKS_H
#define QALGORITHMS_QPEAKS_H

// internal
#include "qalgorithms_datatypes.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <vector>

namespace qAlgorithms
{
    // @todo rework this module so all centroiding / feature construction parts are in a separate block and qPeaks can be used as a library header

    // Optional collective description: resultIndex indexes the appended result
    // vector; start/end and plateau-first maximaIndices index the original input.
    // Bounds are inclusive and cover the observed/fitted-apex hypothesis envelope.
    // A closed expansion may span separately retained protected peaks; only
    // the maxima owned by this unresolved description appear in maximaIndices.
    // Representative-fit areas and overlapping individual areas are not additive.
    // Infeasible means no joint assignment exists; SearchLimit means unproven
    // within the explicit budget. Neither establishes a physical peak count.
    struct UnresolvedPeak
    {
        static constexpr uint8_t Infeasible = 1, SearchLimit = 2;
        size_t resultIndex = 0;
        uint8_t reason = Infeasible;
        uint32_t start = 0, end = 0;
        std::vector<uint32_t> maximaIndices;
    };

    // ### General Implementation ### //

    /// @brief This function serves as the central interface for performing peak detection on an arbitrary set of x and y data.
    /// @param intensity_base Values used to determine height of the peak
    /// @param x_values Values used as width of the peak, same length as y_values. It is a requirement that x values are evenly
    ///                 spaced when applying this function, since it depends on that assumption for its precalculated matric transform.
    ///                 We use "x" because this is either m/z or RT.
    /// @param DF_cum cumulative degrees of freedom to account for interpolations and mean values. This is the cumulative
    //                              sum of the number of measured points for all y. It has the same length as x and y.
    //                              This can be set to nullptr if all df are one.
    /// @param length length of y, x and df if supplied
    /// @param maxscale up to which scale regressions should be attempted. The maximum peak width (in points) is maxscale * 2 + 1.
    //                  If maxscale would exceed the length of x or y, it is set to (the length of either -1) / 2.
    /// @param result All found peaks are appended to this vector.
    /// @param unresolvedGroups nullptr preserves strict selection (also used by centroiding).
    ///        A supplied vector enables explicitly marked collective descriptions
    ///        when a single existing accepted fit covers the entire unresolved
    ///        envelope without replacing protected fits. Metadata is appended;
    ///        fit quantities belong to that representative, not a sum of analytes.
    /// @return the total number of appended peaks. If this is < 0, the function did not execute correctly.
    /// @details Errorcodes:
    ///  0 = no valid peaks were found (this is not an error per se, when supplying ex. a constant y this is the expected result)
    /// -1 = one of y_values, x_values or detectedPeaks was nullptr
    /// -2 = lenght is smaller than minimum required
    /// ### Assumptions: ###  these are not tested @todo
    /// x and degreesOfFreedom_cum increase monotonically.
    /// y has equal variance at every point
    /// there is enough space to write the results
    int32_t qpeaks_find(
        const float *intensity_base,
        const float *x_values,
        const uint16_t *DF_cum,
        const size_t length,
        std::vector<RegressionGauss> *result,
        std::vector<UnresolvedPeak> *unresolvedGroups = nullptr);

    enum class PeakModelKind { PiecewiseQuadratic, Hyperbolic };
    struct PeakFitOptions
    {
        PeakModelKind model = PeakModelKind::PiecewiseQuadratic;
    };

    // Selection is checked before input processing or result mutation.
    // Hyperbolic (and unknown enum values) throw std::invalid_argument in Stage A.
    int32_t qpeaks_find(const float *intensity, const float *x,
        const uint16_t *df, size_t length, std::vector<RegressionGauss> *result,
        const PeakFitOptions &options, std::vector<UnresolvedPeak> *groups = nullptr);

    namespace peak_model
    {
        // Units are supplied explicitly: the adapter cannot infer m/z versus RT.
        enum class AxisUnit { SampleIndex, MassToCharge, RetentionTime, CallerDefined };
        struct PeakProperties
        {
            double position, height, width, area;
            double positionUncertainty, heightUncertainty, areaUncertainty;
            bool valid;
            AxisUnit axisUnit; // width/position: axis; area: intensity*axis
            // Height and its uncertainty retain the caller's intensity unit.
        };
        PeakProperties fromLegacy(const RegressionGauss &, AxisUnit);

        // Intentionally incomplete primary trait: unsupported models have no
        // mathematical or statistical capabilities by accidental fallback.
        template<PeakModelKind> struct Model;
        template<> struct Model<PeakModelKind::PiecewiseQuadratic>
        {
            using Coefficients = RegCoeffs;
            // Historical form space, including a clipped non-apex flank.
            static bool hasLegacyPeakShape(const Coefficients &c)
            {
                const bool apexLeft = c.b1 < 0;
                return !((c.b2 >= 0 && apexLeft) || (c.b3 >= 0 && !apexLeft));
            }
            static void runningFit(const float *, size_t, std::vector<Coefficients> *);
            static double logValue(const Coefficients *, double);
            static double refitAmplitude(const float *, Span_i32, float *, Coefficients *);
            static double apex(const Coefficients *);
            static double width(const Coefficients *, double height, double step);
            static double integral(const Coefficients *, double step, double mse, double *uncertainty);
            static double positionUncertainty(const Coefficients *, double mse);
            static double heightUncertainty(const Coefficients *, double mse);
            static double covarianceProjection(const double gradient[4], size_t scale);
            // Legacy interior-apex gate; integrability is handled by the legacy
            // closed integral, which also supports its historical mixed tails.
            static bool apexHasSupport(float position, const Coefficients *);
        };
    }

    struct XML_File;                           // forward declaration so at least the header does not couple with read_file
    size_t findCentroids(const XML_File *data, // @todo the internal file representation should be more generic
                         const std::vector<uint32_t> *selectedIndices,
                         std::vector<CentroidPeak> *centroids);

    void centroids_to_mzml(const std::filesystem::path *pathSource);

    /// @param intensity_log logarithmy of the intensity values. This is the y axis of the fit. The x axis is required to be equidistant.
    /// @param maxscale maximum scale of a peak that should be attempted to fit.
    /// @param coeffs Sets of coefficients for all possible regressions. They are written out in the order scale = 2, scale = 3, ... , scale = maxscale
    //           Note that not all scales have the same number of regressions
    void findCoefficients(
        const float *intensity_log,
        const size_t length,
        std::vector<RegCoeffs> *coeffs);

    // Multi-trace running regression (same peak shape, trace-specific intercepts).
    // Kept in qPeaks alongside the original single-trace findCoefficients.
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
    void forEachCoefficientsMultiTrace(
        const std::vector<const float *> &logTraces,
        size_t length,
        const SharedShapeFitVisitor &visit);

    // Collecting equivalent of the streaming API. Output order is scale-major,
    // matching findCoefficients; memory is proportional to the result size.
    void findCoefficientsMultiTrace(
        const std::vector<const float *> &logTraces,
        size_t length,
        std::vector<SharedShapeRegCoeffs> *coeffs);

    // Experimental linear folding of several independent mass profiles.
    //
    // log(I_j) = b0[j] + b1[j]*u +
    //            (u < 0 ? b2 : b3)*u*u
    //
    // All profiles may have DIFFERENT u coordinates and point counts. They
    // must use consistent local units, with u=0 near the candidate peak.
    // The curvature breakpoint remains at u=0; an individual slope changes
    // the fitted apex but does NOT translate that breakpoint. Symmetric mode
    // represents translated quadratic peaks exactly.
    enum class FoldingCurvatureMode : uint8_t
    {
        symmetric,
        asymmetric
    };

    struct FoldingProfileView
    {
        const double *u = nullptr;       // local instrument sampling coordinates
        const float *logIntensity = nullptr; // already logarithmized, strictly finite
        size_t length = 0;
    };

    struct FoldingPeakEstimate
    {
        double b0 = 0, b1 = 0;
        double residualSumSquares = 0;
        // Entries of (X^T X)^-1 before multiplication by the pooled log-MSE:
        double varianceB0Factor = 0, varianceB1Factor = 0;
        double covarianceB0B1Factor = 0;
        double covarianceB1B2Factor = 0, covarianceB1B3Factor = 0;
        // Conditional delta-method apex uncertainty. Assumes a concave side
        // and an apex within the observed u range. Does not include anchor,
        // coordinate-calibration, selection or model-shape uncertainty.
        bool hasInteriorApex = false;
        double apexU = 0, apexVarianceFactor = 0, apexStandardError = 0;
    };

    struct FoldingLinearFit
    {
        FoldingCurvatureMode mode = FoldingCurvatureMode::asymmetric;
        double b2 = 0, b3 = 0;
        // Row-major covariance of (b2,b3), divided by the residual MSE.
        // In symmetric mode both coefficients are identical and perfectly
        // correlated, so the expanded 2x2 covariance is rank one.
        std::array<double, 4> curvatureCovarianceFactor = {};
        double residualSumSquares = 0, mse = 0;
        size_t degreesOfFreedom = 0;
        std::vector<FoldingPeakEstimate> peaks;
    };

    // Ordinary least squares with trace-specific intercepts and slopes and
    // shared quadratic curvature(s). Eliminates 2x2 local blocks and solves
    // only a 1x1 (symmetric) or 2x2 (asymmetric) Schur complement.
    // Scratch space is O(number of profiles); no interpolation, nonlinear
    // optimization or growing inverse matrix. Throws std::invalid_argument
    // for invalid/nonfinite data or too few residual degrees of freedom,
    // std::domain_error for rank-deficient designs.
    //
    // Statistical results assume fixed u coordinates, independent,
    // homoscedastic log-domain errors. No correction for correlated profile
    // samples, imperfect peak selection or the physical signal baseline.
    FoldingLinearFit fitPeakFoldingLinear(
        const std::vector<FoldingProfileView> &profiles,
        FoldingCurvatureMode mode = FoldingCurvatureMode::asymmetric);

    // Experimental, scan-local folding of raw mzML profile spectra.
    // The mzML reader supplies raw m/z and intensity arrays; this routine
    // conservatively selects candidate isolated peaks, folds their measured
    // local coordinates, and calls fitPeakFoldingLinear. It is not qCentroids
    // and must not be used as a validated feature list.
    enum class FoldingSamplingAxis : uint8_t
    {
        tofSqrtMz,        // u proportional to sqrt(m/z)
        orbitrapInvSqrtMz // u proportional to -1/sqrt(m/z)
    };

    struct FoldingSpectrumPeak
    {
        size_t anchorIndex = 0, pointCount = 0;
        double anchorMz = 0, anchorIntensity = 0;
        double apexMz = 0, apexUncertaintyMz = 0;
        bool validApex = false;
        FoldingPeakEstimate coefficients;
    };

    struct FoldingSpectrumFit
    {
        bool fitted = false; // false => no valid grouped model (no fake peaks)
        size_t observedMaxima = 0, selectedCandidates = 0;
        double samplingStep = 0;
        FoldingLinearFit fit;
        std::vector<FoldingSpectrumPeak> peaks;
    };

    // The candidate detector uses a small uninterrupted region around each
    // local maximum, requires both sides and rejects nearby shoulders.
    // Sampling step is estimated from the observed mz array (no filling or
    // interpolation). Results with fewer than two candidate peaks have
    // fitted=false. On malformed unsorted mz data a domain_error is raised.
    FoldingSpectrumFit foldSpectrumProfiles(
        const float *mz, const float *intensity, size_t length,
        FoldingSamplingAxis axis,
        FoldingCurvatureMode mode = FoldingCurvatureMode::asymmetric,
        double minimumApexIntensity = 100.0);

    // mutate b0 so that it is optimal for the exponential case if b1, b2 and b3 are identical

    enum invalid // NOLINT @todo rework causes for invalidation
    {
        ok = 0,
        no_apex = 1,
        invalid_apex = 2,
        no_df = 3,
        invalid_apexToEdge = 4, // this will probably be removed
        f_test_fail = 5,
        invalid_quadratic = 6,
        invalid_area = 7,
        invalid_height = 8,
        invalid_chisq = 9,
        none = -1,
    };

    // ### Feature-specific Code ### //

    size_t findFeatures(const std::vector<EIC> *EICs,
                        const std::vector<float> *convertRT, // correct RT corresponding to every scan number
                        std::vector<FeaturePeak> *res);

#define GLOBAL_MINSCALE 2
#define MINLENGTH 5

#include "../external/qalgorithms_matinverse.h"

#ifndef QALGORITHMS_MAXSCALE_PRECOMPILED
    #error Required array of matrix inverses not defined. Ensure that qalgorithms_matinverse.h exists. You can try regenerating it by going into the "external" directory and running 'cc qalgorithms_matinverse.c -o a.exe && ./a.exe'.
#endif

    const size_t maxscale_global = QALGORITHMS_MAXSCALE_PRECOMPILED;

    // the maximum possible number of apexes that can be present in one continuous section
    // we assume a small number for stack allocation
    const size_t max_apex_per_group = 128;

    // This file contains precomputed values for the inverse matrix transpose (X^T X)^(-1).
    // The precalculation is required since it is a very expensive operation, especially for large X.

    // (X^T X)^(-1) is pre-calculated and stored in the array qalgo_matInverse (included via qalgorithms_matinverse.h).
    // Every array element is accessed by the index of its corresponding scale. This means that, while
    // qalgo_matInverse[0] and qalgo_matInverse[1] are defined, they can never be used to calculate anything
    // because the regression function is only defined starting at scale 2.
    // Only six values of the final matrix are required for the simple case, see below:

    // xT is the transpose of the design matrix X.
    // for scale = 2:
    // xT = | 1  1  1  1  1 |    : all ones
    //      |-2 -1  0  1  2 |    : from -scale to scale
    //      | 4  1  0  0  0 |    : x^2 values for x < 0
    //      | 0  0  0  1  4 |    : x^2 values for x > 0

    // A multi-trace design adds one intercept column per aligned trace, but keeps the same three shared shape columns.

    // For m aligned traces with separate b0 and shared b1,b2,b3, the inverse has seven
    // unique nonzero entries for m >= 2 (counting sign changes as equivalent).
    // A1 = 1/n + (A-1/n)/m, A2 = (A-1/n)/m; B,C,D,E,F are divided by m.
    // example for num_peaks = 2:
    // original matrix with the unique values [a, b, c, d, e, f] (six unique values)
    // | A  0  B  B |
    // | 0  C  D -D |
    // | B  D  E  F |
    // | B -D  F  E |

    // new matrix with the unique values [A1, A2, B, C, D, E, F] (seven unique values)
    // | A1  A2  0  B  B |
    // | A2  A1  0  B  B |
    // | 0   0   C  D -D |
    // | B   B   D  E  F |
    // | B   B  -D  F  E |

    // for num_peaks = 3:
    // new matrix with the unique values [A1, A2, B, C, D, E, F] (the same seven unique values)
    // | A1  A2  A2  0  B  B |
    // | A2  A1  A2  0  B  B |
    // | A2  A2  A1  0  B  B |
    // | 0   0   0   C  D -D |
    // | B   B   B   D  E  F |
    // | B   B   B  -D  F  E |

    // Note that no more than seven different values are needed per scale, even for a multidimensional approach.
    // A1 and A2 are derived from the existing cached A, the scale's n=2*scale+1,
    // and the runtime trace count; no second precomputed cache is necessary.
    inline void peak_model::Model<PeakModelKind::PiecewiseQuadratic>::runningFit(
        const float *y, size_t n, std::vector<Coefficients> *out)
    {
        findCoefficients(y, n, out);
    }
} // namespace qAlgorithms

#endif
