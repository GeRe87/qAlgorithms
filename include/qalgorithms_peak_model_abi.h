#ifndef QALGORITHMS_PEAK_MODEL_ABI_H
#define QALGORITHMS_PEAK_MODEL_ABI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QALGORITHMS_PEAK_MODEL_ABI_V1 1u
#define QALGORITHMS_PEAK_MODEL_ENTRY_V1 "qalgorithms_peak_model_v1"

/* Function results. The host owns all output buffers. */
#define QALGORITHMS_PEAK_FAILED 0
#define QALGORITHMS_PEAK_PASSED 1
#define QALGORITHMS_PEAK_UNAVAILABLE 2

typedef struct QAlgorithmsPeakNoiseV1
{
    double log_mse;
    uint64_t residual_df;
    double alpha;
} QAlgorithmsPeakNoiseV1;

typedef struct QAlgorithmsPeakAreaV1
{
    int32_t shape_certified;
    double log_lower_bound;
    double left;
    double right;
} QAlgorithmsPeakAreaV1;

typedef struct QAlgorithmsPeakCoverageV1
{
    double fraction;
    double standard_error;
    double lower_confidence_bound;
} QAlgorithmsPeakCoverageV1;

/* Export qalgorithms_peak_model_v1 with C linkage and default visibility.
 * Its descriptor and function pointers live until the module is unloaded.
 * Each fit is allocated and destroyed in the same module. A non-null fit
 * returned by fit_window is destroyed even if that call reports failure.
 * Log observations and other host buffers are borrowed for each call only.
 * Callbacks must not throw across this boundary. Sample i has coordinate
 * (i - floor(length/2)) * step. Area uses intensity times that coordinate's
 * unit. The same fit procedure is replayed in background rotations. Within
 * one validation the host destroys each fit before creating the next.
 * fit_window fits the supplied log observations by OLS. parameter_count is
 * the fitted rank used to compute residual degrees of freedom.
 */
typedef struct QAlgorithmsPeakModelV1
{
    uint32_t abi_version;
    uint32_t struct_size;
    const char *name;
    uint32_t parameter_count;
    int32_t (*fit_window)(const double *logs, uint64_t length, double step, void **fit);
    void (*destroy_fit)(void *fit);
    int32_t (*has_peak_shape)(const void *fit);
    int32_t (*apex)(const void *fit, double *position);
    int32_t (*predict_log)(const void *fit, double coordinate, double *value);
    int32_t (*form_confidence)(const void *fit, const QAlgorithmsPeakNoiseV1 *noise);
    int32_t (*area_evidence)(const void *fit, const QAlgorithmsPeakNoiseV1 *noise,
        double minimum_area, QAlgorithmsPeakAreaV1 *result);
    int32_t (*window_coverage)(const void *fit, const QAlgorithmsPeakNoiseV1 *noise,
        double left, double right, double minimum_fraction, QAlgorithmsPeakCoverageV1 *result);
} QAlgorithmsPeakModelV1;

typedef const QAlgorithmsPeakModelV1 *(*QAlgorithmsPeakModelEntryV1)(void);

#ifdef __cplusplus
}
#endif

#endif
