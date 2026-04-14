import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log2

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def cpm_normalize_cython(np.ndarray[DTYPE_t, ndim=2] counts):
    cdef int n_genes = counts.shape[0]
    cdef int n_samples = counts.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] result = np.zeros((n_genes, n_samples), dtype=DTYPE)
    cdef double total
    cdef int i, j

    for j in range(n_samples):
        total = 0.0
        for i in range(n_genes):
            total += counts[i, j]
        for i in range(n_genes):
            result[i, j] = (counts[i, j] / total) * 1e6

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def tmm_normalize_cython(np.ndarray[DTYPE_t, ndim=2] counts):
    cdef int n_genes = counts.shape[0]
    cdef int n_samples = counts.shape[1]
    cdef int i, j, ref_idx
    cdef double ref_total, sample_total, weight_sum, ratio_sum, log_r, w
    cdef double ref_i, sample_i, diff, min_diff

    # 纯C循环计算library sizes
    cdef np.ndarray[DTYPE_t, ndim=1] lib_sizes = np.zeros(n_samples, dtype=DTYPE)
    for j in range(n_samples):
        for i in range(n_genes):
            lib_sizes[j] += counts[i, j]

    # 找参考样本（最接近中位数）
    cdef double median_size = np.median(lib_sizes)
    min_diff = 1e18
    ref_idx = 0
    for j in range(n_samples):
        diff = lib_sizes[j] - median_size
        if diff < 0:
            diff = -diff
        if diff < min_diff:
            min_diff = diff
            ref_idx = j

    ref_total = lib_sizes[ref_idx]

    cdef np.ndarray[DTYPE_t, ndim=1] scaling_factors = np.ones(n_samples, dtype=DTYPE)

    for j in range(n_samples):
        sample_total = lib_sizes[j]
        weight_sum = 0.0
        ratio_sum  = 0.0

        for i in range(n_genes):
            ref_i    = counts[i, ref_idx]
            sample_i = counts[i, j]
            if ref_i > 0 and sample_i > 0:
                # 用C的log2，不用numpy
                log_r = (log2(sample_i / sample_total) -
                         log2(ref_i / ref_total))
                w = 1.0 / (ref_i + sample_i)
                ratio_sum  += w * log_r
                weight_sum += w

        if weight_sum > 0:
            scaling_factors[j] = 2.0 ** (ratio_sum / weight_sum)

    # 应用scaling factors
    cdef np.ndarray[DTYPE_t, ndim=2] result = np.zeros((n_genes, n_samples), dtype=DTYPE)
    for j in range(n_samples):
        for i in range(n_genes):
            result[i, j] = counts[i, j] / scaling_factors[j]

    return result
