# ============================================================
# normalize_python.py — 纯Python实现的CPM normalization
# 用于和Cython版本做benchmark对比
# ============================================================

import numpy as np
import time

def cpm_normalize_python(counts):
    """
    CPM normalization: Counts Per Million
    公式: CPM = (count / total_counts_per_sample) * 1e6
    纯Python + numpy实现，没有任何编译优化
    """
    n_genes, n_samples = counts.shape
    result = np.zeros_like(counts, dtype=np.float64)
    
    # 故意用循环而不是向量化，方便和Cython对比
    for j in range(n_samples):
        total = 0.0
        for i in range(n_genes):
            total += counts[i, j]
        for i in range(n_genes):
            result[i, j] = (counts[i, j] / total) * 1e6
    
    return result


def tmm_normalize_python(counts):
    """
    TMM normalization (简化版)
    计算每个样本相对于参考样本的scaling factor
    """
    n_genes, n_samples = counts.shape
    
    # 选参考样本（总count最接近中位数的样本）
    lib_sizes = np.sum(counts, axis=0)
    ref_idx = np.argmin(np.abs(lib_sizes - np.median(lib_sizes)))
    ref = counts[:, ref_idx].astype(np.float64)
    ref_total = np.sum(ref)
    
    scaling_factors = np.ones(n_samples)
    
    for j in range(n_samples):
        sample = counts[:, j].astype(np.float64)
        sample_total = np.sum(sample)
        
        # 过滤低表达基因
        valid = (ref > 0) & (sample > 0)
        ref_v    = ref[valid]
        sample_v = sample[valid]
        
        # 计算log ratio
        log_ratio = np.log2(sample_v / sample_total) - np.log2(ref_v / ref_total)
        
        # 加权平均（TMM核心）
        weight = (1.0 / (ref_v + sample_v))  # 简化权重
        scaling_factors[j] = 2 ** np.sum(weight * log_ratio) / np.sum(weight)
    
    # 应用scaling factors
    result = counts.astype(np.float64) / scaling_factors[np.newaxis, :]
    return result


if __name__ == "__main__":
    # 模拟45样本 × 20000基因的count矩阵
    np.random.seed(42)
    counts = np.random.negative_binomial(10, 0.3, size=(20000, 45)).astype(np.float64)
    
    print("矩阵大小: {} 基因 × {} 样本".format(counts.shape[0], counts.shape[1]))
    print()
    
    # Benchmark CPM
    print("CPM Normalization (Pure Python):")
    times = []
    for i in range(3):
        start = time.time()
        result = cpm_normalize_python(counts)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  第{i+1}次: {elapsed:.4f} 秒")
    print(f"  平均: {np.mean(times):.4f} 秒")
    print()
    
    # Benchmark TMM
    print("TMM Normalization (Pure Python):")
    times = []
    for i in range(3):
        start = time.time()
        result = tmm_normalize_python(counts)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  第{i+1}次: {elapsed:.4f} 秒")
    print(f"  平均: {np.mean(times):.4f} 秒")
