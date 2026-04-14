# ============================================================
# run_benchmark.py — 对比pure Python vs Cython性能
# ============================================================

import numpy as np
import time
import json
import sys
sys.path.insert(0, '.')

from normalize_python import cpm_normalize_python, tmm_normalize_python

# 导入Cython版本（需要先编译）
try:
    from normalize_cython import cpm_normalize_cython, tmm_normalize_cython
    CYTHON_AVAILABLE = True
except ImportError:
    print("Cython版本未编译，请先运行: python setup.py build_ext --inplace")
    CYTHON_AVAILABLE = False


def benchmark(func, counts, n_runs=10):
    """运行n次取平均时间"""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = func(counts)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return np.mean(times), np.std(times), result


def run_all_benchmarks():
    # 模拟真实数据大小：20000基因 × 45样本
    np.random.seed(42)
    counts = np.random.negative_binomial(10, 0.3, size=(20000, 45)).astype(np.float64)
    
    print("=" * 55)
    print("Benchmark: Pure Python vs Cython")
    print("矩阵大小: 20000基因 × 45样本")
    print("=" * 55)
    
    results = {}

    # ── CPM ──────────────────────────────────────────
    print("\n[ CPM Normalization ]")
    
    mean_py, std_py, out_py = benchmark(cpm_normalize_python, counts)
    print(f"  Pure Python : {mean_py:.4f} ± {std_py:.4f} 秒")
    results['cpm_python'] = mean_py

    if CYTHON_AVAILABLE:
        mean_cy, std_cy, out_cy = benchmark(cpm_normalize_cython, counts)
        print(f"  Cython      : {mean_cy:.4f} ± {std_cy:.4f} 秒")
        print(f"  加速比      : {mean_py / mean_cy:.1f}x")
        results['cpm_cython'] = mean_cy
        results['cpm_speedup'] = mean_py / mean_cy

        # 验证结果一致性
        if np.allclose(out_py, out_cy, rtol=1e-5):
            print("  结果验证   : ✓ 一致")
        else:
            print("  结果验证   : ✗ 不一致！")

    # ── TMM ──────────────────────────────────────────
    print("\n[ TMM Normalization ]")
    
    mean_py, std_py, out_py = benchmark(tmm_normalize_python, counts)
    print(f"  Pure Python : {mean_py:.4f} ± {std_py:.4f} 秒")
    results['tmm_python'] = mean_py

    if CYTHON_AVAILABLE:
        mean_cy, std_cy, out_cy = benchmark(tmm_normalize_cython, counts)
        print(f"  Cython      : {mean_cy:.4f} ± {std_cy:.4f} 秒")
        print(f"  加速比      : {mean_py / mean_cy:.1f}x")
        results['tmm_cython'] = mean_cy
        results['tmm_speedup'] = mean_py / mean_cy

    print("\n" + "=" * 55)
    
    # 保存结果供画图用
    with open("benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("结果已保存到 benchmark_results.json")
    
    return results


if __name__ == "__main__":
    run_all_benchmarks()
