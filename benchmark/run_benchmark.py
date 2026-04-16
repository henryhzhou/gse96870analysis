import numpy as np
import time
import json
import sys
sys.path.insert(0, '.')

from normalize_python import cpm_normalize_python, tmm_normalize_python
try:
    from normalize_cython import cpm_normalize_cython, tmm_normalize_cython
    CYTHON_AVAILABLE = True
except ImportError:
    CYTHON_AVAILABLE = False


def benchmark(func, counts, n_runs=15):
    """运行n次，去掉最大最小值，返回trimmed mean和std"""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        _ = func(counts)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    times = sorted(times)[2:-2]  # 去掉最大2个和最小2个
    return np.mean(times), np.std(times)


def run_speedup_benchmark():
    print("=" * 55)
    print("Speedup Benchmark (20000 genes x 45 samples, n=15, trimmed)")
    print("=" * 55)

    np.random.seed(42)
    counts = np.random.negative_binomial(10, 0.3, size=(20000, 45)).astype(np.float64)

    results = {}

    print("\n[ CPM Normalization ]")
    mean_py, std_py = benchmark(cpm_normalize_python, counts)
    print(f"  Pure Python : {mean_py:.4f} ± {std_py:.4f} s")
    results['cpm_python_mean'] = mean_py
    results['cpm_python_std']  = std_py

    if CYTHON_AVAILABLE:
        mean_cy, std_cy = benchmark(cpm_normalize_cython, counts)
        print(f"  Cython      : {mean_cy:.4f} ± {std_cy:.4f} s")
        print(f"  Speedup     : {mean_py/mean_cy:.1f}x")
        results['cpm_cython_mean'] = mean_cy
        results['cpm_cython_std']  = std_cy
        results['cpm_speedup']     = mean_py / mean_cy

    print("\n[ TMM Normalization ]")
    mean_py, std_py = benchmark(tmm_normalize_python, counts)
    print(f"  Pure Python : {mean_py:.4f} ± {std_py:.4f} s")
    results['tmm_python_mean'] = mean_py
    results['tmm_python_std']  = std_py

    if CYTHON_AVAILABLE:
        mean_cy, std_cy = benchmark(tmm_normalize_cython, counts)
        print(f"  Cython      : {mean_cy:.4f} ± {std_cy:.4f} s")
        print(f"  Speedup     : {mean_py/mean_cy:.1f}x")
        results['tmm_cython_mean'] = mean_cy
        results['tmm_cython_std']  = std_cy
        results['tmm_speedup']     = mean_py / mean_cy

    return results


def run_scaling_benchmark():
    print("\n" + "=" * 55)
    print("Scaling Benchmark (varying matrix size, n=10, trimmed)")
    print("=" * 55)

    gene_sizes = [1000, 5000, 10000, 20000, 50000]

    scaling = {
        'gene_sizes': gene_sizes,
        'cpm_python': [], 'cpm_python_std': [],
        'tmm_python': [], 'tmm_python_std': [],
        'cpm_cython': [], 'cpm_cython_std': [],
        'tmm_cython': [], 'tmm_cython_std': [],
    }

    for n_genes in gene_sizes:
        np.random.seed(42)
        counts = np.random.negative_binomial(10, 0.3, size=(n_genes, 45)).astype(np.float64)

        m, s = benchmark(cpm_normalize_python, counts, n_runs=10)
        scaling['cpm_python'].append(m)
        scaling['cpm_python_std'].append(s)

        m, s = benchmark(tmm_normalize_python, counts, n_runs=10)
        scaling['tmm_python'].append(m)
        scaling['tmm_python_std'].append(s)

        if CYTHON_AVAILABLE:
            m, s = benchmark(cpm_normalize_cython, counts, n_runs=10)
            scaling['cpm_cython'].append(m)
            scaling['cpm_cython_std'].append(s)

            m, s = benchmark(tmm_normalize_cython, counts, n_runs=10)
            scaling['tmm_cython'].append(m)
            scaling['tmm_cython_std'].append(s)

        print(f"  {n_genes:6d} genes: CPM_py={scaling['cpm_python'][-1]:.4f}s  TMM_py={scaling['tmm_python'][-1]:.4f}s")

    return scaling


if __name__ == "__main__":
    speedup = run_speedup_benchmark()
    scaling = run_scaling_benchmark()

    with open("benchmark_results.json", "w") as f:
        json.dump({'speedup': speedup, 'scaling': scaling}, f, indent=2)
    print("\nResults saved to benchmark_results.json")
