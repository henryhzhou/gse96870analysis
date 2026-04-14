# ============================================================
# plot_benchmark.py — 生成benchmark可视化图
# ============================================================

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 读取benchmark结果
with open("benchmark_results.json") as f:
    results = json.load(f)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Pure Python vs Cython Benchmark\n(20,000 genes × 45 samples)",
             fontsize=14, fontweight='bold')

colors = {'Python': '#E74C3C', 'Cython': '#2ECC71'}

for ax, (method, key_py, key_cy) in zip(
    axes,
    [("CPM Normalization", "cpm_python", "cpm_cython"),
     ("TMM Normalization", "tmm_python", "tmm_cython")]
):
    py_time = results.get(key_py, 0)
    cy_time = results.get(key_cy, 0)
    speedup = results.get(key_py.replace("python","speedup").replace("cpm_","cpm_").replace("tmm_","tmm_"), 0)

    bars = ax.bar(['Pure Python', 'Cython'],
                  [py_time, cy_time],
                  color=[colors['Python'], colors['Cython']],
                  width=0.5, edgecolor='black', linewidth=0.8)

    # 在bar上标注时间
    for bar, t in zip(bars, [py_time, cy_time]):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 0.01 * max(py_time, cy_time),
                f'{t:.3f}s', ha='center', va='bottom', fontweight='bold')

    # 标注Speedup
    if cy_time > 0:
        ax.text(0.5, 0.95, f'Speedup: {py_time/cy_time:.1f}×',
                transform=ax.transAxes, ha='center', va='top',
                fontsize=13, color='#2C3E50', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#F9F9F9', edgecolor='grey'))

    ax.set_title(method, fontsize=12)
    ax.set_ylabel("Execution Time (s)")
    ax.set_ylim(0, max(py_time, cy_time) * 1.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("benchmark_comparison.pdf", dpi=150, bbox_inches='tight')
plt.savefig("benchmark_comparison.png", dpi=150, bbox_inches='tight')
print("图已保存: benchmark_comparison.pdf / .png")
