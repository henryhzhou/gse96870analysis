# ============================================================
# plot_benchmark.py — Speedup + Scaling plots
# ============================================================

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

with open("benchmark_results.json") as f:
    data = json.load(f)

sp = data['speedup']
sc = data['scaling']

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

colors = {'Python': '#E24B4A', 'Cython': '#1D9E75'}

# ── Plot 1: CPM Speedup with error bars ───────────────────
ax1 = fig.add_subplot(gs[0, 0])
methods = ['Pure Python', 'Cython']
means   = [sp['cpm_python_mean'], sp.get('cpm_cython_mean', 0)]
stds    = [sp['cpm_python_std'],  sp.get('cpm_cython_std', 0)]
bars = ax1.bar(methods, means, color=[colors['Python'], colors['Cython']],
               width=0.5, edgecolor='black', linewidth=0.8)
ax1.errorbar(methods, means, yerr=stds, fmt='none', color='black',
             capsize=6, capthick=1.5, elinewidth=1.5)
for bar, m in zip(bars, means):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(stds)*1.1,
             f'{m:.3f}s', ha='center', va='bottom', fontweight='bold', fontsize=10)
speedup = sp.get('cpm_speedup', 0)
ax1.text(0.5, 0.92, f'Speedup: {speedup:.1f}×', transform=ax1.transAxes,
         ha='center', fontsize=12, fontweight='bold', color='#2C2C2A',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#F1EFE8', edgecolor='grey'))
ax1.set_title('CPM Normalization — Speedup', fontweight='bold')
ax1.set_ylabel('Execution time (s)')
ax1.set_ylim(0, max(means) * 1.4)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.text(0.5, -0.15, 'n=10 runs, mean ± std', transform=ax1.transAxes,
         ha='center', fontsize=9, color='gray', style='italic')

# ── Plot 2: TMM Speedup with error bars ───────────────────
ax2 = fig.add_subplot(gs[0, 1])
means2 = [sp['tmm_python_mean'], sp.get('tmm_cython_mean', 0)]
stds2  = [sp['tmm_python_std'],  sp.get('tmm_cython_std', 0)]
bars2 = ax2.bar(methods, means2, color=[colors['Python'], colors['Cython']],
                width=0.5, edgecolor='black', linewidth=0.8)
ax2.errorbar(methods, means2, yerr=stds2, fmt='none', color='black',
             capsize=6, capthick=1.5, elinewidth=1.5)
for bar, m in zip(bars2, means2):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.05,
             f'{m:.3f}s', ha='center', va='bottom', fontweight='bold', fontsize=10)
speedup2 = sp.get('tmm_speedup', 0)
ax2.text(0.5, 0.92, f'Speedup: {speedup2:.1f}×', transform=ax2.transAxes,
         ha='center', fontsize=12, fontweight='bold', color='#2C2C2A',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#F1EFE8', edgecolor='grey'))
ax2.set_title('TMM Normalization — Speedup', fontweight='bold')
ax2.set_ylabel('Execution time (s)')
ax2.set_ylim(0, max(means2) * 1.4)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.text(0.5, -0.15, 'n=10 runs, mean ± std', transform=ax2.transAxes,
         ha='center', fontsize=9, color='gray', style='italic')

# ── Plot 3: CPM Scaling ────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
gene_sizes = sc['gene_sizes']
ax3.errorbar(gene_sizes, sc['cpm_python'], yerr=sc['cpm_python_std'],
             label='Pure Python', color=colors['Python'], marker='o',
             linewidth=2, capsize=4, markersize=6)
if 'cpm_cython' in sc:
    ax3.errorbar(gene_sizes, sc['cpm_cython'], yerr=sc['cpm_cython_std'],
                 label='Cython', color=colors['Cython'], marker='s',
                 linewidth=2, capsize=4, markersize=6, linestyle='--')
ax3.set_title('CPM Normalization — Scaling', fontweight='bold')
ax3.set_xlabel('Number of genes')
ax3.set_ylabel('Execution time (s)')
ax3.legend(fontsize=10)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k'))
ax3.text(0.5, -0.18, 'n=5 runs per size, 45 samples fixed', transform=ax3.transAxes,
         ha='center', fontsize=9, color='gray', style='italic')

# ── Plot 4: TMM Scaling ────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
ax4.errorbar(gene_sizes, sc['tmm_python'], yerr=sc['tmm_python_std'],
             label='Pure Python', color=colors['Python'], marker='o',
             linewidth=2, capsize=4, markersize=6)
if 'tmm_cython' in sc:
    ax4.errorbar(gene_sizes, sc['tmm_cython'], yerr=sc['tmm_cython_std'],
                 label='Cython', color=colors['Cython'], marker='s',
                 linewidth=2, capsize=4, markersize=6, linestyle='--')
ax4.set_title('TMM Normalization — Scaling', fontweight='bold')
ax4.set_xlabel('Number of genes')
ax4.set_ylabel('Execution time (s)')
ax4.legend(fontsize=10)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k'))
ax4.text(0.5, -0.18, 'n=5 runs per size, 45 samples fixed', transform=ax4.transAxes,
         ha='center', fontsize=9, color='gray', style='italic')

fig.suptitle('Pure Python vs Cython: Speedup & Scaling Analysis\n(RNA-seq count matrix normalization)',
             fontsize=13, fontweight='bold', y=1.01)

plt.savefig('benchmark_full.png', dpi=150, bbox_inches='tight')
plt.savefig('benchmark_full.pdf', dpi=150, bbox_inches='tight')
print("Saved: benchmark_full.png / .pdf")
