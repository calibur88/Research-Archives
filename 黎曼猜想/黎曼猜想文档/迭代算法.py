#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Riemann Zeros Fast Generation
基于 Lambert W 主项公式的直接计算算法
T_n ≈ 2π(n-7/8) / W_0((n-7/8)/e)
"""

import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

try:
    import mpmath as mp
    mp.mp.dps = 50
except ImportError:
    print("请安装 mpmath: pip install mpmath")
    exit()


def lambert_w_seed(n, c=7/8):
    """Lambert W 主项：γ_n ≈ 2π(n-c) / W_0((n-c)/e)"""
    m = n - c
    W_val = lambertw(m / np.e).real
    return 2 * np.pi * m / W_val


def generate_zeros_lambert(start, end):
    """使用 Lambert W 主项公式直接计算第 start 到第 end 个零点虚部"""
    zeros = []
    for n in range(start, end + 1):
        zeros.append(lambert_w_seed(n))
    return np.array(zeros)


def get_true_zeros_mpmath(start, end):
    """使用 mpmath 获取真实零点虚部"""
    zeros = []
    for n in range(start, end + 1):
        z = mp.zetazero(n)
        zeros.append(float(z.imag))
    return np.array(zeros)


def compute_spacing_ratio(zeros):
    """计算相邻能级间距比 r = min(s_n, s_{n+1}) / max(s_n, s_{n+1})"""
    if len(zeros) < 3:
        return 0.0, []
    spacings = np.diff(zeros)
    ratios = []
    for i in range(len(spacings) - 1):
        s1, s2 = spacings[i], spacings[i + 1]
        if max(s1, s2) > 1e-10:
            r = min(s1, s2) / max(s1, s2)
            ratios.append(r)
    mean_r = np.mean(ratios) if ratios else 0.0
    return mean_r, ratios


def interval_statistics(zeros_real, zeros_gen):
    """区间分布统计"""
    min_val = min(zeros_real[0], zeros_gen[0])
    max_val = max(zeros_real[-1], zeros_gen[-1])
    bin_edges = np.linspace(min_val, max_val, 11)
    counts_real, _ = np.histogram(zeros_real, bins=bin_edges)
    counts_gen, _ = np.histogram(zeros_gen, bins=bin_edges)
    return counts_real, counts_gen, bin_edges


def analyze_interval(zeros_real, zeros_gen, interval_name):
    """区间统计分析与输出"""
    print(f"\n{'='*70}")
    print(f"区间分析: {interval_name}")
    print(f"{'='*70}")
    
    print(f"真实零点数量: {len(zeros_real)}")
    print(f"生成零点数量: {len(zeros_gen)}")
    
    print(f"\n数值范围:")
    print(f"  真实: [{zeros_real[0]:.3f}, {zeros_real[-1]:.3f}]")
    print(f"  生成: [{zeros_gen[0]:.3f}, {zeros_gen[-1]:.3f}]")
    
    r_real, ratios_real = compute_spacing_ratio(zeros_real)
    r_gen, ratios_gen = compute_spacing_ratio(zeros_gen)
    
    print(f"\n间距比统计:")
    print(f"  理论GUE ⟨r⟩:     0.602")
    print(f"  真实零点 ⟨r⟩:    {r_real:.4f}")
    print(f"  生成零点 ⟨r⟩:    {r_gen:.4f}")
    print(f"  与GUE理论偏差:   {abs(r_gen - 0.602):.4f}")
    
    mean_spacing_real = np.mean(np.diff(zeros_real))
    mean_spacing_gen = np.mean(np.diff(zeros_gen))
    print(f"\n平均间距:")
    print(f"  真实: {mean_spacing_real:.4f}")
    print(f"  生成: {mean_spacing_gen:.4f}")
    print(f"  相对误差: {abs(mean_spacing_gen - mean_spacing_real)/mean_spacing_real*100:.2f}%")
    
    counts_real, counts_gen, bin_edges = interval_statistics(zeros_real, zeros_gen)
    print(f"\n区间分布统计:")
    for i in range(10):
        print(f"  [{bin_edges[i]:.1f}, {bin_edges[i+1]:.1f}): 真实={counts_real[i]:3d}, 生成={counts_gen[i]:3d}")
    
    mse = np.mean((counts_real - counts_gen) ** 2)
    print(f"\n密度分布MSE: {mse:.2f}")
    
    errors = np.abs(zeros_gen - zeros_real)
    print(f"\n数值精度:")
    print(f"  平均绝对误差: {np.mean(errors):.6f}")
    print(f"  最大绝对误差: {np.max(errors):.6f}")
    print(f"  平均相对误差: {np.mean(errors/zeros_real)*100:.4f}%")
    
    return {
        "r_real": r_real,
        "r_gen": r_gen,
        "ratios_real": ratios_real,
        "ratios_gen": ratios_gen,
        "counts_real": counts_real,
        "counts_gen": counts_gen,
        "bin_edges": bin_edges,
        "zeros_real": zeros_real,
        "zeros_gen": zeros_gen,
        "errors": errors,
    }


def plot_comparison(stats_list, interval_names, filename="zeros_comparison.png"):
    """三区间对比可视化"""
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    fig.suptitle("Riemann Zeros: Lambert W Formula vs mpmath True Zeros", fontsize=14, fontweight="bold")
    
    for idx, (name, stats) in enumerate(zip(interval_names, stats_list)):
        ax1 = axes[idx, 0]
        ax1.set_title(f"{name}\nInterval Distribution", fontsize=10)
        x_pos = np.arange(10)
        width = 0.35
        ax1.bar(x_pos - width/2, stats["counts_real"], width, label="True Zeros", color="blue", alpha=0.7, edgecolor="black")
        ax1.bar(x_pos + width/2, stats["counts_gen"], width, label="Lambert W", color="red", alpha=0.7, edgecolor="black")
        ax1.set_xlabel("Bin Index")
        ax1.set_ylabel("Count")
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        ax2 = axes[idx, 1]
        ax2.set_title("Spacing Ratio Distribution", fontsize=10)
        bins = np.linspace(0, 1, 21)
        if stats["ratios_real"]:
            ax2.hist(stats["ratios_real"], bins=bins, alpha=0.6, label=f'True ⟨r⟩={stats["r_real"]:.3f}', color="blue", density=True)
        if stats["ratios_gen"]:
            ax2.hist(stats["ratios_gen"], bins=bins, alpha=0.6, label=f'Lambert ⟨r⟩={stats["r_gen"]:.3f}', color="red", density=True)
        ax2.axvline(0.602, color="green", linestyle="--", linewidth=2, label="GUE Theory 0.602")
        ax2.axvline(0.386, color="orange", linestyle=":", linewidth=2, label="Poisson 0.386")
        ax2.set_xlabel("Spacing Ratio r")
        ax2.set_ylabel("Probability Density")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
        
        ax3 = axes[idx, 2]
        ax3.set_title("Absolute Error", fontsize=10)
        n_start = int(name.split('-')[0])
        n_vals = np.arange(len(stats["errors"])) + n_start
        ax3.plot(n_vals, stats["errors"], 'g-', lw=1)
        ax3.set_xlabel("n")
        ax3.set_ylabel("Absolute Error")
        ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved: {filename}")


def main():
    print("Loading true zeros from mpmath...")
    print("=" * 70)
    
    results = []
    interval_names = []
    
    print("\nGenerating zeros 1-100 via Lambert W formula...")
    zeros_gen_1 = generate_zeros_lambert(1, 100)
    zeros_real_1 = get_true_zeros_mpmath(1, 100)
    stats1 = analyze_interval(zeros_real_1, zeros_gen_1, "Zeros 1-100")
    results.append(stats1)
    interval_names.append("1-100")
    
    print("\nGenerating zeros 500-1000 via Lambert W formula...")
    zeros_gen_2 = generate_zeros_lambert(500, 1000)
    zeros_real_2 = get_true_zeros_mpmath(500, 1000)
    stats2 = analyze_interval(zeros_real_2, zeros_gen_2, "Zeros 500-1000")
    results.append(stats2)
    interval_names.append("500-1000")
    
    print("\nGenerating zeros 1000-2000 via Lambert W formula...")
    zeros_gen_3 = generate_zeros_lambert(1000, 2000)
    zeros_real_3 = get_true_zeros_mpmath(1000, 2000)
    stats3 = analyze_interval(zeros_real_3, zeros_gen_3, "Zeros 1000-2000")
    results.append(stats3)
    interval_names.append("1000-2000")
    
    plot_comparison(results, interval_names)
    
    print("\n" + "=" * 70)
    print("Sample Comparison:")
    print("=" * 70)
    for name, zr, zg in [
        ("1-100", zeros_real_1, zeros_gen_1),
        ("500-1000", zeros_real_2, zeros_gen_2),
        ("1000-2000", zeros_real_3, zeros_gen_3),
    ]:
        print(f"\n{name}:")
        indices = [0, len(zr) // 2, -1]
        for i in indices:
            print(f"  n={i+1:4d}: True={zr[i]:12.6f}, Lambert={zg[i]:12.6f}, err={abs(zr[i]-zg[i]):.4f}")


if __name__ == "__main__":
    main()