#! python
# -*- coding:utf-8 -*-
###
# --------------------------------------------------------------------------------
# 文件名: test.py
# 创建时间: 2026-02-17 16:59:15 Tue
# 说明:
# 作者: Calibur88
# 主机: LAPTOP-D92A7OL2
# --------------------------------------------------------------------------------
# 最后编辑作者: Calibur88
# 最后修改时间: 2026-02-17 17:40:23 Tue
# --------------------------------------------------------------------------------
# Copyright (c) 2026 Calibur88
# --------------------------------------------------------------------------------
# 更新历史:
# --------------------------------------------------------------------------------
# 时间      		作者		信息
# ----------		---		------------------------------------------------------
###
"""
GUE Polished Riemann Zero Generation
基于希尔伯特-波利亚猜想的物理实现：前向迭代 + GUE谱波利色算符
"""

import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def read_zeros(filename):
    """读取黎曼零点数据文件"""
    zeros = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                try:
                    zeros.append(float(line))
                except ValueError:
                    continue
    return np.array(zeros)


def forward_iteration(seed, n_steps):
    """
    前向迭代公式（von Mangoldt渐近）：
    γ_{n+1} = γ_n + 2π / ln(γ_n / 2π)
    """
    zeros = [seed]
    for _ in range(n_steps - 1):
        gamma_n = zeros[-1]
        gamma_next = gamma_n + 2 * np.pi / np.log(gamma_n / (2 * np.pi))
        zeros.append(gamma_next)
    return np.array(zeros)


def lambert_w_seed(n):
    """
    Lambert W函数计算第n个零点的渐近近似值：
    γ_n ≈ 2π(n - 3/8) / W((n - 3/8)/e)
    """
    m = n - 3 / 8
    W_val = lambertw(m / np.e).real
    return 2 * np.pi * m / W_val


def polish_operator_gue(zeros_bare, seed=42):
    """
    波利色算符：GUE随机矩阵谱映射

    构造N×N高斯酉系综（GUE）哈密顿量 H = (X + X†)/√(2N)，
    提取本征值谱并线性映射至裸谱区间，保持GUE统计特性（间距比⟨r⟩≈0.602）。
    """
    np.random.seed(seed)
    n = len(zeros_bare)

    # GUE随机矩阵构造
    X = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
    H = (X + X.conj().T) / np.sqrt(2 * n)

    # 本征值提取与排序
    eigvals = np.sort(np.real(np.linalg.eigvals(H)))

    # 线性映射至裸谱数值范围
    gue_min, gue_max = eigvals[0], eigvals[-1]
    bare_min, bare_max = zeros_bare[0], zeros_bare[-1]

    if abs(gue_max - gue_min) < 1e-10:
        return zeros_bare

    scale = (bare_max - bare_min) / (gue_max - gue_min)
    offset = bare_min - scale * gue_min

    zeros_phys = scale * eigvals + offset

    return zeros_phys


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

    # 间距比统计
    r_real, ratios_real = compute_spacing_ratio(zeros_real)
    r_gen, ratios_gen = compute_spacing_ratio(zeros_gen)

    print(f"\n间距比统计:")
    print(f"  理论GUE ⟨r⟩:     0.602")
    print(f"  真实零点 ⟨r⟩:    {r_real:.4f}")
    print(f"  生成零点 ⟨r⟩:    {r_gen:.4f}")
    print(f"  与GUE理论偏差:   {abs(r_gen - 0.602):.4f}")

    # 平均间距
    mean_spacing_real = np.mean(np.diff(zeros_real))
    mean_spacing_gen = np.mean(np.diff(zeros_gen))
    print(f"\n平均间距:")
    print(f"  真实: {mean_spacing_real:.4f}")
    print(f"  生成: {mean_spacing_gen:.4f}")
    print(
        f"  相对误差: {abs(mean_spacing_gen - mean_spacing_real)/mean_spacing_real*100:.2f}%"
    )

    # 区间分布
    counts_real, counts_gen, bin_edges = interval_statistics(zeros_real, zeros_gen)

    print(f"\n区间分布统计:")
    for i in range(10):
        print(
            f"  [{bin_edges[i]:.1f}, {bin_edges[i+1]:.1f}): 真实={counts_real[i]:3d}, 生成={counts_gen[i]:3d}"
        )

    mse = np.mean((counts_real - counts_gen) ** 2)
    print(f"\n密度分布MSE: {mse:.2f}")

    return {
        "r_real": r_real,
        "r_gen": r_gen,
        "ratios_real": ratios_real,
        "ratios_gen": ratios_gen,
        "counts_real": counts_real,
        "counts_gen": counts_gen,
        "bin_edges": bin_edges,
    }


def plot_comparison(stats_list, interval_names, filename="zeros_comparison.png"):
    """三区间对比可视化"""
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    fig.suptitle(
        "Riemann Zeros: Bare Spectrum + GUE Polish Operator",
        fontsize=14,
        fontweight="bold",
    )

    for idx, (name, stats) in enumerate(zip(interval_names, stats_list)):
        # 区间密度分布
        ax1 = axes[idx, 0]
        ax1.set_title(f"{name}\nInterval Distribution", fontsize=10)
        x_pos = np.arange(10)
        width = 0.35
        ax1.bar(
            x_pos - width / 2,
            stats["counts_real"],
            width,
            label="True Zeros",
            color="blue",
            alpha=0.7,
            edgecolor="black",
        )
        ax1.bar(
            x_pos + width / 2,
            stats["counts_gen"],
            width,
            label="GUE Polished",
            color="red",
            alpha=0.7,
            edgecolor="black",
        )
        ax1.set_xlabel("Bin Index")
        ax1.set_ylabel("Count")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 间距比分布直方图
        ax2 = axes[idx, 1]
        ax2.set_title("Spacing Ratio Distribution", fontsize=10)
        bins = np.linspace(0, 1, 21)
        if stats["ratios_real"]:
            ax2.hist(
                stats["ratios_real"],
                bins=bins,
                alpha=0.6,
                label=f'True ⟨r⟩={stats["r_real"]:.3f}',
                color="blue",
                density=True,
            )
        if stats["ratios_gen"]:
            ax2.hist(
                stats["ratios_gen"],
                bins=bins,
                alpha=0.6,
                label=f'GUE ⟨r⟩={stats["r_gen"]:.3f}',
                color="red",
                density=True,
            )
        ax2.axvline(
            0.602, color="green", linestyle="--", linewidth=2, label="GUE Theory 0.602"
        )
        ax2.axvline(
            0.386, color="orange", linestyle=":", linewidth=2, label="Poisson 0.386"
        )
        ax2.set_xlabel("Spacing Ratio r")
        ax2.set_ylabel("Probability Density")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        # 间距比统计对比
        ax3 = axes[idx, 2]
        ax3.set_title("Mean Spacing Ratio", fontsize=10)
        categories = ["True", "GUE Gen.", "GUE Theory"]
        values = [stats["r_real"], stats["r_gen"], 0.602]
        colors = ["blue", "red", "green"]
        bars = ax3.bar(categories, values, color=colors, alpha=0.7)
        ax3.set_ylabel("⟨r⟩")
        ax3.set_ylim(0, 1)
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax3.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{val:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )
        ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved: {filename}")


def main():
    """主程序：三区间零点生成与统计验证"""
    print("Loading Riemann zeros data...")
    try:
        zeros_real_all = read_zeros("zeros1.txt")
        print(f"Loaded {len(zeros_real_all)} zeros")
    except FileNotFoundError:
        print("Error: zeros1.txt not found")
        return

    zeros_real_all = np.sort(zeros_real_all)
    results = []
    interval_names = []

    # 区间1: 第1-100个零点
    print("\nGenerating zeros 1-100...")
    seed_1 = 14.134725142
    zeros_gen_1_bare = forward_iteration(seed_1, 100)
    zeros_gen_1 = polish_operator_gue(zeros_gen_1_bare, seed=42)
    zeros_real_1 = zeros_real_all[0:100]
    stats1 = analyze_interval(zeros_real_1, zeros_gen_1, "Zeros 1-100")
    results.append(stats1)
    interval_names.append("1-100")

    # 区间2: 第500-1000个零点
    print("\nGenerating zeros 500-1000...")
    seed_500 = lambert_w_seed(500)
    zeros_gen_2_bare = forward_iteration(seed_500, 501)
    zeros_gen_2 = polish_operator_gue(zeros_gen_2_bare, seed=123)
    zeros_real_2 = zeros_real_all[499:1000]
    stats2 = analyze_interval(zeros_real_2, zeros_gen_2, "Zeros 500-1000")
    results.append(stats2)
    interval_names.append("500-1000")

    # 区间3: 第1000-2000个零点
    print("\nGenerating zeros 1000-2000...")
    seed_1000 = lambert_w_seed(1000)
    zeros_gen_3_bare = forward_iteration(seed_1000, 1001)
    zeros_gen_3 = polish_operator_gue(zeros_gen_3_bare, seed=456)
    zeros_real_3 = zeros_real_all[999:2000]
    stats3 = analyze_interval(zeros_real_3, zeros_gen_3, "Zeros 1000-2000")
    results.append(stats3)
    interval_names.append("1000-2000")

    # 可视化
    plot_comparison(results, interval_names)

    # 样本输出
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
            print(
                f"  n={i+1:4d}: True={zr[i]:12.6f}, GUE={zg[i]:12.6f}, "
                f"err={abs(zr[i]-zg[i]):.4f}"
            )


if __name__ == "__main__":
    main()
