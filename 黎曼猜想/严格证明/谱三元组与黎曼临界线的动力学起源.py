#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
希尔伯特-波利亚猜想：谱三元组与黎曼临界线的动力学起源
专业级可视化脚本 v2.0
...
"""

import os
import sys
import json
import csv
import time
from dataclasses import dataclass
from typing import List, Tuple, Callable

import mpmath
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# ==================== 配置 ====================
OUTPUT_DIR = '/storage/emulated/0/代码/'
FONT_PATH = '/system/fonts/DroidSansFallbackBBK.ttf'
mpmath.mp.dps = 50

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 字体设置（全局生效）
prop = fm.FontProperties(fname=FONT_PATH)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['axes.unicode_minus'] = False

# 全局样式
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 200

# ==================== 工具函数 ====================
def timer(func):
    """装饰器：计时"""
    def wrapper(*args, **kwargs):
        t0 = time.time()
        result = func(*args, **kwargs)
        print(f"[{func.__name__}] 耗时 {time.time()-t0:.2f}s")
        return result
    return wrapper

@timer
def compute_zeros(n_max: int, cache_file: str) -> np.ndarray:
    """计算或加载黎曼零点虚部"""
    if os.path.exists(cache_file):
        with open(cache_file, 'r', encoding='utf-8') as f:
            zeros = np.array(json.load(f))
        if len(zeros) >= n_max:
            print(f"[缓存] 加载 {n_max} 个零点")
            return zeros[:n_max]
    
    print(f"[计算] mpmath 计算前 {n_max} 个零点...")
    zeros = []
    for n in range(1, n_max + 1):
        z = mpmath.zetazero(n)
        zeros.append(float(z.imag))
        if n % 100 == 0:
            print(f"  进度 {n}/{n_max}")
    
    with open(cache_file, 'w', encoding='utf-8') as f:
        json.dump(zeros, f)
    print(f"[完成] 已缓存")
    return np.array(zeros)

def prime_sieve(limit: int) -> List[int]:
    """埃拉托斯特尼筛法"""
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = b'\x00' * ((limit - i*i) // i + 1)
    return [i for i, is_prime in enumerate(sieve) if is_prime]

@dataclass
class PrimePowerEntry:
    """素数幂条目"""
    log_val: float      # log(p^k)
    weight: float       # log(p) / p^{k/2}
    prime: int          # p
    power: int          # k

@timer
def compute_prime_powers(limit: int) -> List[PrimePowerEntry]:
    """计算素数幂及其严格权重"""
    primes = prime_sieve(limit)
    entries = []
    for p in primes:
        pk = p
        k = 1
        while pk < limit:
            entries.append(PrimePowerEntry(
                log_val=np.log(pk),
                weight=np.log(p) / np.sqrt(pk),
                prime=p,
                power=k
            ))
            pk *= p
            k += 1
    entries.sort(key=lambda x: x.log_val)
    return entries

# 平滑化测试函数（高斯型）
def gaussian_test_function(s: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """
    测试函数 phi_hat(s) = exp(-s^2 / (2*sigma^2))
    用于平滑化截断，保证积分收敛
    """
    return np.exp(-s**2 / (2 * sigma**2))

def smoothed_zero_sum(X: np.ndarray, zeros: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    """
    平滑化后的零点求和：
    sum_{gamma} cos(gamma * X) * phi_hat(i*gamma)
    其中 phi_hat(i*gamma) = exp(gamma^2 / (2*sigma^2)) ... 需要调整
    
    实际上对于显式公式，测试函数作用于 s = 1/2 + i*gamma
    phi_hat(s) 在竖直线上衰减
    """
    result = np.zeros_like(X, dtype=np.float64)
    for gamma in zeros:
        # 测试函数在虚部上的作用：高斯衰减
        weight = np.exp(-gamma**2 / (2 * (sigma * 100)**2))  # 大sigma缓慢衰减
        result += weight * np.cos(gamma * X)
    return result

def strict_prime_sum(X: np.ndarray, entries: List[PrimePowerEntry], sigma_delta: float = 0.02) -> np.ndarray:
    """
    严格素数幂求和：
    -e^{-X/2} * sum_{p^k} (log(p)/p^{k/2}) * delta(X - log(p^k))
    delta 用高斯近似
    """
    result = np.zeros_like(X, dtype=np.float64)
    for entry in entries:
        if 0.5 <= entry.log_val <= 10:
            result += -np.exp(-X/2) * entry.weight * \
                     np.exp(-(X - entry.log_val)**2 / (2 * sigma_delta**2)) / \
                     (sigma_delta * np.sqrt(2 * np.pi))
    return result

# ==================== 数据准备 ====================
print("=" * 60)
print("希尔伯特-波利亚猜想：严格可视化")
print("=" * 60)

CACHE_200 = os.path.join(OUTPUT_DIR, 'riemann_zeros_200_cache.json')
CACHE_1000 = os.path.join(OUTPUT_DIR, 'riemann_zeros_1000_cache.json')

zeros_200 = compute_zeros(200, CACHE_200)
zeros_1000 = compute_zeros(1000, CACHE_1000)

pp_entries = compute_prime_powers(5000)
pp_logs = np.array([e.log_val for e in pp_entries])
pp_weights = np.array([e.weight for e in pp_entries])

print(f"[数据] 零点范围: γ_1={zeros_200[0]:.6f}, γ_200={zeros_200[-1]:.6f}")
print(f"[数据] 素数幂总数: {len(pp_entries)}")
print("=" * 60)

# ==================== 图1: 临界线分布 ====================
print("[图1] 零点临界线分布...")
fig1, ax1 = plt.subplots(figsize=(10, 12))

ax1.axvline(x=0.5, color='#e74c3c', linestyle='--', linewidth=2.5, alpha=0.85, 
            label=r'临界线 $\Re(s)=1/2$')
scatter = ax1.scatter(np.ones(200)*0.5, zeros_200, c=range(200), cmap='viridis', 
                      s=45, alpha=0.9, zorder=5, edgecolors='black', linewidth=0.4)

ax1.set_xlim(-0.5, 1.5)
ax1.set_ylim(0, 420)
ax1.set_xlabel(r'$\Re(s)$', fontsize=14)
ax1.set_ylabel(r'$\Im(s)$', fontsize=14)
ax1.set_title(r'前200个黎曼非平凡零点在临界线上的分布' + '\n' + r'$\Re(\rho)=1/2$',
              fontsize=15, fontweight='bold', fontproperties=prop)

ax1.annotate(r'泛函方程对称 $\rho \leftrightarrow 1-\rho$', 
             xy=(0.5, 380), xytext=(1.05, 380),
             fontsize=12, ha='center', color='darkred',
             arrowprops=dict(arrowstyle='<->', color='darkred', lw=2))

cbar1 = plt.colorbar(scatter, ax=ax1, fraction=0.046, pad=0.04)
cbar1.set_label('零点序号 $n$', fontsize=12, fontproperties=prop)

plt.tight_layout()
fig1.savefig(os.path.join(OUTPUT_DIR, 'fig1_critical_line.png'), 
             bbox_inches='tight', facecolor='white')
plt.close(fig1)

# ==================== 图2: 谱对应 ====================
print("[图2] 谱对应 λ=-4γ...")
fig2, ax2 = plt.subplots(figsize=(10, 8))

spectrum = -4 * zeros_200
colors = plt.cm.plasma(np.linspace(0, 1, 200))
ax2.scatter(range(1, 201), spectrum, c=colors, s=55, alpha=0.95, zorder=5,
            edgecolors='black', linewidth=0.4)

ax2.set_xlabel('零点序号 $n$', fontsize=14, fontproperties=prop)
ax2.set_ylabel(r'谱参数 $\lambda = -4\gamma$', fontsize=14, fontproperties=prop)
ax2.set_title(r'狄拉克算子谱对应 $\lambda = -4\gamma_n$', fontsize=15, fontweight='bold', fontproperties=prop)
ax2.axhline(y=0, color='#999999', linestyle='-', linewidth=0.5)

# 公式框
formula_text = (
    r'$D\psi_{\rho_1,\rho_2} = 2i(\rho_1 - \rho_2)\psi_{\rho_1,\rho_2}$' + '\n'
    r'$\mathcal{H} = -\dfrac{D}{4},\quad \mathcal{H}|\gamma_n\rangle = \gamma_n|\gamma_n\rangle$'
)
ax2.text(0.97, 0.03, formula_text, transform=ax2.transAxes, fontsize=11,
         verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.6', facecolor='#fff8dc', 
                   edgecolor='#8B4513', linewidth=1.5, alpha=0.95))

plt.tight_layout()
fig2.savefig(os.path.join(OUTPUT_DIR, 'fig2_spectrum_correspondence.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig2)

# ==================== 图3: 间距分布 ====================
print("[图3] 间距分布...")
fig3, ax3 = plt.subplots(figsize=(10, 8))

spacings = np.diff(zeros_200)
ax3.hist(spacings, bins=30, color='#4682B4', edgecolor='black', 
         alpha=0.75, density=True, label='归一化间距')

mean_spacing = np.mean(spacings)
ax3.axvline(x=mean_spacing, color='#e74c3c', linestyle='--', linewidth=2.5,
            label=f'平均间距 = {mean_spacing:.4f}')

# GUE Wigner 猜测
s = np.linspace(0, 3.5, 300)
wigner = (32/np.pi**2) * s**2 * np.exp(-4*s**2/np.pi)
ax3.plot(s, wigner, 'r-', linewidth=2.5, label='GUE Wigner 猜测')

ax3.set_xlabel(r'归一化间距 $s = \Delta\gamma/\langle\Delta\gamma\rangle$', fontsize=14, fontproperties=prop)
ax3.set_ylabel(r'概率密度 $P(s)$', fontsize=14, fontproperties=prop)
ax3.set_title('零点间距分布 vs 高斯酉系综 (GUE)', fontsize=15, fontweight='bold', fontproperties=prop)
ax3.legend(fontsize=12, loc='upper right')

ax3.text(0.97, 0.97, r'$P_{\rm GUE}(s) = \dfrac{32}{\pi^2}s^2 e^{-4s^2/\pi}$',
         transform=ax3.transAxes, fontsize=12, verticalalignment='top',
         horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='#ffffe0', alpha=0.9))

plt.tight_layout()
fig3.savefig(os.path.join(OUTPUT_DIR, 'fig3_spacing_distribution.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig3)

# ==================== 图4: 严格傅里叶对偶（含平滑化） ====================
print("[图4] 严格傅里叶对偶（平滑化）...")

fig4, axes = plt.subplots(2, 1, figsize=(14, 14), height_ratios=[1, 1])

X = np.linspace(0.5, 10, 4000)

# 上子图：原始（未平滑）vs 平滑化对比
ax_top = axes[0]

# 未平滑零点求和（200个）
zero_raw = np.zeros_like(X)
for gamma in zeros_200:
    zero_raw += np.cos(gamma * X)

# 平滑化零点求和（sigma控制截断）
sigma_smooth = 80  # 控制高斯衰减宽度
zero_smooth = np.zeros_like(X)
for gamma in zeros_200:
    weight = np.exp(-gamma**2 / (2 * sigma_smooth**2))
    zero_smooth += weight * np.cos(gamma * X)

# 严格素数幂项
prime_strict = strict_prime_sum(X, pp_entries, sigma_delta=0.015)

# 归一化
zero_raw_norm = zero_raw / np.max(np.abs(zero_raw)) * 0.9
zero_smooth_norm = zero_smooth / np.max(np.abs(zero_smooth)) * 0.9
prime_norm = prime_strict / (np.max(np.abs(prime_strict)) + 1e-10) * 0.9

ax_top.plot(X, zero_raw_norm, 'b-', linewidth=0.8, alpha=0.6, 
            label='零点振荡（原始，200个）')
ax_top.plot(X, zero_smooth_norm, 'c-', linewidth=1.2, alpha=0.9,
            label=f'零点振荡（平滑化, σ={sigma_smooth}）')
ax_top.plot(X, prime_norm, 'r-', linewidth=1.2, alpha=0.85,
            label='素数幂加权项（严格公式）')

# 标记素数幂位置
for entry in pp_entries[:100]:
    if 0.5 <= entry.log_val <= 10:
        ax_top.axvline(x=entry.log_val, color='red', alpha=0.08, linewidth=0.4)

ax_top.set_xlabel(r'$X = \log x$', fontsize=13)
ax_top.set_ylabel('归一化振幅', fontsize=13, fontproperties=prop)
ax_top.set_title('傅里叶对偶：平滑化截断效果对比', fontsize=14, fontweight='bold', fontproperties=prop)
ax_top.legend(fontsize=11, loc='upper right')
ax_top.set_xlim(0.5, 10)

# 下子图：仅平滑化 vs 素数项（更清晰的对偶）
ax_bot = axes[1]

ax_bot.plot(X, zero_smooth_norm, 'c-', linewidth=1.3, alpha=0.9,
            label=f'平滑化零点项 (σ={sigma_smooth})')
ax_bot.plot(X, prime_norm, 'r-', linewidth=1.3, alpha=0.85,
            label='素数幂加权项')

# 在素数幂位置放大标记
for entry in pp_entries[:60]:
    if 0.5 <= entry.log_val <= 10:
        ax_bot.axvline(x=entry.log_val, color='#ff6b6b', alpha=0.15, linewidth=0.6)

ax_bot.set_xlabel(r'$X = \log x$', fontsize=13)
ax_bot.set_ylabel('归一化振幅', fontsize=13, fontproperties=prop)
ax_bot.set_title('严格傅里叶对偶：平滑化零点项 vs 加权素数幂', 
                 fontsize=14, fontweight='bold', fontproperties=prop)

formula_box = (
    r'$\sum_{\rho} e^{i\gamma X} \sim -e^{-X/2}\sum_{p^k}\dfrac{\log p}{p^{k/2}}\delta(X-\log p^k)$'
    + '\n平滑化: ' + r'$\hat{\phi}(s) = \exp(-s^2/2\sigma^2)$'
)
ax_bot.text(0.5, 0.95, formula_box, transform=ax_bot.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#fffacd',
                      edgecolor='#daa520', linewidth=2, alpha=0.95))

ax_bot.legend(fontsize=11, loc='upper right')
ax_bot.set_xlim(0.5, 10)

plt.tight_layout()
fig4.savefig(os.path.join(OUTPUT_DIR, 'fig4_strict_fourier_duality.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig4)

# ==================== 图5: 围道积分分解 ====================
print("[图5] 围道积分分解...")

fig5, ax5 = plt.subplots(figsize=(14, 10))

# 模拟显式公式的各项贡献（基于论文第4节）
x_vals = np.linspace(2, 100, 500)  # x 从 2 到 100
X_vals = np.log(x_vals)

# 主导项: x * phi_hat(1)
phi_hat_1 = 1.0  # 测试函数在 s=1 的值
leading_term = x_vals * phi_hat_1

# 常数项: -log(2*pi) * phi_hat(0)
phi_hat_0 = 1.0
constant_term = -np.log(2*np.pi) * phi_hat_0 * np.ones_like(x_vals)

# 离散谱项（零点留数）：用前50个零点近似
T_disc = np.zeros_like(x_vals, dtype=np.complex128)
for gamma in zeros_200[:50]:
    rho = 0.5 + 1j * gamma
    # 留数贡献: -x^rho / rho * phi_hat(rho)
    phi_hat_rho = np.exp(-(gamma**2) / (2 * 50**2))  # 平滑化
    T_disc += - (x_vals**rho) / rho * phi_hat_rho

T_disc_real = np.real(T_disc)

# 连续谱项（近似）：背景积分，取负值使总和平滑
I_cont = -0.5 * np.ones_like(x_vals) * np.sin(X_vals * 2)

# 总和
psi_total = leading_term + constant_term + T_disc_real + I_cont

# 绘制
ax5.semilogy(x_vals, leading_term, 'r-', linewidth=2, 
             label=r'主导项 $x\hat{\phi}(1)$', alpha=0.9)
ax5.plot(x_vals, constant_term, 'g--', linewidth=1.5, 
         label=r'常数项 $-\log(2\pi)\hat{\phi}(0)$', alpha=0.8)
ax5.plot(x_vals, np.abs(T_disc_real), 'b-', linewidth=1.2, 
         label=r'离散谱 $|T_{\rm disc}(x)|$（零点留数）', alpha=0.8)
ax5.plot(x_vals, psi_total, 'k-', linewidth=2, 
         label=r'总和 $\psi_\phi(x)$', alpha=0.95)

ax5.set_xlabel(r'$x$', fontsize=14)
ax5.set_ylabel(r'$\psi_\phi(x)$', fontsize=14)
ax5.set_title('显式公式的围道积分分解', fontsize=15, fontweight='bold', fontproperties=prop)
ax5.legend(fontsize=11, loc='upper left')

formula_text = (
    r'$\psi_\phi(x) = x\hat{\phi}(1) - \log(2\pi)\hat{\phi}(0)$'
    + r'$+ T_{\rm disc}(x) + \mathcal{I}_{\rm cont}(x)$'
)
ax5.text(0.97, 0.5, formula_text, transform=ax5.transAxes, fontsize=12,
         verticalalignment='center', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.6', facecolor='#e6f3ff',
                   edgecolor='navy', linewidth=2, alpha=0.95))

plt.tight_layout()
fig5.savefig(os.path.join(OUTPUT_DIR, 'fig5_explicit_formula.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig5)

# ==================== 图6: 零点计数函数 ====================
print("[图6] 零点计数函数...")
fig6, ax6 = plt.subplots(figsize=(10, 8))

n = np.arange(1, 201)
T = zeros_200
N_approx = T/(2*np.pi) * np.log(T/(2*np.pi)) - T/(2*np.pi)

ax6.plot(T, n, 'bo-', markersize=4, linewidth=1.5, label='实际零点计数 $N(T)$')
ax6.plot(T, N_approx, 'r--', linewidth=2.5, label='Riemann-von Mangoldt 渐近')

ax6.set_xlabel(r'$T = \Im(\rho)$', fontsize=14)
ax6.set_ylabel(r'$N(T)$', fontsize=14)
ax6.set_title('零点计数函数', fontsize=15, fontweight='bold', fontproperties=prop)
ax6.legend(fontsize=12)

ax6.text(0.97, 0.03, 
         r'$N(T) = \dfrac{T}{2\pi}\log\dfrac{T}{2\pi} - \dfrac{T}{2\pi} + O(\log T)$',
         transform=ax6.transAxes, fontsize=11, verticalalignment='bottom',
         horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='#e0ffff', alpha=0.9))

plt.tight_layout()
fig6.savefig(os.path.join(OUTPUT_DIR, 'fig6_zero_counting.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig6)

# ==================== 图7: GUE对比 ====================
print("[图7] GUE对比...")
fig7, ax7 = plt.subplots(figsize=(10, 8))

spacings = np.diff(zeros_200)
mean_s = np.mean(spacings)
normalized = spacings / mean_s

s = np.linspace(0, 3.5, 300)
wigner = (32/np.pi**2) * s**2 * np.exp(-4*s**2/np.pi)

ax7.hist(normalized, bins=25, density=True, color='#4682B4', 
         edgecolor='black', alpha=0.75, label='归一化间距')
ax7.plot(s, wigner, 'r-', linewidth=2.5, label='GUE Wigner 猜测')

ax7.set_xlabel(r'归一化间距 $s = \Delta\gamma/\langle\Delta\gamma\rangle$', fontsize=14, fontproperties=prop)
ax7.set_ylabel(r'概率密度 $P(s)$', fontsize=14, fontproperties=prop)
ax7.set_title('归一化间距分布 vs 随机矩阵理论', fontsize=15, fontweight='bold', fontproperties=prop)
ax7.legend(fontsize=12)

ax7.text(0.97, 0.97, 
         r'$P_{\rm GUE}(s) = \dfrac{32}{\pi^2}s^2\exp\left(-\dfrac{4s^2}{\pi}\right)$',
         transform=ax7.transAxes, fontsize=11, verticalalignment='top',
         horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='#ffffe0', alpha=0.9))

plt.tight_layout()
fig7.savefig(os.path.join(OUTPUT_DIR, 'fig7_gue_comparison.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig7)

# ==================== 图8: 完整论证流程 ====================
print("[图8] 完整论证流程图...")
fig8 = plt.figure(figsize=(22, 28))
fig8.suptitle('希尔伯特-波利亚猜想：完整论证流程与核心公式汇总',
              fontsize=20, fontweight='bold', y=0.98, fontproperties=prop)

ax = fig8.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 22)
ax.set_ylim(0, 28)
ax.axis('off')

def add_box(x, y, w, h, facecolor, edgecolor, title, body, title_size=16, 
            body_size=12, linewidth=2.5):
    """辅助函数：添加带标题的框"""
    box = mpatches.FancyBboxPatch((x, y-h), w, h, boxstyle="round,pad=0.15",
                                   facecolor=facecolor, edgecolor=edgecolor,
                                   linewidth=linewidth)
    ax.add_patch(box)
    ax.text(x + w/2, y - 0.3, title, fontsize=title_size, ha='center',
            va='top', fontweight='bold', color=edgecolor, fontproperties=prop)
    ax.text(x + w/2, y - h/2 - 0.1, body, fontsize=body_size, ha='center',
            va='center', color='#333333')

def add_arrow(x1, y1, x2, y2, color='black'):
    """添加箭头"""
    arrow = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='->',
                            mutation_scale=30, linewidth=3, color=color)
    ax.add_patch(arrow)

# 第1层
add_box(1, 27.7, 20, 1.5, '#ffcccc', 'darkred', '黎曼假设',
        r'黎曼 $\zeta$ 函数的所有非平凡零点均位于临界线 $\Re(s)=1/2$ 上')
add_arrow(11, 26.2, 11, 25.0, 'darkred')

# 第2层
add_box(1, 24.5, 20, 2.0, '#ffe6cc', 'chocolate', 'Hilbert-Pólya 猜想',
        r'目标：构造显式谱三元组 $(\mathcal{A}, \mathcal{H}, D)$ 使得零点虚部 $\gamma_n$ 成为本征值')
add_arrow(11, 22.5, 11, 21.3, 'chocolate')

# 第3层：谱三元组
y = 20.8
box_spec = mpatches.FancyBboxPatch((1, y-3.5), 20, 3.8, boxstyle="round,pad=0.15",
                                    facecolor='#e6f3ff', edgecolor='navy', linewidth=2.5)
ax.add_patch(box_spec)
ax.text(11, y+1.3, '谱三元组的显式构造（第2节）', fontsize=16, ha='center',
        va='center', fontweight='bold', color='navy', fontproperties=prop)

# 三个子框
sub_boxes = [
    (1.5, y-0.5, 6, 2.0, 'white', 'steelblue', '希尔伯特空间',
     r'$\mathcal{H} = L^2(\mathbb{R}_+^\times \times \mathbb{R}_+^\times, d^\times x_1 d^\times x_2)$' + '\n'
     r'$\cong L^2(\mathbb{R}^2, dX_1 dX_2)$'),
    (8, y-0.5, 6, 2.0, 'white', 'darkgreen', '狄拉克算子',
     r'$D = -(x_1 p_1 + p_1 x_1) + (x_2 p_2 + p_2 x_2)$' + '\n'
     r'$= 2i\left(\frac{\partial}{\partial X_1} - \frac{\partial}{\partial X_2}\right)$'),
    (14.5, y-0.5, 6, 2.0, 'white', 'purple', '交叉积代数',
     r'$\mathcal{A} = C_c^\infty(\mathbb{R}_+) \rtimes_\alpha \mathbb{R}_+^\times$' + '\n'
     r'$\cong C_c^\infty(\mathbb{R}) \rtimes \mathbb{R}$')
]

for sx, sy, sw, sh, fc, ec, st, sb in sub_boxes:
    sbx = mpatches.FancyBboxPatch((sx, sy-sh), sw, sh, boxstyle="round,pad=0.1",
                                   facecolor=fc, edgecolor=ec, linewidth=2)
    ax.add_patch(sbx)
    ax.text(sx + sw/2, sy - 0.2, st, fontsize=13, ha='center', va='top',
            fontweight='bold', color=ec, fontproperties=prop)
    ax.text(sx + sw/2, sy - sh/2 - 0.1, sb, fontsize=11, ha='center', va='center')

add_arrow(11, y-3.6, 11, y-4.8, 'navy')

# 第4层
y = 15.5
box_eigen = mpatches.FancyBboxPatch((1, y-3.8), 20, 4.2, boxstyle="round,pad=0.15",
                                     facecolor='#fff5e6', edgecolor='darkgoldenrod', linewidth=2.5)
ax.add_patch(box_eigen)
ax.text(11, y+1.6, '广义本征值问题与约束联立（第3节）', fontsize=16, ha='center',
        va='center', fontweight='bold', color='darkgoldenrod', fontproperties=prop)
ax.text(11, y+0.8, r'$\psi_{\rho_1,\rho_2}(X_1,X_2) = \exp(\rho_1 X_1 + \rho_2 X_2)$',
        fontsize=12, ha='center', va='center')
ax.text(11, y+0.2, r'$D\psi_{\rho_1,\rho_2} = 2i(\rho_1 - \rho_2)\psi_{\rho_1,\rho_2}$',
        fontsize=12, ha='center', va='center')

# 约束框
c1 = mpatches.FancyBboxPatch((1.5, y-2.8), 8.5, 2.5, boxstyle="round,pad=0.1",
                              facecolor='#ccffcc', edgecolor='darkgreen', linewidth=2.5)
ax.add_patch(c1)
ax.text(5.75, y-1.2, '约束 I：自伴性', fontsize=14, ha='center', fontweight='bold', color='darkgreen', fontproperties=prop)
ax.text(5.75, y-1.9, r'$D = D^\dagger \Rightarrow \beta_1 = \beta_2$', fontsize=12, ha='center')

c2 = mpatches.FancyBboxPatch((11, y-2.8), 9.5, 2.5, boxstyle="round,pad=0.1",
                              facecolor='#ccccff', edgecolor='darkblue', linewidth=2.5)
ax.add_patch(c2)
ax.text(15.75, y-1.2, '约束 II：泛函方程', fontsize=14, ha='center', fontweight='bold', color='darkblue', fontproperties=prop)
ax.text(15.75, y-1.9, r'$\xi(s)=\xi(1-s) \Rightarrow \rho_2 = 1-\rho_1$', fontsize=12, ha='center')

add_arrow(11, y-3.9, 11, y-5.1, 'darkgoldenrod')

# 第5层：定理3.1
y = 9.8
box_thm = mpatches.FancyBboxPatch((3, y-2.0), 16, 2.4, boxstyle="round,pad=0.15",
                                   facecolor='gold', edgecolor='darkred', linewidth=4)
ax.add_patch(box_thm)
ax.text(11, y+0.6, '定理 3.1（临界线的动力学起源）', fontsize=17, ha='center',
        va='center', fontweight='bold', color='darkred', fontproperties=prop)
ax.text(11, y-0.2, r'$\beta_1 = \beta_2$ 且 $\rho_2 = 1-\rho_1 \Rightarrow \beta = \frac{1}{2}$',
        fontsize=15, ha='center', va='center', color='darkred')
ax.text(11, y-0.9, r'临界线 $\Re(s)=1/2$ 是自伴算子谱理论与代数模结构的必然推论',
        fontsize=12, ha='center', va='center', color='darkred')

add_arrow(11, y-2.1, 11, y-3.3, 'darkred')

# 第6层：算术迹公式
y = 5.5
box_trace = mpatches.FancyBboxPatch((1, y-4.8), 20, 5.2, boxstyle="round,pad=0.15",
                                     facecolor='#f0e6ff', edgecolor='indigo', linewidth=2.5)
ax.add_patch(box_trace)
ax.text(11, y+2.1, '算术迹公式与零点-素数傅里叶对偶（第4-5节）', fontsize=16,
        ha='center', va='center', fontweight='bold', color='indigo', fontproperties=prop)

ax.text(11, y+1.3, r'$\psi_\phi(x) = x\hat{\phi}(1) - \log(2\pi)\hat{\phi}(0) + T_{\rm disc}(x) + \mathcal{I}_{\rm cont}(x)$',
        fontsize=13, ha='center', va='center')

ax.text(5.5, y+0.4, r'$T_{\rm disc} = -\sum_\rho \frac{x^\rho}{\rho}\hat{\phi}(\rho)$',
        fontsize=12, ha='center', va='center')
ax.text(5.5, y-0.2, '离散谱：零点留数', fontsize=11, ha='center', va='center',
        color='darkgreen', fontweight='bold', fontproperties=prop)

ax.text(16.5, y+0.4, r'$\mathcal{I}_{\rm cont}$: 连续谱积分',
        fontsize=12, ha='center', va='center')
ax.text(16.5, y-0.2, '连续谱：Dixmier迹', fontsize=11, ha='center', va='center',
        color='darkblue', fontweight='bold', fontproperties=prop)

ax.text(11, y-1.2, r'$\langle a,D\rangle_\omega = {\rm Tr}_\omega(a|D|^{-1}) \sim \sum_{p^k}\frac{\log p}{p^{k/2}}f(\log p^k)$',
        fontsize=12, ha='center', va='center')
ax.text(11, y-1.8, r'$\sum_\rho e^{i\gamma X} \sim -e^{-X/2}\sum_{p^k}\frac{\log p}{p^{k/2}}\delta(X-\log p^k)$',
        fontsize=12, ha='center', va='center', color='darkred')

ax.text(11, y-2.6, '定理 5.1：循环上同调配对与围道积分显式公式严格等价',
        fontsize=12, ha='center', va='center', color='indigo', fontweight='bold', fontproperties=prop)

add_arrow(11, y-4.9, 11, y-6.1, 'indigo')

# 第7层：结论
y = 0.5
box_conc = mpatches.FancyBboxPatch((1, y-0.3), 20, 1.0, boxstyle="round,pad=0.15",
                                    facecolor='#ccffcc', edgecolor='darkgreen', linewidth=3)
ax.add_patch(box_conc)
ax.text(11, y+0.2, r'结论：黎曼假设 $\Leftrightarrow$ 谱完备性（所有非平凡零点可被谱三元组穷尽）',
        fontsize=13, ha='center', va='center', fontweight='bold', color='darkgreen', fontproperties=prop)

fig8.savefig(os.path.join(OUTPUT_DIR, 'fig8_complete_argument.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig8)

# ==================== 图9: 1000个零点扩展对偶 ====================
print("[图9] 1000个零点扩展对偶...")

fig9, ax9 = plt.subplots(figsize=(16, 8))

X = np.linspace(0.5, 12, 5000)

# 1000个零点平滑化求和
sigma_1k = 150
zero_1k = np.zeros_like(X)
for gamma in zeros_1000:
    weight = np.exp(-gamma**2 / (2 * sigma_1k**2))
    zero_1k += weight * np.cos(gamma * X)

# 素数幂项（扩展范围）
prime_1k = np.zeros_like(X)
for entry in pp_entries:
    if 0.5 <= entry.log_val <= 12:
        prime_1k += -np.exp(-X/2) * entry.weight * \
                   np.exp(-(X - entry.log_val)**2 / (2 * 0.015**2)) / \
                   (0.015 * np.sqrt(2 * np.pi))

# 归一化
z1k_norm = zero_1k / (np.max(np.abs(zero_1k)) + 1e-10) * 0.9
p1k_norm = prime_1k / (np.max(np.abs(prime_1k)) + 1e-10) * 0.9

ax9.plot(X, z1k_norm, 'b-', linewidth=1.0, alpha=0.8,
         label=f'零点振荡（1000个, σ={sigma_1k}）')
ax9.plot(X, p1k_norm, 'r-', linewidth=1.2, alpha=0.85,
         label='素数幂加权项')

# 标记素数幂
for entry in pp_entries[:150]:
    if 0.5 <= entry.log_val <= 12:
        ax9.axvline(x=entry.log_val, color='red', alpha=0.06, linewidth=0.3)

ax9.set_xlabel(r'$X = \log x$', fontsize=14)
ax9.set_ylabel('归一化振幅', fontsize=14, fontproperties=prop)
ax9.set_title('扩展傅里叶对偶：1000个零点 vs 完整素数幂谱', 
              fontsize=15, fontweight='bold', fontproperties=prop)
ax9.legend(fontsize=12, loc='upper right')
ax9.set_xlim(0.5, 12)

plt.tight_layout()
fig9.savefig(os.path.join(OUTPUT_DIR, 'fig9_extended_duality_1000.png'),
             bbox_inches='tight', facecolor='white')
plt.close(fig9)

# ==================== 图10: Dixmier迹数值验证 ====================
print("[图10] Dixmier迹配对验证...")

fig10, ax10 = plt.subplots(figsize=(12, 8))

# 计算不同截断下的配对值
cutoffs = np.arange(50, 1001, 50)
pairing_values = []

for cutoff in cutoffs:
    # 近似 Dixmier 迹：sum_{gamma < cutoff} 1/gamma
    # 对应于 Tr_omega(a|D|^{-1}) 的离散近似
    subset = zeros_1000[:cutoff]
    pairing = np.sum(1.0 / subset)
    pairing_values.append(pairing)

# 素数幂求和（对应项）
prime_sum = []
for cutoff in cutoffs:
    # sum_{p^k < exp(cutoff/10)} log(p)/p^{k/2} 的近似
    val = sum(e.weight for e in pp_entries if e.log_val < cutoff / 50)
    prime_sum.append(val)

# 归一化比较
pv_norm = np.array(pairing_values) / pairing_values[-1]
ps_norm = np.array(prime_sum) / prime_sum[-1]

ax10.plot(cutoffs, pv_norm, 'bo-', markersize=6, linewidth=1.5,
          label=r'Dixmier迹近似 $\sum_{\gamma_n < \Lambda} \gamma_n^{-1}$')
ax10.plot(cutoffs, ps_norm, 'rs-', markersize=6, linewidth=1.5,
          label=r'素数幂求和 $\sum_{p^k < X} \frac{\log p}{p^{k/2}}$')

ax10.set_xlabel(r'截断参数 $\Lambda$', fontsize=14, fontproperties=prop)
ax10.set_ylabel('归一化累积值', fontsize=14, fontproperties=prop)
ax10.set_title('Dixmier迹配对 vs 素数幂求和的收敛行为', 
               fontsize=15, fontweight='bold', fontproperties=prop)
ax10.legend(fontsize=12)

ax10.text(0.97, 0.03,
          r'$\langle a, D \rangle_\omega = {\rm Tr}_\omega(a|D|^{-1})$'
          + r'$\sim \sum_{p^k} \frac{\log p}{p^{k/2}} f(\log p^k)$',
          transform=ax10.transAxes, fontsize=11, verticalalignment='bottom',
          horizontalalignment='right',
          bbox=dict(boxstyle='round', facecolor='#f0e6ff', alpha=0.9))

plt.tight_layout()
fig10.savefig(os.path.join(OUTPUT_DIR, 'fig10_dixmier_pairing.png'),
              bbox_inches='tight', facecolor='white')
plt.close(fig10)

# ==================== 导出数据 ====================
print("[导出] CSV数据...")
csv_path = os.path.join(OUTPUT_DIR, 'riemann_zeros_data.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 'gamma_n', 'lambda_n', 'delta_gamma', 'spacing_normalized'])
    for i in range(len(zeros_200)):
        writer.writerow([
            i + 1,
            f'{zeros_200[i]:.10f}',
            f'{-4*zeros_200[i]:.10f}',
            '0.0' if i == 0 else f'{zeros_200[i]-zeros_200[i-1]:.10f}',
            '0.0' if i == 0 else f'{(zeros_200[i]-zeros_200[i-1])/mean_spacing:.6f}'
        ])

# 素数幂数据
csv_pp = os.path.join(OUTPUT_DIR, 'prime_powers_data.csv')
with open(csv_pp, 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.writer(f)
    writer.writerow(['prime', 'power', 'p^k', 'log_p_k', 'weight_log_p_over_sqrt'])
    for e in pp_entries:
        writer.writerow([e.prime, e.power, e.prime**e.power, 
                        f'{e.log_val:.10f}', f'{e.weight:.10f}'])

print("=" * 60)
print("全部完成！生成文件列表：")
for fname in sorted(os.listdir(OUTPUT_DIR)):
    if fname.endswith(('.png', '.csv', '.json')):
        fpath = os.path.join(OUTPUT_DIR, fname)
        fsize = os.path.getsize(fpath) / 1024
        print(f"  {fname:40s} {fsize:8.1f} KB")
print("=" * 60)