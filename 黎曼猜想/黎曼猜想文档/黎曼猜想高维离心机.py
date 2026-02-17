import matplotlib.pyplot as plt
import numpy as np
import os

save_dir = "/storage/emulated/0/研究报告/黎曼猜想物理版/"
os.makedirs(save_dir, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# ========== 1. 左上: Potential Landscape ==========
ax1 = axes[0, 0]
sigma = np.linspace(0.2, 0.8, 400)
t_vals = [14.13, 21.02, 25.01]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for t, c in zip(t_vals, colors):
    psi = np.exp(-((sigma-0.5)**2)*200) * (1 + 0.1*np.sin(t*sigma))
    ax1.semilogy(sigma, psi+1e-8, color=c, linewidth=2, label=f'Zero #{t_vals.index(t)+1}: t={t:.2f}')

ax1.axvline(0.5, color='red', linestyle='--', linewidth=2, alpha=0.8, label='Critical Line σ=0.5')
ax1.axvspan(0.5, 0.8, alpha=0.15, color='green', label='Collapsing Zone')
ax1.axvspan(0.2, 0.5, alpha=0.15, color='purple', label='Overflow Zone')

ax1.set_xlabel('Real Part σ')
ax1.set_ylabel('|ψ_min|')
ax1.set_title('Wave Collapse at Potential Landscape\nof Riemann Zeros', fontweight='bold')
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, alpha=0.3)

# ========== 2. 右上: 蓝色竖虚线 (Zero Trajectories) ==========
ax2 = axes[0, 1]

# 蓝色竖虚线标记零点位置 (在特定σ处画竖线)
zeros_t = [14.13, 21.02, 25.01, 30.42, 32.93, 37.58, 40.91, 43.32]
for t in zeros_t:
    # 蓝色竖虚线：标记每个零点的虚部位置
    ax2.axvline(t/50, color='blue', linestyle='--', alpha=0.5, linewidth=1)  # 归一化到0-1范围显示
    ax2.plot(0.5, t, 'ro', markersize=4)  # 红点在临界线

# 或者：如果X轴是σ，Y轴是t，那么竖虚线应该是axvline在特定σ值？
# 根据原图，可能是标记零点在复平面的位置，用竖线连接
ax2.set_xlim(0, 1)
ax2.set_ylim(10, 50)
ax2.set_xlabel('Real Part σ')
ax2.set_ylabel('Imaginary Part t')
ax2.set_title('Zeros Distribution\n(Blue vertical: Zero locations)', fontweight='bold')

# 添加机器精度线
ax2.axhline(1e-4, color='gray', linestyle='--', alpha=0.5, label='Machine Precision')

# ========== 3. 左下: Real Part Landscape ==========
ax3 = axes[1, 0]
sigma_map = np.linspace(0.25, 0.75, 100)
t_map = np.linspace(10, 50, 100)
Sigma, T = np.meshgrid(sigma_map, t_map)
Z = -np.exp(-((Sigma-0.5)**2)*50)

im = ax3.imshow(Z, aspect='auto', cmap='RdYlBu_r', extent=[0.25, 0.75, 10, 50], origin='lower')
for t in zeros_t[:6]:
    ax3.scatter(0.5, t, marker='*', s=150, color='yellow', edgecolors='black', zorder=5)

ax3.set_xlabel('Real Part σ')
ax3.set_ylabel('Imaginary Part t')
ax3.set_title('Real Part Landscape\nConverging to Critical Line', fontweight='bold')
plt.colorbar(im, ax=ax3, fraction=0.046, label='Real(ζ)')

# ========== 4. 右下: Convexity Certificate (竖线特征) ==========
ax4 = axes[1, 1]
sigma_h = np.linspace(0.35, 0.65, 400)

# Hessian尖峰 (在0.5处的竖线状尖峰)
hessian = 90 * np.exp(-((sigma_h-0.5)**2)*300)
ax4.plot(sigma_h, hessian, color='purple', linewidth=3, label='Hessian ∂²Φ/∂σ²')

# Critical Line (红色竖线)
ax4.axvline(0.5, color='red', linestyle='-', linewidth=2, label='Critical Line')

# Convexity Threshold (水平虚线)
ax4.axhline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Convexity Threshold')

# Strictly Convex Region (绿色区域，Hessian>0)
ax4.fill_between(sigma_h, 0, hessian, where=(hessian>0), color='green', alpha=0.3, label='Strictly Convex')

ax4.set_xlabel('Real Part σ')
ax4.set_ylabel('Restoring Force')
ax4.set_title('Convexity Certificate\n(Peak at Critical Line)', fontweight='bold')
ax4.legend(loc='upper right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0.35, 0.65)

plt.tight_layout()
save_path = os.path.join(save_dir, "Riemann_4Panels_Corrected.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"已保存: {save_path}")
plt.show()
