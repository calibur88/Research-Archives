import matplotlib.pyplot as plt
import numpy as np
import os

# 创建保存目录
save_dir = "/storage/emulated/0/研究报告/黎曼猜想物理版/"
os.makedirs(save_dir, exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), dpi=150)

# ========== 左图：零点轨迹相图 ==========
ax1.set_xlim(0, 60)
ax1.set_ylim(0, 9)

# 生成前几个黎曼零点的虚部
zeros_t = [14.13, 21.02, 25.01, 30.42, 32.93, 37.58, 40.91, 43.32, 48.00, 52.97, 56.45]

# 绘制水平轨迹线（表示不同零点随时间的演化）
for i, t in enumerate(zeros_t):
    # 实际零点（黄色星星，位于σ=0.5）
    ax1.scatter(t, 0.5, marker='*', s=300, color='yellow', 
                edgecolors='black', linewidths=1, zorder=5)
    
    # 绘制水平细线表示轨迹收敛到0.5
    sigma_line = np.linspace(0.1, 0.9, 50)
    # 添加微小波动模拟收敛过程
    noise = 0.02 * np.sin(t * 0.5) * np.exp(-abs(sigma_line - 0.5)*10)
    ax1.plot([t]*50, sigma_line + noise, color='lightblue', 
             alpha=0.3, linewidth=0.8)

# 添加"瞬时锁定"零点（青色星星，位于更高σ值，表示初始状态）
instant_zeros_t = [15, 22, 26, 31]
for t in instant_zeros_t:
    ax1.scatter(t, 0.5, marker='*', s=200, color='cyan', 
                edgecolors='black', linewidths=1, alpha=0.8, zorder=5)

# 临界线（绿色虚线）
ax1.axhline(0.5, color='green', linestyle='--', linewidth=2.5, 
            label='Critical Line σ=0.5 (Attractor)')

ax1.set_xlabel('Imaginary Part t', fontsize=12, fontweight='bold')
ax1.set_ylabel('Real Part σ', fontsize=12, fontweight='bold')
ax1.set_title('Phase Space of Abyss Flow\nZeros Converging to Attractor', 
              fontsize=13, fontweight='bold')
ax1.legend(loc='upper left')
ax1.grid(True, alpha=0.3)

# ========== 右图：吸引子势阱热力图 ==========
ax2.set_xlim(10, 65)
ax2.set_ylim(0.1, 0.9)

# 创建网格
t_grid = np.linspace(10, 65, 200)
sigma_grid = np.linspace(0.1, 0.9, 100)
T, Sigma = np.meshgrid(t_grid, sigma_grid)

# 计算"深渊流"势能：在σ=0.5处有最小值（吸引子）
# 使用高斯势阱 + 随t的微小调制
potential = ((Sigma - 0.5)**2) * 100 + 0.1 * np.sin(T * 0.2)

# 反转颜色：黑色=低能量（吸引子），黄色=高能量
im = ax2.imshow(potential, aspect='auto', cmap='hot_r',  # _r表示反转
                extent=[10, 65, 0.1, 0.9], origin='lower',
                vmin=0, vmax=5)

# 标记零点位置（青色星星）
for t in zeros_t[:8]:
    ax2.scatter(t, 0.5, marker='*', s=200, color='cyan', 
                edgecolors='white', linewidths=1, zorder=5, label='Zeros (Instant Lock)')

ax2.axhline(0.5, color='white', linestyle='--', linewidth=2, alpha=0.8)

ax2.set_xlabel('Imaginary Part t', fontsize=12, fontweight='bold')
ax2.set_ylabel('Real Part σ', fontsize=12, fontweight='bold')
ax2.set_title('Attracting Potential Landscape\n(Black=Minimum Energy)', 
              fontsize=13, fontweight='bold')
ax2.legend(loc='upper right')

# 添加颜色条
cbar = plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
cbar.set_label('Potential Energy |ψ|', rotation=270, labelpad=20)

plt.tight_layout()
save_path = os.path.join(save_dir, "Abyss_Flow_Attraction_Space.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"双联图已保存至: {save_path}")
plt.show()
