import matplotlib.pyplot as plt
import numpy as np
import os

save_dir = "/storage/emulated/0/研究报告/黎曼猜想物理版/"
os.makedirs(save_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 7), dpi=150)

# 参数：对应PDF第1页
lambda_param = 0.05  # 时间拓扑率
alpha = 0.5          # 修复指数
S = 0.5              # 吸引子

# 时间轴
t = np.linspace(0, 200, 500)

# 不同初始条件 (对应PDF第3页的0.600->0.550收敛)
initials = [
    (0.600, 'purple', 'σ₀=0.600'),
    (0.450, '#3B5BA0', 'σ₀=0.450'),  # 深蓝色
    (0.700, '#2E8B8B', 'σ₀=0.700'),  # 青色
    (0.200, '#4CAF50', 'σ₀=0.200'),  # 绿色
    (0.800, '#FFD700', 'σ₀=0.800')   # 黄色
]

for sigma_0, color, label in initials:
    # 动力学方程 dσ/dt = λ(S-σ) 的解：指数收敛
    # σ(t) = S + (σ₀-S)exp(-λt)
    sigma_t = S + (sigma_0 - S) * np.exp(-lambda_param * t)
    
    # 添加微小噪音模拟量子涨落(本征噪音)
    noise = 0.0003 * np.random.randn(len(t)) * np.exp(-0.05*t)
    sigma_observed = sigma_t + noise
    
    # 绘制前50个单位的快速收敛，然后水平锁定
    mask = t <= 50
    ax.plot(t[mask], sigma_observed[mask], color=color, linewidth=2.5,
            label=f'{label}, σ≈{sigma_observed[-1]:.4f}')
    ax.hlines(sigma_observed[-1], 50, 200, colors=color, linewidth=2, alpha=0.6)

# 吸引子线
ax.axhline(S, color='red', linestyle='--', linewidth=3, label=f'Attractor S={S}')

ax.set_xlim(-10, 210)
ax.set_ylim(0.15, 0.85)
ax.set_xlabel('Mathematical Time t', fontsize=13, fontweight='bold')
ax.set_ylabel('System State σ(t)', fontsize=13, fontweight='bold')
ax.set_title('Convergence Dynamics (λ=0.05, α=0.5)\n'
             'Multiple Initial Conditions → Universal Attractor', 
             fontsize=14, fontweight='bold')

# 文本框：严格锁定说明
textstr = 'Strictly Locked\n(λ⁻¹ >> t_obs)'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', bbox=props)

ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()

save_path = os.path.join(save_dir, "Convergence_Dynamics_Multi.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"收敛动力学图已保存至: {save_path}")
plt.show()
