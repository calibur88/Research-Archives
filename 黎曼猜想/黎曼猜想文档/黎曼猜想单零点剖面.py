import matplotlib.pyplot as plt
import numpy as np
import os

save_dir = "/storage/emulated/0/研究报告/黎曼猜想物理版/"
os.makedirs(save_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 7), dpi=150)

# 数据：第一个零点的虚部
t_zero = 14.1347

# X轴：实部 σ (0到1)
sigma = np.linspace(0.1, 0.9, 500)

# 应力场：在σ=0.5处为尖峰（三角形/帐篷状），向两侧线性衰减
# 这是"修复力" R(σ-S) 的绝对值 |σ - 0.5| 的倒数或线性锥形
peak_height = 90  # 峰值高度

# 创建三角形尖峰
stress_field = peak_height * (1 - 2 * np.abs(sigma - 0.5))
stress_field = np.maximum(stress_field, 0)  # 底部截断到0

# 绘制应力场（紫色粗线）
ax.plot(sigma, stress_field, color='purple', linewidth=4, 
        label=f'Stress Field @ t={t_zero:.4f}')

# 临界线（红色虚线）
ax.axvline(0.5, color='#FF6B6B', linestyle='--', linewidth=3, 
           label='Critical Line σ=0.5')

# 区域标记
ax.axvspan(0.5, 0.6, alpha=0.2, color='green', label='Collapsing Zone')
ax.axvspan(0.4, 0.5, alpha=0.2, color='blue', label='Overflow Zone')

# 设置坐标轴
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-5, 100)
ax.set_xlabel('Real Part σ', fontsize=14, fontweight='bold')
ax.set_ylabel('Stress Intensity |R(σ-S)|', fontsize=14, fontweight='bold')
ax.set_title('One-Dimensional Wave Collapse at the Critical Line\n'
             'Abyss Stress Field @ σ=0.5 (The "Attractor" Mechanism)', 
             fontsize=13, fontweight='bold')

# 图例
ax.legend(loc='lower center', fontsize=11, framealpha=0.9)

# 添加网格
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# 添加注释
ax.annotate('Peak Stress\n(Maximum Restoring Force)', 
            xy=(0.5, 90), xytext=(0.65, 75),
            arrowprops=dict(arrowstyle='->', color='red', lw=2),
            fontsize=11, color='red', fontweight='bold')

ax.annotate('Linear Decay\n(Collapsing)', 
            xy=(0.6, 45), xytext=(0.75, 30),
            arrowprops=dict(arrowstyle='->', color='green', lw=1.5),
            fontsize=10, color='green')

ax.annotate('Linear Decay\n(Overflow)', 
            xy=(0.4, 45), xytext=(0.15, 30),
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.5),
            fontsize=10, color='blue')

plt.tight_layout()
save_path = os.path.join(save_dir, "Stress_Field_Profile_t14.1347.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"应力场剖面图已保存至: {save_path}")
plt.show()
