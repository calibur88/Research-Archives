
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# 前10个黎曼零点虚部
gammas = np.array([
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832
])

# 放大10倍
scale = 10
center_x = 5  # 0.5 * 10
scaled_gammas = scale * gammas

fig, ax = plt.subplots(figsize=(16, 10))

# 绘制坐标轴
ax.axhline(y=0, color='black', linewidth=2, zorder=1)
ax.axvline(x=0, color='black', linewidth=2, zorder=1)

# 绘制新的临界线（原0.5变为5）
ax.axvline(x=center_x, color='red', linewidth=3, linestyle='--', 
          label=f'Critical Line (scaled): Re(s) = {center_x}', zorder=2)

# 绘制放大后的临界带边界（原0和1变为0和10）
ax.axvline(x=0, color='blue', linewidth=1.5, linestyle=':', alpha=0.7)
ax.axvline(x=10, color='blue', linewidth=1.5, linestyle=':', alpha=0.7)
ax.fill_betweenx([-max(scaled_gammas)*1.2, max(scaled_gammas)*1.2], 0, 10, 
                alpha=0.1, color='lightblue', label='Scaled Critical Strip (0,10)')

# 为每个零点绘制同心圆和交点
colors = plt.cm.viridis(np.linspace(0, 1, len(gammas)))

for i, (gamma_orig, gamma_scaled) in enumerate(zip(gammas, scaled_gammas)):
    # 绘制圆 (圆心(5,0), 半径10γ)
    circle = Circle((center_x, 0), gamma_scaled, fill=False, 
                   color=colors[i], linewidth=2, alpha=0.8)
    ax.add_patch(circle)
    
    # 计算与实轴交点
    left_x = center_x - gamma_scaled   # 5 - 10γ
    right_x = center_x + gamma_scaled  # 5 + 10γ
    
    # 绘制交点（仅每几个标注避免拥挤）
    if i < 5:  # 只标注前5个避免图形太乱
        ax.plot(left_x, 0, 'o', color=colors[i], markersize=8, markeredgecolor='black', markeredgewidth=1)
        ax.plot(right_x, 0, 's', color=colors[i], markersize=8, markeredgecolor='black', markeredgewidth=1)
        
        # 标注数值
        ax.annotate(f'{left_x:.1f}', xy=(left_x, 0), xytext=(left_x-15, -8-2*i),
                   fontsize=9, color=colors[i], fontweight='bold')
        ax.annotate(f'{right_x:.1f}', xy=(right_x, 0), xytext=(right_x+2, -8-2*i),
                   fontsize=9, color=colors[i], fontweight='bold')
        
        # 绘制半径线到上半零点位置
        ax.plot([center_x, center_x], [0, gamma_scaled], color=colors[i], 
               linewidth=1.5, linestyle='-', alpha=0.6)
        ax.plot(center_x, gamma_scaled, '^', color=colors[i], markersize=7, 
               markeredgecolor='black', markeredgewidth=0.5)

# 绘制共同圆心
ax.plot(center_x, 0, 'r*', markersize=20, markerfacecolor='red', 
       markeredgecolor='darkred', markeredgewidth=2, label='Common Center (5,0)', zorder=10)

# 设置坐标轴
max_gamma = max(scaled_gammas)
ax.set_xlim(center_x - max_gamma - 20, center_x + max_gamma + 20)
ax.set_ylim(-max_gamma - 10, max_gamma + 10)

ax.set_xlabel('Re(s) - Real Part (×10 scaled)', fontsize=13, fontweight='bold')
ax.set_ylabel('Im(s) - Imaginary Part (×10 scaled)', fontsize=13, fontweight='bold')
ax.set_title('Concentric Circles from Scaled Zeros: 10s = (10-10s)\n' + 
            f'Center: ({center_x}, 0), Radius = 10×γ for each zero', 
            fontsize=14, fontweight='bold')

# 添加图例和信息
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)

# 添加说明文本
info_text = (
    "Scaling: s → 10s\n"
    f"Original center: (0.5, 0) → Scaled center: ({center_x}, 0)\n"
    f"Radius for γ={gammas[0]:.2f}: {scaled_gammas[0]:.2f}\n"
    f"Intersections: {center_x} ± 10γ\n"
    f"Leftmost point: {center_x - max(scaled_gammas):.1f}\n"
    f"Rightmost point: {center_x + max(scaled_gammas):.1f}"
)
ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('/storage/emulated/0/研究报告/黎曼猜想/riemann_concentric_circles_scaled.png', dpi=150, bbox_inches='tight')
plt.show()

# 输出数值表格
print("=" * 80)
print("同心圆分析：s → 10s 缩放变换 (10s = 10 - 10s)")
print("=" * 80)
print(f"{'Zero#':<6} {'Original γ':<15} {'Scaled Radius(10γ)':<20} {'Left(5-10γ)':<18} {'Right(5+10γ)':<18}")
print("-" * 80)

for i, gamma in enumerate(gammas, 1):
    radius = 10 * gamma
    left = 5 - radius
    right = 5 + radius
    print(f"{i:<6} {gamma:<15.4f} {radius:<20.4f} {left:<18.4f} {right:<18.4f}")

print("-" * 80)
print(f"共同圆心：({center_x}, 0)  [原 (0.5, 0) 放大10倍]")
print(f"左侧交点规律：5 - 10γ（均为负值，随着γ增大向左远离）")
print(f"右侧交点规律：5 + 10γ（均为正值，范围从 {5+10*gammas[0]:.1f} 到 {5+10*gammas[-1]:.1f}）")
print("=" * 80)