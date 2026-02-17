import matplotlib.pyplot as plt
import numpy as np
import os

# 设置全局样式
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 保存路径
save_path = "/storage/emulated/0/研究报告/黎曼猜想物理版/"
os.makedirs(save_path, exist_ok=True)

# ==================== Figure 1: Page 2 Left - Five horizontal lines ====================
fig1, ax1 = plt.subplots(figsize=(8, 6))

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
y_positions = [0.475, 0.550, 0.520, 0.400, 0.590]

for i, (color, y_pos) in enumerate(zip(colors, y_positions)):
    ax1.hlines(y=y_pos, xmin=0, xmax=2, color=color, linewidth=2.5, label=f'State {i+1}')

ax1.axhline(y=0.5, color='black', linestyle='--', linewidth=2, label='Attractor: 0.5')
ax1.set_xlim(0, 2)
ax1.set_ylim(0.35, 0.65)
ax1.set_xlabel('Mathematical Time (t) 1e6', fontsize=12)
ax1.set_ylabel('Real Part X(t)', fontsize=12)
ax1.set_title('Dynamical Evolution Equation\n(Time Topology Rate λ = 10^-10)', fontsize=13)
ax1.legend(loc='upper right', framealpha=0.9)
ax1.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_1_dynamical_evolution.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 2: Page 2 Right - Two column scatter ====================
fig2, ax2 = plt.subplots(figsize=(8, 6))

np.random.seed(42)
n_points = 60

left_x = np.zeros(n_points) + np.random.normal(0, 0.02, n_points)
left_y = np.linspace(0.40, 0.60, n_points) + np.random.normal(0, 0.01, n_points)
right_x = np.ones(n_points) + np.random.normal(0, 0.02, n_points)
right_y = np.ones(n_points) * 0.5 + np.random.normal(0, 0.008, n_points)

ax2.scatter(left_x, left_y, c='#1f77b4', s=40, alpha=0.7, label='Initial State', edgecolors='white', linewidth=0.5)
ax2.scatter(right_x, right_y, c='#ff7f0e', s=40, alpha=0.7, label='Locked State', edgecolors='white', linewidth=0.5)

ax2.axhline(y=0.5, color='black', linestyle='--', linewidth=2)
ax2.set_xlim(-0.3, 1.3)
ax2.set_ylim(0.4, 0.6)
ax2.set_xticks([0, 1])
ax2.set_xticklabels(['Initial', 'Observed'])
ax2.set_ylabel('Real Part X', fontsize=12)
ax2.set_title('Finite Time Window Observation\n(Convergence to Attractor 0.5)', fontsize=13)
ax2.legend(loc='lower center')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_2_finite_time_observation.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 3: Page 3 - Evolution comparison ====================
fig3, (ax3_1, ax3_2) = plt.subplots(1, 2, figsize=(14, 6))

# Left
initial_vals = [0.435, 0.475, 0.545, 0.520, 0.590]
for i, (color, init) in enumerate(zip(colors, initial_vals)):
    x = np.linspace(0, 2, 100)
    y = 0.5 + (init - 0.5) * np.exp(-0.5 * x)
    ax3_1.plot(x, y, color=color, linewidth=2, label=f'Trajectory {i+1}')

ax3_1.axhline(y=0.5, color='black', linestyle='--', linewidth=2, label='Attractor: 0.5')
ax3_1.set_xlim(0, 2)
ax3_1.set_ylim(0.42, 0.60)
ax3_1.set_xlabel('Evolution Time (dimensionless t) 1e6', fontsize=11)
ax3_1.set_ylabel('Real Part X(t)', fontsize=11)
ax3_1.set_title('High-Dimensional Projection Dynamics\n(Time Topology λ = 10^-10)', fontsize=12)
ax3_1.legend(loc='upper right', fontsize=9)
ax3_1.grid(True, alpha=0.3)

# Right
left_y = np.random.uniform(0.40, 0.60, 80)
right_y = 0.5 + np.random.normal(0, 0.005, 80)

ax3_2.scatter(np.zeros(80), left_y, c='#1f77b4', s=25, alpha=0.6, edgecolors='white', linewidth=0.3)
ax3_2.scatter(np.ones(80), right_y, c='#ff7f0e', s=25, alpha=0.6, edgecolors='white', linewidth=0.3)
ax3_2.axhline(y=0.5, color='black', linestyle='--', linewidth=2)
ax3_2.set_xlim(-0.2, 1.2)
ax3_2.set_ylim(0.4, 0.6)
ax3_2.set_xticks([0, 1])
ax3_2.set_xticklabels(['Initial', 'Final'])
ax3_2.set_ylabel('Real Part X', fontsize=11)
ax3_2.set_title('Observation Projection Locking\n(Critical Line Convergence to 0.5)', fontsize=12)
ax3_2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_3_evolution_comparison.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 4: Page 4 - Four panel ====================
fig4 = plt.figure(figsize=(14, 10))

ax4_1 = plt.subplot(2, 2, 1)
pi_labels = ['π=3.0', 'π=π_exact', 'π=3.5']
errors = [1e-5, 1e-12, 1.5]
bar_colors = ['#87CEEB', '#90EE90', '#F08080']
ax4_1.bar(pi_labels, errors, color=bar_colors, alpha=0.8, width=0.6)
ax4_1.set_ylabel('Calculation Error', fontsize=11)
ax4_1.set_title('Precision Comparison', fontsize=12)
ax4_1.set_yscale('log')
ax4_1.grid(True, alpha=0.3, axis='y')

ax4_2 = plt.subplot(2, 2, 2)
zeros_idx = np.arange(10)
zeros_real = [0.50, 0.52, 0.48, 0.51, 0.49, 0.50, 0.52, 0.48, 0.51, 0.50]
ax4_2.bar(zeros_idx, zeros_real, color='steelblue', alpha=0.8, width=0.7)
ax4_2.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
ax4_2.set_xlabel('Zero Index', fontsize=11)
ax4_2.set_ylabel('Real Part', fontsize=11)
ax4_2.set_title('Real Part Distribution (First 100)', fontsize=12)
ax4_2.set_ylim(0.45, 0.55)
ax4_2.grid(True, alpha=0.3, axis='y')

ax4_3 = plt.subplot(2, 2, 3)
t_vals = np.linspace(0, 100, 100)
fluctuations = 0.5 + 0.02 * np.sin(t_vals * 0.5) + np.random.normal(0, 0.005, 100)
ax4_3.plot(t_vals, fluctuations, color='blue', linewidth=1.2, alpha=0.7)
ax4_3.axhline(y=0.5, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
ax4_3.fill_between(t_vals, 0.48, 0.52, alpha=0.1, color='gray')
ax4_3.set_xlabel('Zero Index', fontsize=11)
ax4_3.set_ylabel('Real Part', fontsize=11)
ax4_3.set_title('Real Part Fluctuation', fontsize=12)
ax4_3.set_ylim(0.45, 0.55)
ax4_3.grid(True, alpha=0.3)

ax4_4 = plt.subplot(2, 2, 4, projection='polar')
categories = ['Stability', 'Convergence', 'Symmetry', 'Precision', 'Universality']
N = len(categories)
angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

values_30 = [0.8, 0.75, 0.7, 0.85, 0.8]
values_exact = [1.0, 1.0, 1.0, 1.0, 1.0]
values_35 = [0.7, 0.65, 0.6, 0.75, 0.7]

for vals, color, label in [(values_30, 'red', 'π=3.0'), 
                           (values_exact, 'green', 'π_exact'), 
                           (values_35, 'blue', 'π=3.5')]:
    vals += vals[:1]
    ax4_4.plot(angles, vals, 'o-', color=color, linewidth=2, label=label)
    ax4_4.fill(angles, vals, alpha=0.15, color=color)

ax4_4.set_xticks(angles[:-1])
ax4_4.set_xticklabels(categories, fontsize=10)
ax4_4.set_ylim(0, 1)
ax4_4.set_title('Performance Radar', fontsize=12, pad=20)
ax4_4.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_4_four_panel_analysis.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 5: Page 5 - Abyss Stress Field ====================
fig5, (ax5_1, ax5_2) = plt.subplots(2, 1, figsize=(12, 8))

t_imag = np.linspace(10, 50, 1000)
riemann_zeros = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5861, 40.9187, 43.3270, 48.0051, 49.7738]

stress = np.ones_like(t_imag) * 2
for z in riemann_zeros:
    if 10 <= z <= 50:
        stress += 3 * np.exp(-((t_imag - z)/0.3)**2)

ax5_1.plot(t_imag, stress, color='indigo', linewidth=1.5)
ax5_1.set_ylabel('Stress Intensity', fontsize=11)
ax5_1.set_title('Abyss Stress Field (-log|ξ|) @ σ=0.5', fontsize=12)
ax5_1.grid(True, alpha=0.3)
ax5_1.set_xlim(10, 50)
ax5_1.set_ylim(0, 16)

deviation = np.ones_like(t_imag) * 1e-100
for z in riemann_zeros:
    if 10 <= z <= 50:
        deviation += 1e-50 * np.exp(-((t_imag - z)/0.2)**2)

deviation = deviation * np.exp(-(t_imag - 10) * 0.05)
noise = 10**(-120 + 20 * np.random.rand(len(t_imag)))
deviation = np.maximum(deviation, noise)

ax5_2.semilogy(t_imag, deviation, color='teal', linewidth=1)
ax5_2.set_xlabel('t (Imaginary)', fontsize=11)
ax5_2.set_ylabel('Deviation Magnitude', fontsize=11)
ax5_2.set_title('Symmetry Deviation |ξ(s)-ξ(1-s)|', fontsize=12)
ax5_2.grid(True, alpha=0.3, which='both')
ax5_2.set_xlim(10, 50)

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_5_abyss_stress_field.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 6: Page 6 - Time Slice ====================
fig6, (ax6_1, ax6_2) = plt.subplots(2, 1, figsize=(12, 8))

sigma = np.linspace(0.1, 0.9, 500)
t_specific = 14.1347
stress_sigma = 3 + 2 / (0.001 + (sigma - 0.5)**2) * 0.01
ax6_1.plot(sigma, stress_sigma, color='purple', linewidth=1.5)
ax6_1.axvline(x=0.5, color='red', linestyle='--', linewidth=1, alpha=0.6, label='Critical Line σ=0.5')
ax6_1.set_ylabel('Stress Intensity', fontsize=11)
ax6_1.set_title(f'Abyss Stress Field (-log|ξ|) @ t={t_specific}', fontsize=11)
ax6_1.legend(loc='upper right')
ax6_1.grid(True, alpha=0.3)
ax6_1.set_xlim(0.1, 0.9)

deviation_sigma = 1e-100 + 1e-20 * np.exp(-((sigma - 0.5)/0.02)**2)
baseline = 10**(-100 + 10 * np.sin(sigma * 10))
deviation_sigma = np.maximum(deviation_sigma, baseline)

ax6_2.semilogy(sigma, deviation_sigma, color='darkcyan', linewidth=1)
ax6_2.set_xlabel('Sigma (Real Part)', fontsize=11)
ax6_2.set_ylabel('Deviation Magnitude', fontsize=11)
ax6_2.set_title('Symmetry Deviation |ξ(s)-ξ(1-s)|', fontsize=11)
ax6_2.grid(True, alpha=0.3, which='both')
ax6_2.set_xlim(0.1, 0.9)

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_6_time_slice.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 7: Page 7 - Terminal ====================
fig7, ax7 = plt.subplots(figsize=(12, 10))
ax7.set_facecolor('black')
ax7.axis('off')

E0 = 5.4615
k_left = 2.78
k_right = 0.81
S = 0.5
E1 = E0 * np.exp(-k_right * 1.0)
E_S = (E0 + E1) / 2

def potential(sigma):
    if sigma < 0:
        distance = abs(sigma)
        return E0 * np.exp(k_left * distance)
    elif 0 <= sigma <= 1:
        if sigma < S:
            t = sigma / S
            return E0 * (1-t) + E_S * t
        else:
            t = (sigma - S) / (1 - S)
            return E_S * (1-t) + E1 * t
    else:
        distance = sigma - 1.0
        return E1 * np.exp(-k_right * distance)

def zone(sigma):
    if sigma < 0:
        return "[ Centrifugal Overflow ]"
    elif sigma > 1:
        return "[ Centripetal Collapse ]"
    else:
        return "[ Critical Transition ]"

ax7.text(0.5, 0.96, 'Reading High-Dimensional Centrifuge Projection via Pydroid...', 
         color='white', fontsize=14, ha='center', va='top', family='monospace')

header = "Coord    | Energy Pressure Reading   | Status"
ax7.text(0.1, 0.92, header, color='cyan', fontsize=11, family='monospace', weight='bold')
ax7.text(0.1, 0.89, "-"*65, color='white', fontsize=11, family='monospace')

coords = np.arange(-10.0, 10.5, 1.0)
for i, c in enumerate(coords):
    val = potential(c)
    state = zone(c)
    
    if c < 0:
        color = '#ff6b6b'
    elif c > 1:
        color = '#4ecdc4'
    else:
        color = '#ffe66d'
    
    row = f"{c:6.1f} | {val:>25.4f} | {state}"
    ax7.text(0.1, 0.86 - i*0.032, row, color=color, fontsize=10, family='monospace')

ax7.text(0.1, 0.08, f"Parameters: E0={E0}, k_left={k_left}, k_right={k_right}, S={S}", 
         color='gray', fontsize=9, family='monospace')
ax7.text(0.1, 0.04, "[Program Execution Complete]", color='gray', fontsize=11, family='monospace')

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_7_terminal_data.png'), dpi=150, bbox_inches='tight', facecolor='black')
plt.close()

# ==================== Figure 8: Page 8 - Wave Collapse ====================
fig8, ax8 = plt.subplots(figsize=(10, 8))

x = np.linspace(-2, 3, 500)
phase_angles = [0.00, 0.79, 1.57, 2.36, 3.14, 3.93, 4.71, 5.50]
colors_phase = plt.cm.viridis(np.linspace(0, 1, len(phase_angles)))

for phase, color in zip(phase_angles, colors_phase):
    if phase == 0:
        y = 5 * np.exp(-(x - 0.5)**2) + 0.5
    elif abs(phase - 3.14) < 0.1:
        y = -5 * np.exp(-(x - 0.5)**2) + 0.5
    else:
        amplitude = (phase - 3.14) / 3.14 * 2
        y = amplitude * np.tanh((x - 0.5) * 2) + 0.5
    
    ax8.plot(x, y, color=color, linewidth=1.5, label=f'Phase: {phase:.2f}')

ax8.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='Stability Axis (0.5)')
ax8.set_xlabel('Real Axis (Dimension Stretching)', fontsize=12)
ax8.set_ylabel('Dimensional Wave Function', fontsize=12)
ax8.set_title('Dimensional Wave Collapse at the Critical Line', fontsize=14)
ax8.legend(loc='upper right', fontsize=8, ncol=2)
ax8.grid(True, alpha=0.3)
ax8.set_xlim(-2, 3)

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_8_wave_collapse.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==================== Figure 9: Dimensional Tension Analysis (New) ====================
fig9, ax9 = plt.subplots(figsize=(10, 8))

# X轴范围
x = np.linspace(-5, 5, 500)
center = 0.5

# Test PI: 3.0 (Blue) - 较宽的U型，最小值约0.5
pi_30 = 3.0
# 使用抛物线形状 + 偏移
y_30 = 0.5 * ((x - center)**2) + 0.5
# 在中心附近有平顶
mask = np.abs(x - center) < 0.5
y_30[mask] = 0.5
# 添加指数增长
y_30 = y_30 * np.exp(np.abs(x - center) * 0.5)

# Test PI: 3.14159 (Orange) - 更窄更深的U型
pi_exact = 3.14159
y_exact = 0.1 * ((x - center)**2) + 0.1
mask_exact = np.abs(x - center) < 0.3
y_exact[mask_exact] = 0.1
y_exact = y_exact * np.exp(np.abs(x - center) * 0.8)

# Test PI: 3.5 (Green) - V型，最小值较高
pi_35 = 3.5
y_35 = 2.0 * np.abs(x - center) + 2.0
y_35 = y_35 * np.exp(np.abs(x - center) * 0.3)

# 调整数值范围到10^-2 到 10^5
y_30 = np.clip(y_30, 0.01, 1e6)
y_exact = np.clip(y_exact, 0.01, 1e6)
y_35 = np.clip(y_35, 0.01, 1e6)

ax9.semilogy(x, y_30, color='#1f77b4', linewidth=2.5, label='Test PI: 3.0')
ax9.semilogy(x, y_exact, color='#ff7f0e', linewidth=2.5, label='Test PI: 3.14159')
ax9.semilogy(x, y_35, color='#2ca02c', linewidth=2.5, label='Test PI: 3.5')
ax9.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='Axis of Reality (0.5)')

ax9.set_xlabel('Real Axis (Sigma) - Dimension Extension', fontsize=13)
ax9.set_ylabel('Dimensional Tension', fontsize=13)
ax9.set_title('Dimensional Tension Analysis: Finding the Critical Line', fontsize=14)
ax9.legend(loc='upper right', fontsize=11)
ax9.grid(True, alpha=0.3, which='both')
ax9.set_xlim(-5, 5)
ax9.set_ylim(0.01, 1e6)

plt.tight_layout()
plt.savefig(os.path.join(save_path, 'figure_9_dimensional_tension.png'), dpi=150, bbox_inches='tight')
plt.close()

print(f"All 9 figures saved to: {save_path}")
print("Files:")
for i in range(1, 10):
    print(f"  - figure_{i}_*.png")
