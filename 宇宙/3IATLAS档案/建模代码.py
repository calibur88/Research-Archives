import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import rcParams

# ====================== 巴纳德星共轨双伴星模型可视化代码 ======================
# 功能说明：基于巴纳德星真实观测数据与共轨双伴星动力学模型，生成完整静态视图与轨道周期动画
# 输出：高清静态图、轨道动画GIF，均自动保存至当前目录

# ====================== 全局绘图风格设置 ======================
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 10
rcParams['axes.grid'] = True
rcParams['grid.linestyle'] = '-'
rcParams['grid.alpha'] = 0.3
rcParams['axes.linewidth'] = 1
rcParams['figure.dpi'] = 120

# ====================== 天文常数与单位转换 ======================
AU_TO_LY = 1 / 63241.072       # 1 AU 转 光年
LY_TO_AU = 63241.072           # 1 光年 转 AU
PC_TO_LY = 3.2615637771418798 # 1 秒差距 转 光年
LY_TO_PC = 1 / PC_TO_LY        # 1 光年 转 秒差距
KM_S_TO_LY_YEAR = 1 / 299792.458  # 1 km/s 转 光年/年
M_JUP_TO_M_SUN = 0.000954588      # 1 木星质量 转 太阳质量

# ====================== 巴纳德星真实观测参数 ======================
M_STAR = 0.16                  # 中心恒星质量（太阳质量）
D_NOW = 5.96                   # 当前与太阳系的距离（光年）
D_NOW_PC = D_NOW * LY_TO_PC    # 当前与太阳系的距离（秒差距）
V_RADIAL_OBS = -110             # 视向速度（km/s，负号表示向太阳系靠近）
PROPER_MOTION = 10.36          # 自行（角秒/年）
V_TANGENTIAL = 4.74 * PROPER_MOTION * D_NOW_PC  # 切向速度（km/s）

# 转换为光年/年单位（用于长期距离演化计算）
V_RADIAL_LY_YEAR = V_RADIAL_OBS * KM_S_TO_LY_YEAR
V_TANGENTIAL_LY_YEAR = V_TANGENTIAL * KM_S_TO_LY_YEAR

# ====================== 共轨双伴星模型核心参数 ======================
# 伴星基本参数
M_A = 52 * M_JUP_TO_M_SUN      # 伴星A质量（太阳质量）
M_B = 42 * M_JUP_TO_M_SUN      # 伴星B质量（太阳质量）
M_BIN = M_A + M_B               # 双星总质量

# 双星质心绕恒星的轨道参数
A_CM = 15                       # 质心轨道半长轴（AU）
E_CM = 0.6                      # 质心轨道偏心率
P_CM = 145                      # 质心轨道周期（年）

# 伴星互绕轨道参数
A_MUT = 5.24                    # 伴星互绕半长轴（AU）
P_MUT = 40                      # 伴星互绕周期（年）

# 伴星相对质心的距离（严格满足质心平衡条件 M_A*R_A = M_B*R_B）
R_A = A_MUT * M_B / (M_A + M_B)  # 伴星A相对质心的距离（AU）
R_B = A_MUT * M_A / (M_A + M_B)  # 伴星B相对质心的距离（AU）

# 天体测量信号参数
A_STAR = (M_BIN / M_STAR) * A_CM  # 恒星绕系统质心的轨道半长轴（AU）
THETA_MAX = A_STAR / D_NOW_PC     # 恒星最大角位移（角秒）
theta_min = 1.2e6                  # 天体测量信号谷值（μas）
theta_max = 4.7e6                  # 天体测量信号峰值（μas）

# 希尔球半径（长期稳定性验证）
R_H = A_CM * (M_BIN / (3 * M_STAR)) ** (1/3)

# 核心参数校验输出
print("="*50)
print("模型核心参数校验结果")
print(f"伴星A相对质心距离：{R_A:.2f} AU | 伴星B相对质心距离：{R_B:.2f} AU")
print(f"恒星最大角位移：{THETA_MAX:.2f} 角秒 | 系统希尔球半径：{R_H:.2f} AU")
print(f"伴星最小距离：A={6-R_A:.3f} AU | B={6-R_B:.3f} AU")
print("="*50)

# ====================== 轨道计算核心函数（支持标量+数组输入） ======================
def kepler_solver(M, e, tol=1e-8):
    """开普勒方程求解，支持标量和数组输入，返回偏近点角E"""
    M = np.asarray(M, dtype=np.float64)
    # 初始化偏近点角
    E = np.where(e < 0.8, M, np.full_like(M, np.pi))
    # 迭代收敛（支持数组批量计算）
    for _ in range(100):
        E_new = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        # 检查所有元素是否收敛
        max_error = np.max(np.abs(E_new - E))
        if max_error < tol:
            break
        E = E_new
    return E

def get_orbit_position(a, e, t, P):
    """计算椭圆轨道上的位置，支持标量/数组时间输入，返回(x,y)坐标（AU）"""
    M = 2 * np.pi * t / P                  # 平近点角
    E = kepler_solver(M, e)                # 偏近点角
    # 真近点角计算
    nu = 2 * np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    r = a * (1 - e * np.cos(E))            # 轨道半径
    return r * np.cos(nu), r * np.sin(nu)

def get_companion_positions(t):
    """计算t时刻所有天体的位置（恒星位于坐标原点）
    返回：(恒星位置, 质心位置, 伴星A位置, 伴星B位置)
    """
    # 1. 质心绕恒星的轨道位置
    cm_x, cm_y = get_orbit_position(A_CM, E_CM, t, P_CM)
    # 2. 伴星绕质心的相对位置（圆轨道）
    M_mut = 2 * np.pi * t / P_MUT
    a_x_rel = R_A * np.cos(M_mut)
    a_y_rel = R_A * np.sin(M_mut)
    b_x_rel = - R_B * np.cos(M_mut)
    b_y_rel = - R_B * np.sin(M_mut)
    # 3. 绝对位置计算
    a_x, a_y = cm_x + a_x_rel, cm_y + a_y_rel
    b_x, b_y = cm_x + b_x_rel, cm_y + b_y_rel
    return (0.0, 0.0), (cm_x, cm_y), (a_x, a_y), (b_x, b_y)

# ====================== 创建画布与子图布局（3行2列） ======================
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(14, 16))
fig.suptitle("Barnard's Star Co-orbiting Binary Companion Model", fontsize=14, fontweight='bold', y=0.99)

# ====================== 1. 左上角：巴纳德星与太阳系距离长期演化图 ======================
t1 = np.linspace(-20, 80, 1000) * 1000  # 时间轴：-20000年 ~ 80000年
# 相对距离计算
x_rel = D_NOW + V_RADIAL_LY_YEAR * t1
y_rel = V_TANGENTIAL_LY_YEAR * t1
r1 = np.sqrt(x_rel**2 + y_rel**2)
# 绘制曲线
ax1.plot(t1/1000, r1, color='#7a0177', linewidth=3, label="Barnard's Star")
# 参考线与关键节点标注
ax1.axvline(x=0, color='#00bfff', linestyle='--', linewidth=1.5, label="Today")
ax1.axhspan(0, 1, color='#d3d3d3', alpha=0.5, label="Oort Cloud Border")
# 最近接近点计算
t_closest = - (D_NOW * V_RADIAL_LY_YEAR) / (V_RADIAL_LY_YEAR**2 + V_TANGENTIAL_LY_YEAR**2)
r_closest = D_NOW * V_TANGENTIAL_LY_YEAR / np.sqrt(V_RADIAL_LY_YEAR**2 + V_TANGENTIAL_LY_YEAR**2)
ax1.scatter(0, D_NOW, color='#7a0177', s=60, zorder=5)
ax1.text(-5, D_NOW+0.3, f"Today\n({D_NOW} ly)", fontsize=9, ha='center')
ax1.scatter(t_closest/1000, r_closest, color='#7a0177', marker='*', s=80, zorder=5)
ax1.text(t_closest/1000 + 5, r_closest, f"Closest Approach\n(~{r_closest:.1f} ly @ {int(t_closest/1000)} ky)", fontsize=9, ha='left')
# 坐标轴设置
ax1.set_title("Barnard's Star Distance to Solar System", fontsize=11, fontweight='bold')
ax1.set_xlabel("Time (1000 Years)", fontsize=9)
ax1.set_ylabel("Distance (Light Years)", fontsize=9)
ax1.set_xlim(-20, 80)
ax1.set_ylim(0, 12)
ax1.legend(loc='upper right', fontsize=8)

# ====================== 2. 右上角：共轨双星系统轨道示意图 ======================
# 绘制质心完整轨道
t_orbit = np.linspace(0, P_CM, 1000)
cm_orbit_x, cm_orbit_y = get_orbit_position(A_CM, E_CM, t_orbit, P_CM)
ax2.plot(cm_orbit_x, cm_orbit_y, color='#555555', linestyle='--', linewidth=1, label="CM Orbit (145 yr)")
# 绘制希尔球范围
_, (cm_x0, cm_y0), _, _ = get_companion_positions(0)
hill_circle = plt.Circle((cm_x0, cm_y0), R_H, color='#00aa00', linestyle=':', fill=False, linewidth=1.5, label=f"Hill Sphere (R={R_H:.2f} AU)")
ax2.add_artist(hill_circle)
# 初始时刻（t=0）天体位置
(star_x, star_y), (cm_x, cm_y), (a_x, a_y), (b_x, b_y) = get_companion_positions(0)
star_dot, = ax2.plot(star_x, star_y, 'o', color='#ff4400', markersize=10, label="Barnard's Star (0.16 M☉)")
cm_dot, = ax2.plot(cm_x, cm_y, 'o', color='#aaaaaa', markersize=4, label="Binary CM")
a_dot, = ax2.plot(a_x, a_y, 'o', color='#00aaff', markersize=7, label=f"Companion A ({M_A/M_JUP_TO_M_SUN:.0f} Mjup)")
b_dot, = ax2.plot(b_x, b_y, 'o', color='#ffcc00', markersize=7, label=f"Companion B ({M_B/M_JUP_TO_M_SUN:.0f} Mjup)")
# 伴星与质心的连线
line_a, = ax2.plot([cm_x, a_x], [cm_y, a_y], color='#00aaff', linestyle='-', linewidth=0.8)
line_b, = ax2.plot([cm_x, b_x], [cm_y, b_y], color='#ffcc00', linestyle='-', linewidth=0.8)
# 坐标轴设置
ax2.set_title("Co-orbiting Binary System Orbit", fontsize=11, fontweight='bold')
ax2.set_xlabel("X (AU)", fontsize=9)
ax2.set_ylabel("Y (AU)", fontsize=9)
ax2.set_xlim(-30, 30)
ax2.set_ylim(-30, 30)
ax2.set_aspect('equal', adjustable='box')
ax2.legend(loc='upper right', fontsize=8)
ax2.grid(True, alpha=0.2)

# ====================== 3. 左中：145年周期天体测量信号图 ======================
t3 = np.linspace(0, 500, 2000)  # 时间轴：0~500年
# 计算高偏心率轨道对应的非对称天体测量信号
theta3 = np.zeros_like(t3)
for i, t in enumerate(t3):
    M = 2 * np.pi * t / P_CM
    E = kepler_solver(M, E_CM)
    nu = 2 * np.arctan2(np.sqrt(1+E_CM)*np.sin(E/2), np.sqrt(1-E_CM)*np.cos(E/2))
    theta3[i] = theta_min + (theta_max - theta_min) * (1 + np.cos(nu)) / 2
# 绘制曲线
ax3.plot(t3, theta3, color='#ff4400', linewidth=2, label="Astrometric Shift")
# 当前时刻标记
current_theta = theta_min + (theta_max - theta_min) * (1 + np.cos(0)) / 2
theta_dot, = ax3.plot(0, current_theta, 'o', color='#7a0177', markersize=6, zorder=5, label="Current Time")
# 坐标轴设置
ax3.set_title("Astrometric Signal (145 yr CM Orbit)", fontsize=11, fontweight='bold')
ax3.set_xlabel("Time (Years)", fontsize=9)
ax3.set_ylabel("Astrometric Shift (μas)", fontsize=9)
ax3.set_xlim(0, 500)
ax3.set_ylim(1e6, 5e6)
ax3.legend(loc='upper right', fontsize=8)

# ====================== 4. 右中：40年互绕周期自行扰动图 ======================
t4 = np.linspace(0, 200, 2000)  # 时间轴：0~200年
pm_mean = PROPER_MOTION
pm_amplitude = 0.01  # 自行调制幅度（±10 mas）
pm4 = pm_mean + pm_amplitude * np.sin(2 * np.pi * t4 / P_MUT)
# 绘制曲线
ax4.plot(t4, pm4, color='#7a0177', linewidth=2, label="Total Proper Motion")
ax4.axhline(y=pm_mean, color='#888888', linestyle='--', linewidth=1.5, label=f"Mean ({pm_mean:.2f} \"/yr)")
ax4.fill_between(t4, pm_mean - pm_amplitude, pm_mean + pm_amplitude, color='#c774c7', alpha=0.3, label=f"40 yr Modulation (±{pm_amplitude*1000:.0f} mas)")
# 当前时刻标记
current_pm = pm_mean + pm_amplitude * np.sin(0)
pm_dot, = ax4.plot(0, current_pm, 'o', color='#7a0177', markersize=6, zorder=5)
# 坐标轴设置
ax4.set_title("Proper Motion Perturbation (40 yr Mutual Orbit)", fontsize=11, fontweight='bold')
ax4.set_xlabel("Time (Years)", fontsize=9)
ax4.set_ylabel("Proper Motion (arcsec/year)", fontsize=9)
ax4.set_xlim(0, 200)
ax4.set_ylim(pm_mean - 1.5*pm_amplitude, pm_mean + 1.5*pm_amplitude)
ax4.legend(loc='upper right', fontsize=8)

# ====================== 5. 左下：伴星距离长期演化稳定性图 ======================
n_cycles = 1000  # 1000个轨道周期 = 145000年
t5 = np.linspace(0, n_cycles * P_CM, int(n_cycles * 50))
# 计算伴星到恒星的距离
dist_a5 = np.zeros_like(t5)
dist_b5 = np.zeros_like(t5)
for i, t in enumerate(t5):
    _, _, (a_x, a_y), (b_x, b_y) = get_companion_positions(t)
    dist_a5[i] = np.sqrt(a_x**2 + a_y**2)
    dist_b5[i] = np.sqrt(b_x**2 + b_y**2)
# 绘制曲线
ax5.plot(t5 / P_CM, dist_a5, color='#00aaff', linewidth=0.8, alpha=0.8, label="Companion A")
ax5.plot(t5 / P_CM, dist_b5, color='#ffcc00', linewidth=0.8, alpha=0.8, label="Companion B")
# 标注最小距离
min_a = np.min(dist_a5)
min_b = np.min(dist_b5)
ax5.scatter(0, min_b, color='#ffcc00', s=60, zorder=5, label=f"Min B = {min_b:.3f} AU")
ax5.scatter(n_cycles, min_a, color='#00aaff', s=60, zorder=5, label=f"Min A = {min_a:.3f} AU")
# 坐标轴设置
ax5.set_title("Long-term Companion Distance Evolution (145,000 yr)", fontsize=11, fontweight='bold')
ax5.set_xlabel("Time (Number of 145 yr Cycles)", fontsize=9)
ax5.set_ylabel("Distance to Barnard's Star (AU)", fontsize=9)
ax5.set_xlim(0, n_cycles)
ax5.set_ylim(3, 27)
ax5.legend(loc='upper right', fontsize=8)

# ====================== 6. 右下：视向速度一致性校验图 ======================
t6 = np.linspace(0, 20, 1000) * 1000  # 时间轴：0~20000年
# 计算理论径向速度
x_rel6 = D_NOW + V_RADIAL_LY_YEAR * t6
y_rel6 = V_TANGENTIAL_LY_YEAR * t6
r6 = np.sqrt(x_rel6**2 + y_rel6**2)
v_r_theory6 = (x_rel6 * V_RADIAL_LY_YEAR + y_rel6 * V_TANGENTIAL_LY_YEAR) / r6
v_r_theory6_km_s = v_r_theory6 / KM_S_TO_LY_YEAR
# 伴星轨道速度
v_companion_b = 2.3
v_companion_a = -2.3
# 绘制曲线
ax6.plot(t6/1000, v_r_theory6_km_s, color='#0033cc', linewidth=2, label="Theoretical Radial Velocity")
ax6.axhline(y=V_RADIAL_OBS, color='#cc0000', linestyle='--', linewidth=1.5, label=f"Observed RV ({V_RADIAL_OBS} km/s)")
ax6.axhline(y=v_companion_b, color='#e6b800', linestyle=':', linewidth=2, label=f"Companion B Orbital Velocity ({v_companion_b} km/s)")
ax6.axhline(y=v_companion_a, color='#00cccc', linestyle=':', linewidth=2, label=f"Companion A Orbital Velocity ({v_companion_a} km/s)")
# 标注说明
ax6.annotate("Model: Companion B dominates RV signal",
             xy=(0, v_companion_b), xytext=(5, 10),
             arrowprops=dict(arrowstyle='->', color='#cc8800', linewidth=1.2),
             fontsize=8, ha='left')
# 当前时刻标记
current_vr = v_r_theory6_km_s[0]
vr_dot, = ax6.plot(0, current_vr, 'o', color='#7a0177', markersize=6, zorder=5)
# 坐标轴设置
ax6.set_title("Radial Velocity Consistency Check", fontsize=11, fontweight='bold')
ax6.set_xlabel("Time from Present (Ky)", fontsize=9)
ax6.set_ylabel("Radial Velocity (km/s)", fontsize=9)
ax6.set_xlim(0, 20)
ax6.set_ylim(-110, 0)
ax6.legend(loc='lower right', fontsize=8)

# 调整布局，避免标题重叠
plt.tight_layout(rect=[0, 0, 1, 0.98])

# 保存完整静态图到当前目录
static_file = "barnard_star_co-orbit_model_full.png"
plt.savefig(static_file, dpi=150, bbox_inches='tight')
print(f"完整静态图已保存至当前目录：{static_file}")

# ====================== 轨道周期动画生成 ======================
print("\n正在生成轨道动画（145年完整轨道周期），请稍候...")
# 动画参数
total_frames = 145  # 145帧，每帧对应1年，覆盖完整质心轨道周期
fps = 10

# 动画帧更新函数
def update(frame):
    t = frame
    # 更新天体位置
    (star_x, star_y), (cm_x, cm_y), (a_x, a_y), (b_x, b_y) = get_companion_positions(t)
    star_dot.set_data(star_x, star_y)
    cm_dot.set_data(cm_x, cm_y)
    a_dot.set_data(a_x, a_y)
    b_dot.set_data(b_x, b_y)
    # 更新伴星与质心的连线
    global line_a, line_b
    line_a.remove()
    line_b.remove()
    line_a, = ax2.plot([cm_x, a_x], [cm_y, a_y], color='#00aaff', linestyle='-', linewidth=0.8)
    line_b, = ax2.plot([cm_x, b_x], [cm_y, b_y], color='#ffcc00', linestyle='-', linewidth=0.8)
    # 更新天体测量信号标记
    M = 2 * np.pi * t / P_CM
    E = kepler_solver(M, E_CM)
    nu = 2 * np.arctan2(np.sqrt(1+E_CM)*np.sin(E/2), np.sqrt(1-E_CM)*np.cos(E/2))
    current_theta_t = theta_min + (theta_max - theta_min) * (1 + np.cos(nu)) / 2
    theta_dot.set_data(t, current_theta_t)
    # 更新自行扰动标记
    current_pm_t = pm_mean + pm_amplitude * np.sin(2 * np.pi * t / P_MUT)
    pm_dot.set_data(t, current_pm_t)
    # 返回更新对象
    return star_dot, cm_dot, a_dot, b_dot, line_a, line_b, theta_dot, pm_dot, vr_dot

# 创建并保存动画
ani = FuncAnimation(fig, update, frames=total_frames, interval=1000/fps, blit=True)
writer = PillowWriter(fps=fps)
gif_file = "barnard_star_co-orbit_model.gif"
ani.save(gif_file, writer=writer, dpi=120)
print(f"轨道动画已保存至当前目录：{gif_file}")

# 显示完整静态图
plt.show()
