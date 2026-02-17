#! python
# -*- coding:utf-8 -*-
###
# --------------------------------------------------------------------------------
# 文件名: 图.py
# 创建时间: 2026-02-17 21:42:57 Tue
# 说明:
# 作者: Calibur88
# 主机: LAPTOP-D92A7OL2
# --------------------------------------------------------------------------------
# 最后编辑作者: Calibur88
# 最后修改时间: 2026-02-17 21:46:06 Tue
# --------------------------------------------------------------------------------
# Copyright (c) 2026 Calibur88
# --------------------------------------------------------------------------------
# 更新历史:
# --------------------------------------------------------------------------------
# 时间      		作者		信息
# ----------		---		------------------------------------------------------
###
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(42)
fig = plt.figure(figsize=(16, 12))

# ===== Exp 3: 波动真实 vs 平滑幻觉 =====
ax1 = plt.subplot(2, 2, 1)
time = np.linspace(0, 300, 300)
smooth = 0.6 - 0.1 * (1 - np.exp(-0.01 * time))
noise = np.cumsum(np.random.normal(0, 0.003, 300))
real = smooth + noise * np.exp(-0.005 * time) - np.linspace(0, 0.02, 300)

ax1.plot(time, smooth, "--", color="#5B9BD5", linewidth=2, label="Smooth (Discarded)")
ax1.plot(time, real, "-", color="#ED7D31", linewidth=1.5, label="Real (With Xi)")
ax1.axhline(y=0.5, color="black", linestyle=":", linewidth=2, label="Attractor")
ax1.set_title(
    "Exp 3: Fluctuation is Real, Smooth is Illusion", fontsize=12, fontweight="bold"
)
ax1.set_xlabel("Time")
ax1.set_ylabel("State")
ax1.legend()
ax1.grid(True, alpha=0.3)

# ===== Exp 4: 12维球面投影 =====
ax2 = plt.subplot(2, 2, 2)
np.random.seed(123)
points = np.random.randn(400, 12)
points = points / np.linalg.norm(points, axis=1, keepdims=True)
proj = points @ (np.random.randn(12, 2) * 0.3)

ax2.scatter(
    proj[:, 0],
    proj[:, 1],
    c="#5F9EA0",
    alpha=0.6,
    s=30,
    edgecolors="white",
    linewidth=0.5,
)
ax2.axvline(x=0, color="red", linestyle="--", linewidth=2, label="Equator")
ax2.set_title("Exp 4: 12D Sphere → 2D Line Clustering", fontsize=12, fontweight="bold")
ax2.set_xlabel("Projection X")
ax2.set_ylabel("Projection Y")
ax2.legend()
ax2.grid(True, alpha=0.3)

# ===== Exp 1: 刚性吸引子 =====
ax3 = plt.subplot(2, 2, 3)
t = np.linspace(0, 150, 200)
colors = plt.cm.viridis(np.linspace(0, 1, 12))

for ic in np.linspace(0.1, 0.9, 12):
    traj = 0.5 + (ic - 0.5) * np.exp(-0.15 * t)
    noise = np.random.normal(0, 0.002, 200) * np.exp(-0.1 * t)
    ax3.plot(t, traj + noise, color=colors[len(ax3.lines) // 2], linewidth=2, alpha=0.8)

ax3.axhline(y=0.5, color="red", linestyle="--", linewidth=3, label="Eigenaxis = 0.5")
ax3.set_title(
    "Exp 1: Rigid Attractor - All Roads Lead to 1/2", fontsize=12, fontweight="bold"
)
ax3.set_xlabel("Meta-Time")
ax3.set_ylabel("State")
ax3.legend()
ax3.grid(True, alpha=0.3)

# ===== Exp 2: 观测锁定效应 =====
ax4 = plt.subplot(2, 2, 4)
ax4.axis("off")
inner = ax4.get_subplotspec().subgridspec(2, 2, hspace=0.3, wspace=0.3)

lambdas = [0.1, 0.02, 0.005, 0.001]
data = [
    [0.6, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515],
    [0.6, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50],
    [0.6, 0.502, 0.502, 0.502, 0.502, 0.502, 0.502, 0.502, 0.502],
    [0.6, 0.52, 0.502, 0.502, 0.502, 0.502, 0.502, 0.502, 0.502],
]

for idx, (lam, vals) in enumerate(zip(lambdas, data)):
    ax = fig.add_subplot(inner[idx // 2, idx % 2])
    ax.plot(range(9), vals, "o-", color="#2E75B6", markersize=5)
    ax.axhline(y=0.5, color="red", linestyle="--", linewidth=1.5)
    ax.set_title(f"λ = {lam} (Sparse Obs)", fontsize=10)
    ax.set_ylim(0.48, 0.65)
    ax.grid(True, alpha=0.3)

fig.suptitle(
    "Exp 2: Observational Locking Effect", fontsize=14, fontweight="bold", y=0.52
)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("高维投影.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.show()
