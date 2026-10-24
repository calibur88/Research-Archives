#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-Leaf 24D Leech Lattice Centrifuge - Configurable GIF Generator
File: leech_centrifuge_configurable.py
Save Path: /storage/emulated/0/研究报告/宇宙/

========== CONFIGURATION ZONE (Edit Here) ==========
"""

# --- Rotation Settings ---
TOTAL_ROTATIONS = 3               # 转几圈？ (推荐: 2-5圈)
FRAMES_PER_ROTATION = 60          # 每圈多少帧？ (推荐: 24-72，越多越流畅)
ROTATION_SPEED = 360 / FRAMES_PER_ROTATION  # 每帧转多少度 (自动计算，别改)

# --- Animation Timing ---
INTERVAL_MS = 100                 # 每帧间隔(毫秒)，越小转得越快 (推荐: 50-200)
FPS = 30                          # 输出GIF帧率 (推荐: 10-30)

# --- Visual Effects ---
BREATHING_SPEED = 20              # 叶片呼吸闪烁速度 (推荐: 10-40)
ENERGY_FLOW_SPEED = 8             # 外圈能量流动速度 (推荐: 5-15)
PULSE_SPEED = 30                  # 星核脉动速度 (推荐: 20-40)

# --- Data Settings ---
MAX_PRIME = 3000                  # 计算素数上限 (推荐: 1000-5000)
N_SHOW_TRIPLETS = 8               # 每个叶片显示几个三素数 (推荐: 5-12)

# --- Output ---
OUTPUT_FILENAME = "Centrifuge_Rotation.gif"  # 输出文件名

"""
========== END OF CONFIGURATION ==========
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.use('Agg')
import os
import sys

# Create save directory
SAVE_DIR = "/storage/emulated/0/研究报告/宇宙/"
os.makedirs(SAVE_DIR, exist_ok=True)
GIF_PATH = os.path.join(SAVE_DIR, OUTPUT_FILENAME)

# Auto-calculated parameters
TOTAL_FRAMES = int(TOTAL_ROTATIONS * FRAMES_PER_ROTATION)

print("=" * 60)
print("3-Leaf 24D Leech Centrifuge - Configurable Generator")
print("=" * 60)
print(f"Settings:")
print(f"  Rotations: {TOTAL_ROTATIONS}")
print(f"  Frames/Rotation: {FRAMES_PER_ROTATION}")
print(f"  Total Frames: {TOTAL_FRAMES}")
print(f"  Rotation Speed: {ROTATION_SPEED:.1f} deg/frame")
print(f"  Interval: {INTERVAL_MS} ms")
print(f"  Output: {GIF_PATH}")
print("-" * 60)

# ==================== 1. Basic Algorithms ====================

def sieve_primes(n):
    """Sieve of Eratosthenes"""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i*i:n+1:i] = [False] * len(sieve[i*i:n+1:i])
    return [i for i in range(2, n + 1) if sieve[i]]

def get_prime_triplets(n_max):
    """Get Prime Triplets"""
    primes = sieve_primes(n_max)
    prime_set = set(primes)
    triplets = []
    for p in primes:
        if p > n_max - 6:
            break
        if (p+2) in prime_set and (p+6) in prime_set:
            triplets.append((p, p+2, p+6))
        elif (p+4) in prime_set and (p+6) in prime_set:
            triplets.append((p, p+4, p+6))
    return triplets

# Generate data
print("Generating prime data...")
triplets = get_prime_triplets(MAX_PRIME)
print(f"Loaded Triplets: {len(triplets)} groups")
print("Initializing graphics...")

# ==================== 2. Graphics Setup ====================

fig, ax = plt.subplots(figsize=(10, 10), dpi=100)
fig.patch.set_facecolor('black')
ax.set_facecolor('black')
ax.set_xlim(-12, 12)
ax.set_ylim(-12, 12)
ax.set_aspect('equal')
ax.axis('off')

# Static background: 24 concentric circles
for i in range(1, 25):
    circle = Circle((0, 0), i*0.45, fill=False, color='navy', alpha=0.15, linewidth=0.8)
    ax.add_patch(circle)

# Outer ring: 168 symbolic dots
n_dots = 168
outer_theta = np.linspace(0, 2*np.pi, n_dots, endpoint=False)
outer_r = 11.5
outer_dots = []
for t in outer_theta:
    dot, = ax.plot(outer_r*np.cos(t), outer_r*np.sin(t), '.', 
                   color='gold', markersize=3, alpha=0.6)
    outer_dots.append(dot)

# Center star
center_star, = ax.plot(0, 0, 'r*', markersize=25, markerfacecolor='red',
                       markeredgecolor='yellow', markeredgewidth=2, zorder=10)

# 3 Blades
blade_patches = []
blade_colors = ['#FF1744', '#00E5FF', '#76FF03']  # Red, Cyan, Green

for i in range(3):
    blade = Polygon([(0,0), (0,0), (0,0)], facecolor=blade_colors[i], 
                   alpha=0.7, edgecolor='white', linewidth=2)
    ax.add_patch(blade)
    blade_patches.append(blade)

# Triplet markers
triplet_markers = []
n_show = min(N_SHOW_TRIPLETS, len(triplets))
for i in range(3):
    leaf_markers = []
    for j in range(n_show):
        marker, = ax.plot(0, 0, 'o', color=blade_colors[i], markersize=0)
        leaf_markers.append(marker)
    triplet_markers.append(leaf_markers)

# Titles
title = ax.text(0, 11, '3-LEAF 24D LEECH CENTRIFUGE', 
               ha='center', fontsize=14, fontweight='bold', color='cyan')
subtitle = ax.text(0, -11, '[ Invincible Electric Fan - 3-Body Interference ]', 
                  ha='center', fontsize=12, color='yellow')
status_text = ax.text(-11, -10, '', fontsize=10, color='lime', fontweight='bold')

# Legend
legend_text = (
    f"Prime Triplet Separator\n"
    f"Rotations: {TOTAL_ROTATIONS}\n"
    f"Speed: {ROTATION_SPEED:.0f} deg/f\n"
    f"Load: {len(triplets)} groups"
)
ax.text(8.5, -10, legend_text, fontsize=9, color='white',
       bbox=dict(boxstyle='round', facecolor='darkslategray', alpha=0.8))

# ==================== 3. Animation Function ====================

def init():
    """Initialize"""
    return blade_patches + outer_dots + [center_star] + sum(triplet_markers, [])

def animate(frame):
    """Animation frame"""
    # Rotation
    base_angle = np.radians(frame * ROTATION_SPEED)
    angles = [base_angle, base_angle + 2*np.pi/3, base_angle + 4*np.pi/3]
    
    for i, (angle, blade_patch) in enumerate(zip(angles, blade_patches)):
        # Blade geometry
        blade_length = 9
        tip_x = blade_length * np.cos(angle)
        tip_y = blade_length * np.sin(angle)
        
        width_angle = np.pi/10
        base_r = 1.2
        x1 = base_r * np.cos(angle - width_angle)
        y1 = base_r * np.sin(angle - width_angle)
        x2 = tip_x
        y2 = tip_y
        x3 = base_r * np.cos(angle + width_angle)
        y3 = base_r * np.sin(angle + width_angle)
        
        blade_patch.set_xy([(x1, y1), (x2, y2), (x3, y3)])
        
        # Breathing effect
        alpha = 0.4 + 0.4 * np.sin(np.radians(frame * BREATHING_SPEED + i * 120))
        blade_patch.set_alpha(alpha)
        
        # Update triplet markers
        for j, marker in enumerate(triplet_markers[i]):
            if j < len(triplets):
                r_pos = 3 + j * 0.9
                offset_angle = angle + (j % 3 - 1) * 0.05
                x = r_pos * np.cos(offset_angle)
                y = r_pos * np.sin(offset_angle)
                marker.set_data([x], [y])
                marker.set_markersize(7 - j*0.5)
                marker.set_alpha(0.8)
    
    # Energy flow effect
    for idx, dot in enumerate(outer_dots):
        phase = np.radians(frame * ENERGY_FLOW_SPEED + idx * 360/168)
        intensity = 0.3 + 0.7 * (0.5 + 0.5 * np.sin(phase))
        dot.set_alpha(intensity)
        if intensity > 0.9:
            dot.set_markersize(5)
            dot.set_color('yellow')
        else:
            dot.set_markersize(3)
            dot.set_color('gold')
    
    # Star pulsing
    pulse = 20 + 6 * np.sin(np.radians(frame * PULSE_SPEED))
    center_star.set_markersize(pulse)
    
    # Status update
    current_cycle = frame // FRAMES_PER_ROTATION + 1
    current_deg = (frame * ROTATION_SPEED) % 360
    status_text.set_text(f'Frame: {frame}/{TOTAL_FRAMES}\n'
                        f'Cycle: {current_cycle}/{TOTAL_ROTATIONS}\n'
                        f'Deg: {current_deg:.0f}\n'
                        f'Status: Stable')
    
    return blade_patches + outer_dots + [center_star] + sum(triplet_markers, []) + [status_text]

# ==================== 4. Generate GIF ====================

print(f"Generating {TOTAL_FRAMES} frames...")
print("Please wait...")

try:
    anim = FuncAnimation(fig, animate, init_func=init,
                        frames=TOTAL_FRAMES, interval=INTERVAL_MS, blit=False)
    
    anim.save(GIF_PATH, writer='pillow', fps=FPS, dpi=100)
    plt.close()
    
    print("-" * 60)
    print("SUCCESS!")
    print(f"Location: {GIF_PATH}")
    print(f"File Size: {os.path.getsize(GIF_PATH)/1024:.1f} KB")
    print("-" * 60)
    print("Parameters:")
    print(f"  Total Frames: {TOTAL_FRAMES}")
    print(f"  Duration: ~{TOTAL_FRAMES * INTERVAL_MS / 1000:.1f} seconds")
    print(f"  Cycles: {TOTAL_ROTATIONS}")
    print("=" * 60)
    
except Exception as e:
    print(f"FAILED: {str(e)}")
    print("Tip: Reduce TOTAL_FRAMES or FRAMES_PER_ROTATION if out of memory")
    sys.exit(1)

print(f"\nCheck: Internal Storage/研究报告/宇宙/{OUTPUT_FILENAME}")
