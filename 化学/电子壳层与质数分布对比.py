#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electron Shells vs Prime Distribution - Showing Incomplete Shell 10
Save path: /storage/emulated/0/研究报告/化学/
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import matplotlib.animation as animation
import os

# ==================== Configuration ====================
SAVE_DIR = "/storage/emulated/0/研究报告/化学/"
os.makedirs(SAVE_DIR, exist_ok=True)
OUTPUT_FILE = os.path.join(SAVE_DIR, "Electron_vs_Prime_Shell10_Incomplete.gif")

MAX_N = 10
FRAMES_PER_SHELL = 12
ENDING_FRAMES = 240
FPS = 20
INTERVAL = 50

# ==================== Prime Generation ====================
def sieve_primes(limit):
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = [False] * (((limit - i*i) // i) + 1)
    return [i for i in range(2, limit + 1) if sieve[i]]

def first_n_primes(n):
    if n < 6:
        limit = 20
    else:
        limit = int(n * (np.log(n) + np.log(np.log(n)) + 2)) + 10000
    all_primes = sieve_primes(limit)
    return all_primes[:n]

print("Generating primes...")
primes = first_n_primes(600)
print(f"Generated {len(primes)} primes")

shell_capacities = [2 * n**2 for n in range(1, MAX_N + 1)]
total_capacity = sum(shell_capacities)

# Assign primes to shells
shell_prime_indices = []
idx = 0
for cap in shell_capacities:
    end_idx = min(idx + cap, len(primes))
    shell_prime_indices.append(list(range(idx, end_idx)))
    idx = end_idx

# Calculate missing count for shell 10
filled_shell10 = len(shell_prime_indices[9])
missing_shell10 = shell_capacities[9] - filled_shell10

print(f"Shell 10: {filled_shell10}/{shell_capacities[9]} filled, {missing_shell10} vacancies")

# ==================== Setup Figure ====================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), dpi=100)
fig.patch.set_facecolor('black')

for ax in (ax1, ax2):
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_facecolor('black')

colors_electron = plt.cm.tab10(np.linspace(0, 1, MAX_N))
colors_prime = plt.cm.viridis(np.linspace(0, 1, MAX_N))

# Storage
shell_rings_left = []
shell_rings_right = []
electron_dots = [[] for _ in range(MAX_N)]
prime_dots = [[] for _ in range(MAX_N)]  # Filled primes
vacancy_dots = [[] for _ in range(MAX_N)]  # ADDED: Vacancy markers for shell 10
labels_left = []
labels_right = []

for n in range(1, MAX_N + 1):
    radius = n * 1.2
    
    # Rings
    ring_left = Circle((0, 0), radius, fill=False, 
                       color=colors_electron[n-1], linewidth=2, alpha=0, linestyle='--')
    ax1.add_patch(ring_left)
    shell_rings_left.append(ring_left)
    
    ring_right = Circle((0, 0), radius, fill=False, 
                        color=colors_prime[n-1], linewidth=2, alpha=0, linestyle=':')
    ax2.add_patch(ring_right)
    shell_rings_right.append(ring_right)
    
    # Labels with vacancy info for shell 10
    label_l = ax1.text(radius*0.7, radius*0.7, f'n={n}\n2n²={shell_capacities[n-1]}', 
                       fontsize=9, color='white', ha='center', fontweight='bold', alpha=0)
    labels_left.append(label_l)
    
    if n == 10:
        label_text = f'Shell {n}\n{filled_shell10}/{shell_capacities[n-1]} primes\n⚠ {missing_shell10} vacancies'
        label_color = 'red'
    else:
        actual_count = len(shell_prime_indices[n-1])
        label_text = f'Shell {n}\n{actual_count}/{shell_capacities[n-1]} primes'
        label_color = 'white'
    
    label_r = ax2.text(radius*0.7, -radius*0.7, label_text, 
                       fontsize=9, color=label_color, ha='center', fontweight='bold', alpha=0)
    labels_right.append(label_r)

# Center nuclei
nucleus_left = ax1.plot(0, 0, 'o', color='yellow', markersize=20, 
                        markeredgecolor='orange', markeredgewidth=2, alpha=0)[0]
nucleus_right = ax2.plot(0, 0, 'o', color='red', markersize=20, 
                         markeredgecolor='pink', markeredgewidth=2, alpha=0)[0]

# Titles
title = fig.text(0.5, 0.95, '', ha='center', fontsize=16, color='white', fontweight='bold')
subtitle = fig.text(0.5, 0.02, '', ha='center', fontsize=10, color='gray')

# ==================== Animation ====================
def animate(frame):
    total_build_frames = MAX_N * FRAMES_PER_SHELL
    
    if frame < total_build_frames:
        current_shell = frame // FRAMES_PER_SHELL
        phase = frame % FRAMES_PER_SHELL
        
        for n in range(current_shell + 1):
            radius = (n + 1) * 1.2
            capacity = shell_capacities[n]
            
            # Fade in rings
            if n == current_shell:
                alpha = min(phase / 3, 1) * 0.8 if phase < 3 else 0.8
            else:
                alpha = 0.8
            
            shell_rings_left[n].set_alpha(alpha)
            shell_rings_right[n].set_alpha(alpha)
            
            # Nuclei
            if n == 0 and current_shell == 0:
                if 3 <= phase < 6:
                    nucleus_left.set_alpha((phase - 3) / 3)
                    nucleus_right.set_alpha((phase - 3) / 3)
            elif n == 0:
                nucleus_left.set_alpha(1)
                nucleus_right.set_alpha(1)
            
            # Electrons (left side)
            if n == current_shell and phase >= 6:
                progress = min(1.0, (phase - 6) / 4)
                n_electrons = min(capacity, 50)
                n_show = int(n_electrons * progress)
                
                for dot in electron_dots[n]:
                    dot.remove()
                electron_dots[n] = []
                
                if n_show > 0:
                    angles = np.linspace(0, 2*np.pi, n_electrons, endpoint=False)[:n_show]
                    for angle in angles:
                        x = radius * np.cos(angle)
                        y = radius * np.sin(angle)
                        dot = ax1.plot(x, y, 'o', color=colors_electron[n], markersize=4, alpha=0.9)[0]
                        electron_dots[n].append(dot)
            
            # Primes (right side) - filled positions
            if n == current_shell and phase >= 6:
                progress = min(1.0, (phase - 6) / 4)
                prime_indices = shell_prime_indices[n]
                n_show = int(len(prime_indices) * progress)
                
                for dot in prime_dots[n]:
                    dot.remove()
                prime_dots[n] = []
                
                if n_show > 0:
                    for i in range(n_show):
                        p = primes[prime_indices[i]]
                        angle = (np.log(p) * 2) % (2 * np.pi)
                        x = radius * np.cos(angle)
                        y = radius * np.sin(angle)
                        dot = ax2.plot(x, y, 'o', color=colors_prime[n], markersize=5, alpha=0.9)[0]
                        prime_dots[n].append(dot)
                
                # For shell 10, show vacancies (hollow circles) - appear after filled ones
                if n == 9 and phase >= 9:  # Shell 10 is index 9
                    vac_progress = (phase - 9) / 3
                    n_vac_show = int(min(50, missing_shell10) * vac_progress)  # Show max 50 vacancies to avoid clutter
                    
                    for dot in vacancy_dots[n]:
                        dot.remove()
                    vacancy_dots[n] = []
                    
                    if n_vac_show > 0:
                        # Place vacancies at regular intervals between filled primes
                        vac_angles = np.linspace(0, 2*np.pi, n_vac_show, endpoint=False)
                        for angle in vac_angles:
                            x = radius * np.cos(angle)
                            y = radius * np.sin(angle)
                            # Hollow circle for vacancy
                            dot = ax2.plot(x, y, 'o', color='red', markersize=5, alpha=0.6, 
                                          markerfacecolor='none', markeredgewidth=1)[0]
                            vacancy_dots[n].append(dot)
            
            # Labels
            if n == current_shell and phase >= 11:
                label_alpha = (phase - 11) / 3
                labels_left[n].set_alpha(label_alpha)
                labels_right[n].set_alpha(label_alpha)
    
    else:
        # Rotation phase
        ending_progress = (frame - total_build_frames) / ENDING_FRAMES
        total_angle = 2 * 2 * np.pi  # 2 rotations
        current_rotation = ending_progress * total_angle
        
        for n in range(MAX_N):
            radius = (n + 1) * 1.2
            
            # Rotate electrons
            for i, dot in enumerate(electron_dots[n]):
                if dot.get_alpha() > 0:
                    angles = np.linspace(0, 2*np.pi, min(shell_capacities[n], 50), endpoint=False)
                    if i < len(angles):
                        angle = angles[i] + current_rotation
                        x = radius * np.cos(angle)
                        y = radius * np.sin(angle)
                        dot.set_xdata([x])
                        dot.set_ydata([y])
            
            # Rotate filled primes
            for i, dot in enumerate(prime_dots[n]):
                if dot.get_alpha() > 0 and i < len(shell_prime_indices[n]):
                    p = primes[shell_prime_indices[n][i]]
                    base_angle = (np.log(p) * 2) % (2 * np.pi)
                    angle = base_angle + current_rotation * 0.5
                    x = radius * np.cos(angle)
                    y = radius * np.sin(angle)
                    dot.set_xdata([x])
                    dot.set_ydata([y])
            
            # Rotate vacancies (shell 10 only)
            if n == 9:  # Shell 10
                for i, dot in enumerate(vacancy_dots[n]):
                    if dot.get_alpha() > 0:
                        vac_angles = np.linspace(0, 2*np.pi, len(vacancy_dots[n]), endpoint=False)
                        if i < len(vac_angles):
                            angle = vac_angles[i] + current_rotation * 0.5
                            x = radius * np.cos(angle)
                            y = radius * np.sin(angle)
                            dot.set_xdata([x])
                            dot.set_ydata([y])
    
    # Update titles
    if frame < total_build_frames:
        current_n = frame // FRAMES_PER_SHELL + 1
        if current_n == 10:
            title.set_text(f'Building Shell 10 (INCOMPLETE: only {filled_shell10}/{shell_capacities[9]} primes)')
        else:
            title.set_text(f'Electron Shells vs Prime Distribution - Building Shell {current_n}/{MAX_N}')
        subtitle.set_text('Left: Bohr Model (2n² capacity) | Right: Prime Distribution by Shell Capacity')
    else:
        title.set_text('Shell 10: Only 30/200 primes (15% filled) - 170 vacancies')
        subtitle.set_text(f'Shells 1-9: Fully filled (570 primes) | Shell 10: Incomplete (30/200) | Total: {len(primes)} primes')
    
    return []

# ==================== Render ====================
total_frames = MAX_N * FRAMES_PER_SHELL + ENDING_FRAMES
print(f"Rendering {total_frames} frames...")

anim = animation.FuncAnimation(
    fig, animate, frames=total_frames,
    interval=INTERVAL, blit=False, repeat=False
)

print(f"Saving to: {OUTPUT_FILE}")
writer = animation.PillowWriter(fps=FPS)
anim.save(OUTPUT_FILE, writer=writer)
plt.close()

print(f"Done! File size: {os.path.getsize(OUTPUT_FILE) / 1024:.1f} KB")
print(f"Note: Shell 10 shows {missing_shell10} vacancies as hollow red circles")
