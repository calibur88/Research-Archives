#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
异型流形假说 - 脉冲星计时残差奇点检测
====================================================
方法：
1. 高阶导数检验：计算频率偏移的二阶、三阶导数，检测不连续点或δ尖峰
2. 变分模态分解（VMD）：将残差分解为IMF，检查频谱离散性（谐波结构）
3. 相位锁定检验：跨星突变时刻同步性分析
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
from sklearn.neighbors import KernelDensity
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# 配置
# ============================================================
DATA_DIR = "/storage/emulated/0/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband"
WINDOW_SIZE = 50
STEP_SIZE = 25
SIGMA_THRESHOLD = 3.0
VMD_ALPHA = 2000
VMD_TAU = 0.1
VMD_K = 5
VMD_ITER = 500
MAX_PULSARS = 50

# ============================================================
# 数据读取
# ============================================================
def parse_tim_file(filepath):
    mjds, residuals = [], []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('C'):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                mjd = float(parts[2])
                res = float(parts[3])
                mjds.append(mjd)
                residuals.append(res)
            except ValueError:
                continue
    if len(mjds) < 50:
        return None, None
    sort_idx = np.argsort(mjds)
    return np.array(mjds)[sort_idx], np.array(residuals)[sort_idx]


def extract_frequency_offsets(mjd, residual, window_size=WINDOW_SIZE, step=STEP_SIZE):
    if len(mjd) < window_size:
        return None, None, None
    times, offsets = [], []
    for i in range(0, len(mjd) - window_size, step):
        win_mjd = mjd[i:i+window_size]
        win_res = residual[i:i+window_size]
        x = win_mjd - win_mjd[0]
        A = np.vstack([x, np.ones_like(x)]).T
        slope, _ = np.linalg.lstsq(A, win_res, rcond=None)[0]
        delta_f_nHz = slope / 86400.0 * 1e9
        times.append(np.mean(win_mjd))
        offsets.append(delta_f_nHz)
    if len(times) < 5:
        return None, None, None
    return np.array(times), np.array(offsets), None


def high_order_derivative_analysis(times, offsets, sigma=SIGMA_THRESHOLD):
    if len(times) < 5:
        return [], []
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    offsets = offsets[sort_idx]
    unique_mask = np.concatenate(([True], np.diff(times) > 1e-8))
    times = times[unique_mask]
    offsets = offsets[unique_mask]
    if len(times) < 5:
        return [], []
    t_uniform = np.linspace(times[0], times[-1], len(times)*2)
    try:
        spline = UnivariateSpline(times, offsets, s=0.01, k=3)
    except:
        spline = UnivariateSpline(times, offsets, s=0.1, k=3)
    y_uniform = spline(t_uniform)
    dt = t_uniform[1] - t_uniform[0]
    d2 = np.gradient(np.gradient(y_uniform, dt), dt)
    d3 = np.gradient(d2, dt)
    std_d2 = np.std(np.abs(d2))
    if std_d2 < 1e-12:
        return [], []
    dist = max(1, int(0.5 / dt))
    peaks, _ = find_peaks(np.abs(d2), height=sigma*std_d2, distance=dist)
    peak_times = t_uniform[peaks] if len(peaks) > 0 else []
    std_d3 = np.std(np.abs(d3))
    peak_times3 = []
    if std_d3 > 1e-12:
        peaks3, _ = find_peaks(np.abs(d3), height=sigma*std_d3, distance=dist)
        peak_times3 = t_uniform[peaks3] if len(peaks3) > 0 else []
    all_peaks = np.concatenate([np.array(peak_times), np.array(peak_times3)])
    if len(all_peaks) == 0:
        return [], []
    all_peaks = np.unique(np.round(all_peaks, 2))
    return all_peaks.tolist(), None


# ============================================================
# VMD（简化实现，避免外部依赖）
# ============================================================
def vmd(signal, alpha=2000, tau=0.1, K=5, tol=1e-7):
    N = len(signal)
    f = np.fft.fftfreq(N, 1.0)[:N//2]
    u = np.random.randn(K, N)
    omega = np.zeros(K)
    signal_hat = np.fft.fft(signal)
    for it in range(VMD_ITER):
        u_hat = np.zeros((K, N), dtype=complex)
        for i in range(K):
            u_hat[i, :] = np.fft.fft(u[i, :])
        for i in range(K):
            num, den = 0.0, 0.0
            for j in range(1, N//2):
                if f[j] > 0:
                    num += f[j] * np.abs(u_hat[i, j])**2
                    den += np.abs(u_hat[i, j])**2
            if den > 0:
                omega[i] = num / den
        u_new = np.zeros_like(u)
        for i in range(K):
            residual = np.zeros(N, dtype=complex)
            for k in range(K):
                if k != i:
                    residual += u_hat[k, :]
            for j in range(1, N//2):
                if f[j] > 0:
                    numerator = (1 + alpha * (f[j] - omega[i])**2)**-1
                    u_hat[i, j] = numerator * (signal_hat[j] - residual[j])
            u_new[i, :] = np.real(np.fft.ifft(u_hat[i, :]))
        diff = np.mean([np.mean((u_new[i] - u[i])**2) for i in range(K)])
        u = u_new
        if diff < tol:
            break
    return u, omega


def vmd_analysis(times, offsets, K=VMD_K):
    if len(offsets) < 10:
        return None
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    offsets = offsets[sort_idx]
    unique_mask = np.concatenate(([True], np.diff(times) > 1e-8))
    times = times[unique_mask]
    offsets = offsets[unique_mask]
    if len(times) < 5:
        return None
    t_uniform = np.linspace(times[0], times[-1], len(times))
    try:
        spline = UnivariateSpline(times, offsets, s=0.01, k=3)
    except:
        spline = UnivariateSpline(times, offsets, s=0.1, k=3)
    y_uniform = spline(t_uniform)
    try:
        u, omega = vmd(y_uniform, alpha=VMD_ALPHA, tau=VMD_TAU, K=K)
    except:
        return None
    N = len(y_uniform)
    dt = t_uniform[1] - t_uniform[0]
    freqs = np.fft.fftfreq(N, dt)[:N//2]
    powers = []
    for i in range(K):
        fft_vals = np.fft.fft(u[i, :])[:N//2]
        powers.append(np.abs(fft_vals)**2)
    discrete_modes = []
    for i in range(K):
        p = powers[i]
        p_norm = p / (np.max(p) + 1e-12)
        dist = max(1, int(5 / (freqs[1]-freqs[0]))) if len(freqs) > 1 else 5
        peaks, _ = find_peaks(p_norm, height=0.2, distance=dist)
        if len(peaks) > 0 and p_norm[peaks[0]] > 0.3:
            discrete_modes.append(i)
    harmonic_detected = len(discrete_modes) >= 2
    return {
        'harmonic_detected': harmonic_detected,
        'discrete_modes': discrete_modes,
        'u': u,
        'freqs': freqs,
        'powers': powers
    }


def cross_pulsar_synchronization(mutation_times_list, time_window=10.0):
    all_times = []
    for times in mutation_times_list:
        all_times.extend(times)
    if len(all_times) < 5:
        return None
    all_times = np.array(all_times).reshape(-1, 1)
    kde = KernelDensity(bandwidth=time_window, kernel='gaussian')
    kde.fit(all_times)
    t_grid = np.linspace(np.min(all_times)-20, np.max(all_times)+20, 1000).reshape(-1, 1)
    log_dens = kde.score_samples(t_grid)
    dens = np.exp(log_dens)
    peaks, _ = find_peaks(dens, height=0.1*np.max(dens), distance=10)
    cluster_times = t_grid[peaks].flatten() if len(peaks) > 0 else []
    clusters_with_pulsars = []
    for ct in cluster_times:
        participating = []
        for i, times in enumerate(mutation_times_list):
            if np.any(np.abs(np.array(times) - ct) < time_window):
                participating.append(i)
        if participating:
            clusters_with_pulsars.append((ct, participating))
    if not clusters_with_pulsars:
        return None
    return {
        'cluster_times': cluster_times,
        'clusters_info': clusters_with_pulsars,
        'n_clusters': len(cluster_times)
    }


def main():
    print("=" * 70)
    print("异型流形假说 - 奇点检测（高阶导数 + VMD + 跨星同步）")
    print("=" * 70)
    print(f"数据目录: {DATA_DIR}\n")
    if not os.path.isdir(DATA_DIR):
        print(f"错误：目录不存在 -> {DATA_DIR}")
        sys.exit(1)
    tim_files = []
    for root, _, files in os.walk(DATA_DIR):
        for f in files:
            if f.endswith('.tim') and 'excise' not in f:
                tim_files.append(os.path.join(root, f))
    if len(tim_files) == 0:
        for root, _, files in os.walk(DATA_DIR):
            for f in files:
                if f.endswith('.tim'):
                    tim_files.append(os.path.join(root, f))
    print(f"找到 {len(tim_files)} 个 .tim 文件\n")
    if len(tim_files) == 0:
        sys.exit(1)
    if MAX_PULSARS > 0 and len(tim_files) > MAX_PULSARS:
        tim_files = tim_files[:MAX_PULSARS]
        print(f"限制处理前 {MAX_PULSARS} 颗脉冲星\n")
    all_mutation_times = []
    vmd_results = []
    psr_names = []
    for idx, filepath in enumerate(tim_files):
        psr_name = os.path.basename(filepath).replace('.tim', '')
        print(f"[{idx+1}/{len(tim_files)}] 处理: {psr_name}...", end='')
        mjd, residual = parse_tim_file(filepath)
        if mjd is None:
            print(" 跳过（数据不足）")
            continue
        times, offsets, _ = extract_frequency_offsets(mjd, residual)
        if times is None:
            print(" 跳过（无法提取频率）")
            continue
        peaks, _ = high_order_derivative_analysis(times, offsets)
        n_peaks = len(peaks)
        print(f" 突变点={n_peaks}", end='')
        all_mutation_times.append(peaks)
        vmd_res = None
        if len(offsets) > 20:
            vmd_res = vmd_analysis(times, offsets)
        harmonic = vmd_res['harmonic_detected'] if vmd_res else False
        print(f" VMD谐波={'Y' if harmonic else 'N'}", end='')
        if vmd_res:
            vmd_results.append(vmd_res)
        psr_names.append(psr_name)
        print()
    print(f"\n成功分析 {len(psr_names)} 颗脉冲星")
    print("\n" + "=" * 70)
    print("跨星同步检验结果")
    print("=" * 70)
    sync_res = None
    non_empty = [lst for lst in all_mutation_times if len(lst) > 0]
    if len(non_empty) < 2:
        print("只有一颗或没有脉冲星有突变，无法进行跨星同步检验")
    else:
        sync_res = cross_pulsar_synchronization(non_empty, time_window=10.0)
        if sync_res and sync_res['n_clusters'] > 0:
            print(f"发现 {sync_res['n_clusters']} 个同步簇")
            for ct, participants in sync_res['clusters_info']:
                print(f"  簇在 MJD {ct:.1f}: 包含 {len(participants)} 颗脉冲星")
        else:
            print("未发现显著同步簇")
    if vmd_results:
        harmonic_count = sum(1 for r in vmd_results if r['harmonic_detected'])
        print(f"\nVMD谐波检测: {harmonic_count}/{len(vmd_results)} 颗脉冲星显示谐波结构")
    print("\n" + "=" * 70)
    print("判决一（奇点检测版本）")
    print("=" * 70)
    total_mutations = sum(len(t) for t in all_mutation_times)
    avg_mutations = total_mutations / len(all_mutation_times) if all_mutation_times else 0
    has_sync = sync_res and sync_res['n_clusters'] > 0 if sync_res else False
    harmonic_ratio = harmonic_count/len(vmd_results)*100 if vmd_results else 0.0
    print(f"平均每颗脉冲星突变点数: {avg_mutations:.2f}")
    print(f"跨星同步簇: {'是' if has_sync else '否'}")
    print(f"VMD谐波比例: {harmonic_ratio:.1f}% (基于 {len(vmd_results)} 颗星)")
    if avg_mutations > 1.0 and has_sync:
        print("✅ 支持流型脉动：存在显著突变且跨星同步")
    elif avg_mutations > 1.0 and not has_sync:
        print("⚠️ 部分支持：存在突变但无跨星同步，可能来自局域噪声")
    elif harmonic_ratio > 30.0:
        print("⚠️ 部分支持：VMD检测到谐波结构，但突变证据不足")
    else:
        print("❌ 不支持流型脉动：突变稀少且未同步，无VMD谐波证据")
    print("\n详细结果已显示，如需进一步分析请调整参数重新运行。")


if __name__ == "__main__":
    main()