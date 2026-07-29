#! python
# -*- coding:utf-8 -*-
###
# --------------------------------------------------------------------------------
# 文件名: 脉冲星计时残差奇点检测.py
# 创建时间: 2026-06-25 15:07:35 Thu
# 说明: 
# 作者: Calibur88
# 主机: LAPTOP-CGR9BDFR
# --------------------------------------------------------------------------------
# 最后编辑作者: Calibur88
# 最后修改时间: 2026-06-25 15:26:13 Thu
# --------------------------------------------------------------------------------
# Copyright (c) 2026 Calibur88
# --------------------------------------------------------------------------------
# 更新历史:
# --------------------------------------------------------------------------------
# 时间      		作者		信息
# ----------		---		------------------------------------------------------
###

"""
脉冲星计时残差奇点检测 — 完整分析流水线
============================================================
功能：
  1. 坐标提取（自动缓存）
  2. 全量数据分析：突变点、VMD谐波、跨星同步、天球各向异性
  3. 随机区间采样（前50颗星，多次随机子区间，检验稳定性）
  4. 随机批次采样（从全样本中随机抽取多批，检验复现性）
  5. 生成独立汇总图（每种方法一张）和综合汇总图
  6. 终端打印详细统计信息
"""

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
from sklearn.neighbors import KernelDensity
from scipy.stats import binomtest
import warnings
warnings.filterwarnings("ignore")

# ============================================================
#  配置（请根据实际路径修改）
# ============================================================
DATA_DIR = "./NANOGrav15yr/narrowband"      # 数据根目录（含 tim/ 和 par/）
OUTPUT_DIR = "./output_summary"             # 输出目录
COORD_JSON = "psr_coords.json"              # 坐标缓存文件
OUTPUT_PREFIX = "summary"                   # 输出文件名前缀
OUTPUT_JSON = "psr_analysis_results.json"   # 详细结果

# 分析参数
WINDOW_SIZE = 50
STEP_SIZE = 25
SIGMA_THRESHOLD = 3.0
VMD_ALPHA = 2000
VMD_TAU = 0.1
VMD_K = 5
VMD_ITER = 500
TIME_WINDOW = 10.0      # 同步簇时间窗口（天）

# 随机采样参数
RANDOM_INTERVAL_SAMPLES = 30    # 每颗星随机区间采样次数
RANDOM_INTERVAL_NSTARS = 50     # 参与随机区间采样的星数（前N颗）
BATCH_SIZE = 30                 # 每批脉冲星数量
N_BATCHES = 5                   # 批次数
SEED = 42

np.random.seed(SEED)

# ============================================================
#  中文字体设置（Windows 兼容）
# ============================================================
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'WenQuanYi Micro Hei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================
#  工具函数：坐标提取（从 .par 文件）
# ============================================================
OBLIQUITY_DEG = 23.4392911
OBLIQUITY_RAD = np.radians(OBLIQUITY_DEG)

def parse_par_value(par_path, key):
    try:
        with open(par_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith(key):
                    parts = line.split()
                    if len(parts) >= 2:
                        return parts[1].strip('"\'')
    except:
        pass
    return None

def hms_to_deg(hms_str):
    try:
        parts = hms_str.replace(':', ' ').split()
        if len(parts) == 3:
            h, m, s = float(parts[0]), float(parts[1]), float(parts[2])
            return (h + m/60 + s/3600) * 15.0
    except:
        return None

def dms_to_deg(dms_str):
    try:
        sign = 1.0
        clean = dms_str.strip()
        if clean.startswith('-'):
            sign = -1.0
            clean = clean[1:]
        parts = clean.replace(':', ' ').split()
        if len(parts) == 3:
            d, m, s = float(parts[0]), float(parts[1]), float(parts[2])
            return sign * (d + m/60 + s/3600)
    except:
        return None

def ecliptic_to_equatorial(elong_deg, elat_deg):
    lam = np.radians(elong_deg)
    beta = np.radians(elat_deg)
    eps = OBLIQUITY_RAD
    sin_delta = np.sin(beta)*np.cos(eps) + np.cos(beta)*np.sin(eps)*np.sin(lam)
    delta = np.arcsin(np.clip(sin_delta, -1.0, 1.0))
    y = np.sin(lam)*np.cos(eps) - np.tan(beta)*np.sin(eps)
    x = np.cos(lam)
    alpha = np.arctan2(y, x)
    ra_deg = np.degrees(alpha) % 360.0
    dec_deg = np.degrees(delta)
    return ra_deg, dec_deg

def get_pulsar_coords(par_path):
    raj = parse_par_value(par_path, 'RAJ')
    decj = parse_par_value(par_path, 'DECJ')
    if raj and decj:
        ra = hms_to_deg(raj)
        dec = dms_to_deg(decj)
        if ra is not None and dec is not None:
            return ra, dec, "RAJ/DECJ"
    elong = parse_par_value(par_path, 'ELONG')
    elat = parse_par_value(par_path, 'ELAT')
    if elong and elat:
        try:
            ra, dec = ecliptic_to_equatorial(float(elong), float(elat))
            return ra, dec, "ELONG/ELAT"
        except:
            pass
    lam = parse_par_value(par_path, 'LAMBDA')
    beta = parse_par_value(par_path, 'BETA')
    if lam and beta:
        try:
            ra, dec = ecliptic_to_equatorial(float(lam), float(beta))
            return ra, dec, "LAMBDA/BETA"
        except:
            pass
    return None, None, "FAIL"

def extract_coordinates(data_root, coord_json):
    print("\n🔭 正在提取脉冲星坐标...")
    par_index = {}
    for root, _, files in os.walk(data_root):
        for f in files:
            if f.endswith('.par'):
                full = os.path.join(root, f)
                key = f[:-4]
                par_index[key] = full
    coord_map = {}
    for key, path in par_index.items():
        ra, dec, src = get_pulsar_coords(path)
        if ra is not None:
            coord_map[key] = {'RA': float(ra), 'DEC': float(dec), 'source': src}
    with open(coord_json, 'w') as f:
        json.dump(coord_map, f, indent=2)
    print(f"✅ 提取了 {len(coord_map)} 颗脉冲星坐标 -> {coord_json}")
    return coord_map

# ============================================================
#  数据读取与频率偏移提取
# ============================================================
def parse_tim(filepath):
    mjds, resids = [], []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('C'):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                mjds.append(float(parts[2]))
                resids.append(float(parts[3]))
            except:
                continue
    if len(mjds) < 50:
        return None, None
    idx = np.argsort(mjds)
    return np.array(mjds)[idx], np.array(resids)[idx]

def extract_freq(mjd, resid, ws=WINDOW_SIZE, step=STEP_SIZE):
    if len(mjd) < ws:
        return None, None
    times, offsets = [], []
    for i in range(0, len(mjd) - ws, step):
        wm, wr = mjd[i:i+ws], resid[i:i+ws]
        x = wm - wm[0]
        A = np.vstack([x, np.ones_like(x)]).T
        slope, _ = np.linalg.lstsq(A, wr, rcond=None)[0]
        times.append(np.mean(wm))
        offsets.append(slope / 86400.0 * 1e9)
    if len(times) < 5:
        return None, None
    return np.array(times), np.array(offsets)

# ============================================================
#  高阶导数突变检测
# ============================================================
def deriv_analysis(times, offsets, sigma=SIGMA_THRESHOLD):
    if len(times) < 5:
        return [], None
    idx = np.argsort(times)
    times, offsets = times[idx], offsets[idx]
    mask = np.concatenate(([True], np.diff(times) > 1e-8))
    times, offsets = times[mask], offsets[mask]
    if len(times) < 5:
        return [], None
    tu = np.linspace(times[0], times[-1], len(times)*2)
    try:
        spline = UnivariateSpline(times, offsets, s=0.01, k=3)
    except:
        spline = UnivariateSpline(times, offsets, s=0.1, k=3)
    yu = spline(tu)
    dt = tu[1] - tu[0]
    d2 = np.gradient(np.gradient(yu, dt), dt)
    d3 = np.gradient(d2, dt)
    std2, std3 = np.std(np.abs(d2)), np.std(np.abs(d3))
    dist = max(1, int(0.5/dt))
    peaks = []
    if std2 > 1e-12:
        p, _ = find_peaks(np.abs(d2), height=sigma*std2, distance=dist)
        peaks.extend(tu[p].tolist())
    if std3 > 1e-12:
        p3, _ = find_peaks(np.abs(d3), height=sigma*std3, distance=dist)
        peaks.extend(tu[p3].tolist())
    return np.unique(np.round(peaks, 2)).tolist(), None

# ============================================================
#  VMD 分析（含镜像延拓）
# ============================================================
def vmd(signal, alpha=VMD_ALPHA, tau=VMD_TAU, K=VMD_K, tol=1e-7, max_iter=VMD_ITER):
    N = len(signal)
    if N < 10:
        return None, None, None
    f_m = np.concatenate([signal[::-1], signal, signal[::-1]])
    Nm = len(f_m)
    f_hat = np.fft.fft(f_m)
    freqs = np.fft.fftfreq(Nm, 1.0)
    omega = np.linspace(0, 0.5, K)
    u_hat = np.zeros((K, Nm), dtype=complex)
    for k in range(K):
        mask = (np.abs(freqs) >= k*0.5/K) & (np.abs(freqs) < (k+1)*0.5/K)
        u_hat[k, mask] = f_hat[mask] / (np.sum(mask) + 1e-12)
    lam = np.zeros(Nm, dtype=complex)
    for _ in range(max_iter):
        for k in range(K):
            sum_u = np.sum(u_hat, axis=0) - u_hat[k, :]
            denom = 1.0 + alpha * (freqs - omega[k])**2
            u_hat[k, :] = (f_hat - sum_u + lam/tau) / denom
        for k in range(K):
            num = np.sum(np.abs(freqs) * np.abs(u_hat[k, :])**2)
            den = np.sum(np.abs(u_hat[k, :])**2)
            if den > 1e-12:
                omega[k] = num / den
        sum_u = np.sum(u_hat, axis=0)
        lam += tau * (f_hat - sum_u)
        if np.sum(np.abs(f_hat - sum_u)**2) / np.sum(np.abs(f_hat)**2) < tol:
            break
    u = np.zeros((K, Nm))
    for k in range(K):
        u[k, :] = np.real(np.fft.ifft(u_hat[k, :]))
    return u[:, N:N+Nm//3], omega, u_hat

def vmd_spectrum(u, dt, K=VMD_K):
    N = u.shape[1]
    freqs = np.fft.fftfreq(N, dt)[:N//2]
    powers, pfreqs, ppows = [], [], []
    for k in range(K):
        fftv = np.fft.fft(u[k, :])[:N//2]
        power = np.abs(fftv)**2
        powers.append(power)
        if len(power) > 5:
            p, _ = find_peaks(power, height=0.05*np.max(power), distance=max(1, len(power)//20))
            pfreqs.append(freqs[p] if len(p) else np.array([]))
            ppows.append(power[p] if len(p) else np.array([]))
        else:
            pfreqs.append(np.array([]))
            ppows.append(np.array([]))
    return freqs, powers, pfreqs, ppows

def vmd_analysis(times, offsets, K=VMD_K):
    if len(offsets) < 20:
        return None
    idx = np.argsort(times)
    times, offsets = times[idx], offsets[idx]
    mask = np.concatenate(([True], np.diff(times) > 1e-8))
    times, offsets = times[mask], offsets[mask]
    if len(times) < 10:
        return None
    tu = np.linspace(times[0], times[-1], len(times))
    try:
        spline = UnivariateSpline(times, offsets, s=0.01, k=3)
    except:
        spline = UnivariateSpline(times, offsets, s=0.1, k=3)
    yu = spline(tu)
    dt = tu[1] - tu[0]
    u, omega, _ = vmd(yu, alpha=VMD_ALPHA, tau=VMD_TAU, K=K, max_iter=VMD_ITER)
    if u is None:
        return None
    freqs, powers, pfreqs, ppows = vmd_spectrum(u, dt, K)
    all_peaks = []
    for k in range(K):
        if len(pfreqs[k]) > 0:
            all_peaks.extend(pfreqs[k].tolist())
    harmonic = False
    if len(all_peaks) >= 3:
        all_peaks = np.sort(np.unique(np.round(all_peaks, 4)))
        if len(all_peaks) >= 3:
            diffs = np.diff(all_peaks)
            if len(diffs) >= 2:
                cv = np.std(diffs) / (np.mean(diffs) + 1e-12)
                if cv < 0.5:
                    harmonic = True
    discrete = sum(1 for k in range(K) if len(pfreqs[k]) > 0 and np.max(powers[k]) > 0.1)
    if discrete >= 2:
        harmonic = True
    return {'harmonic': harmonic, 'discrete': discrete}

# ============================================================
#  跨星同步检验
# ============================================================
def sync_analysis(mutation_times_list, tw=TIME_WINDOW):
    all_t = []
    for t in mutation_times_list:
        all_t.extend(t)
    if len(all_t) < 5:
        return None
    all_t = np.array(all_t).reshape(-1, 1)
    kde = KernelDensity(bandwidth=tw, kernel='gaussian')
    kde.fit(all_t)
    t_grid = np.linspace(np.min(all_t)-20, np.max(all_t)+20, 1000).reshape(-1, 1)
    dens = np.exp(kde.score_samples(t_grid))
    peaks, _ = find_peaks(dens, height=0.1*np.max(dens), distance=10)
    ctimes = t_grid[peaks].flatten() if len(peaks) else []
    clusters = []
    for ct in ctimes:
        part = []
        for i, t in enumerate(mutation_times_list):
            if len(t) > 0 and np.any(np.abs(np.array(t) - ct) < tw):
                part.append(i)
        if part:
            clusters.append((ct, part))
    if not clusters:
        return None
    return {
        'ctimes': ctimes,
        'clusters': clusters,
        'n': len(ctimes),
        'dens': dens,
        't_grid': t_grid.flatten()
    }

# ============================================================
#  天球各向异性检验
# ============================================================
def rayleigh_test(ra_deg, dec_deg):
    ra_rad = np.radians(ra_deg)
    dec_rad = np.radians(dec_deg)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    R_vec = np.array([np.sum(x), np.sum(y), np.sum(z)])
    R = np.linalg.norm(R_vec) / len(ra_deg)
    N = len(ra_deg)
    p_value = np.exp(-N * R**2)
    return R, p_value

# ============================================================
#  随机区间采样（蒙特卡洛稳定性测试）
# ============================================================
def random_interval_test(mjd, resid, n_samples=RANDOM_INTERVAL_SAMPLES):
    """对单颗星数据随机抽取子区间，计算突变点数均值"""
    if len(mjd) < 200:
        return None
    counts = []
    for _ in range(n_samples):
        # 随机起始点和长度
        total = len(mjd)
        if total < 100:
            continue
        max_start = total - 100
        if max_start < 0:
            continue
        start = np.random.randint(0, max_start)
        length = np.random.randint(80, min(150, total - start))
        if length < 50:
            continue
        t, o = extract_freq(mjd[start:start+length], resid[start:start+length], ws=30, step=15)
        if t is not None and len(t) > 5:
            peaks, _ = deriv_analysis(t, o, sigma=SIGMA_THRESHOLD)
            counts.append(len(peaks))
        else:
            counts.append(0)
    if not counts:
        return None
    return {'mean': np.mean(counts), 'std': np.std(counts), 'n': len(counts)}

# ============================================================
#  随机批次采样
# ============================================================
def batch_analysis(tim_files, coord_map, batch_size=BATCH_SIZE, n_batches=N_BATCHES):
    """从所有脉冲星中随机抽取多批，每批分析统计"""
    print("\n🔀 开始随机批次采样...")
    # 筛选可用的星（有坐标且数据有效）
    valid_stars = []
    for fp in tim_files:
        psr = os.path.basename(fp).replace('.tim', '')
        if psr not in coord_map:
            continue
        mjd, resid = parse_tim(fp)
        if mjd is None:
            continue
        valid_stars.append(fp)
    if len(valid_stars) < batch_size:
        print(f"⚠️  有效星数 ({len(valid_stars)}) 少于批次大小 {batch_size}，减少批次大小")
        batch_size = len(valid_stars) // 2
        if batch_size < 5:
            return None

    batch_results = []
    for b in range(n_batches):
        # 随机抽取
        batch = np.random.choice(valid_stars, size=batch_size, replace=False).tolist()
        n_muts_batch = []
        vmd_flags_batch = []
        for fp in batch:
            psr = os.path.basename(fp).replace('.tim', '')
            mjd, resid = parse_tim(fp)
            if mjd is None:
                continue
            t, o = extract_freq(mjd, resid)
            if t is None:
                continue
            peaks, _ = deriv_analysis(t, o)
            n_muts_batch.append(len(peaks))
            vres = vmd_analysis(t, o) if len(o) > 20 else None
            harmonic = vres['harmonic'] if vres else False
            vmd_flags_batch.append(harmonic)
        if n_muts_batch:
            batch_results.append({
                'batch_id': b+1,
                'n_stars': len(n_muts_batch),
                'avg_mut': np.mean(n_muts_batch),
                'vmd_ratio': sum(vmd_flags_batch)/len(vmd_flags_batch)*100 if vmd_flags_batch else 0
            })
        else:
            batch_results.append({'batch_id': b+1, 'n_stars': 0, 'avg_mut': 0, 'vmd_ratio': 0})
    print(f"✅ 完成 {n_batches} 批采样")
    return batch_results

# ============================================================
#  可视化函数
# ============================================================
def plot_mutation_histogram(n_muts, save_path):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.hist(n_muts, bins=15, color='#3498db', edgecolor='black', alpha=0.8)
    ax.axvline(np.mean(n_muts), color='red', linestyle='--', linewidth=2,
               label=f'均值 = {np.mean(n_muts):.2f}')
    ax.set_xlabel('突变点数量', fontsize=12)
    ax.set_ylabel('脉冲星数', fontsize=12)
    ax.set_title('突变点数量分布', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 突变分布图已保存: {save_path}")

def plot_vmd_pie(vmd_flags, save_path):
    harmonic_count = sum(vmd_flags)
    non_harmonic = len(vmd_flags) - harmonic_count
    fig, ax = plt.subplots(figsize=(9, 7))
    wedges, texts, autotexts = ax.pie(
        [harmonic_count, non_harmonic],
        labels=['检测到谐波', '未检测到'],
        colors=['#2ecc71', '#e74c3c'],
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 12},
        wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}
    )
    for a in autotexts:
        a.set_fontsize(14)
        a.set_fontweight('bold')
    ax.set_title(f'VMD谐波检测 (占比 {harmonic_count/len(vmd_flags)*100:.1f}%)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ VMD饼图已保存: {save_path}")

def plot_sync_density(sync_res, save_path):
    fig, ax = plt.subplots(figsize=(10, 7))
    if sync_res:
        ax.fill_between(sync_res['t_grid'], sync_res['dens'], alpha=0.4, color='#3498db')
        ax.plot(sync_res['t_grid'], sync_res['dens'], color='#2980b9', linewidth=2)
        for ct, part in sync_res['clusters']:
            if len(part) >= 3:
                ax.axvline(x=ct, color='#e74c3c', alpha=0.6, linewidth=2)
                ax.annotate(f'{len(part)}星', xy=(ct, max(sync_res['dens'])*0.9),
                            fontsize=9, color='#c0392b', fontweight='bold')
        ax.set_xlabel('MJD', fontsize=12)
        ax.set_ylabel('突变事件密度', fontsize=12)
        ax.set_title(f'跨星同步簇 (共 {sync_res["n"]} 个簇)', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, '无显著同步簇', ha='center', va='center', fontsize=14, transform=ax.transAxes)
        ax.set_title('跨星同步检验', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 同步密度图已保存: {save_path}")

def plot_sky_map(ras, decs, n_muts, rayleigh_p, save_path):
    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw={'projection': 'aitoff'})
    ras_arr = np.array(ras)
    decs_arr = np.array(decs)
    lon = np.radians(ras_arr - 180)
    lat = np.radians(decs_arr)
    colors = []
    sizes = []
    for n in n_muts:
        if n >= 10:
            colors.append('#e74c3c')
            sizes.append(60)
        elif n >= 5:
            colors.append('#f39c12')
            sizes.append(40)
        else:
            colors.append('#3498db')
            sizes.append(20)
    sc = ax.scatter(lon, lat, c=colors, s=sizes, alpha=0.8, edgecolors='black', linewidth=0.5)
    ax.set_title(f'天球分布 (Rayleigh p = {rayleigh_p:.2e})', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3)
    legend_elements = [
        Patch(facecolor='#e74c3c', edgecolor='black', label='突变点 ≥10'),
        Patch(facecolor='#f39c12', edgecolor='black', label='突变点 5-9'),
        Patch(facecolor='#3498db', edgecolor='black', label='突变点 <5')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 天球分布图已保存: {save_path}")

def plot_random_interval(interval_results, save_path):
    """绘制随机区间采样结果（前50颗星的均值分布）"""
    valid = [r for r in interval_results if r is not None]
    if not valid:
        print("⚠️  无有效的随机区间结果，跳过绘图")
        return
    means = [r['mean'] for r in valid]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.hist(means, bins=15, color='#9b59b6', edgecolor='black', alpha=0.8)
    ax.axvline(np.mean(means), color='red', linestyle='--', linewidth=2,
               label=f'全局均值 = {np.mean(means):.2f}')
    ax.set_xlabel('随机区间突变点均值', fontsize=12)
    ax.set_ylabel('脉冲星数', fontsize=12)
    ax.set_title(f'随机区间采样稳定性 (前{RANDOM_INTERVAL_NSTARS}颗星)', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 随机区间结果图已保存: {save_path}")

def plot_batch_results(batch_results, save_path):
    """绘制批次采样结果"""
    if not batch_results:
        return
    batch_ids = [f"B{r['batch_id']}" for r in batch_results]
    avg_muts = [r['avg_mut'] for r in batch_results]
    vmd_ratios = [r['vmd_ratio'] for r in batch_results]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    ax1.bar(batch_ids, avg_muts, color='#e67e22', edgecolor='black', alpha=0.85)
    ax1.axhline(np.mean(avg_muts), color='red', ls='--', lw=2, label=f'均值 {np.mean(avg_muts):.2f}')
    ax1.set_xlabel('批次')
    ax1.set_ylabel('平均突变点数')
    ax1.set_title('各批次平均突变点')
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax2.bar(batch_ids, vmd_ratios, color='#1abc9c', edgecolor='black', alpha=0.85)
    ax2.axhline(np.mean(vmd_ratios), color='red', ls='--', lw=2, label=f'均值 {np.mean(vmd_ratios):.1f}%')
    ax2.set_xlabel('批次')
    ax2.set_ylabel('VMD谐波检测率 (%)')
    ax2.set_title('各批次VMD谐波率')
    ax2.legend()
    ax2.grid(alpha=0.3)
    plt.suptitle('随机批次采样结果', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 批次结果图已保存: {save_path}")

def plot_combined(psr_names, n_muts, vmd_flags, mut_times_list, ras, decs, sync_res, rayleigh_p, save_path):
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(n_muts, bins=15, color='#3498db', edgecolor='black', alpha=0.8)
    ax1.axvline(np.mean(n_muts), color='red', linestyle='--', linewidth=2, label=f'均值 = {np.mean(n_muts):.2f}')
    ax1.set_xlabel('突变点数量')
    ax1.set_ylabel('脉冲星数')
    ax1.set_title('突变点数量分布')
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax2 = fig.add_subplot(gs[0, 1])
    harmonic_count = sum(vmd_flags)
    non_harmonic = len(vmd_flags) - harmonic_count
    wedges, texts, autotexts = ax2.pie(
        [harmonic_count, non_harmonic],
        labels=['检测到谐波', '未检测到'],
        colors=['#2ecc71', '#e74c3c'],
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 12},
        wedgeprops={'edgecolor': 'black', 'linewidth': 1.5}
    )
    for a in autotexts:
        a.set_fontsize(13)
        a.set_fontweight('bold')
    ax2.set_title(f'VMD谐波检测 (占比 {harmonic_count/len(vmd_flags)*100:.1f}%)')
    ax3 = fig.add_subplot(gs[1, 0])
    if sync_res:
        ax3.fill_between(sync_res['t_grid'], sync_res['dens'], alpha=0.4, color='#3498db')
        ax3.plot(sync_res['t_grid'], sync_res['dens'], color='#2980b9', linewidth=2)
        for ct, part in sync_res['clusters']:
            if len(part) >= 3:
                ax3.axvline(x=ct, color='#e74c3c', alpha=0.6, linewidth=2)
                ax3.annotate(f'{len(part)}星', xy=(ct, max(sync_res['dens'])*0.9),
                             fontsize=9, color='#c0392b', fontweight='bold')
        ax3.set_xlabel('MJD')
        ax3.set_ylabel('突变事件密度')
        ax3.set_title(f'跨星同步簇 (共 {sync_res["n"]} 个簇)')
        ax3.grid(alpha=0.3)
    else:
        ax3.text(0.5, 0.5, '无显著同步簇', ha='center', va='center', fontsize=14, transform=ax3.transAxes)
        ax3.set_title('跨星同步检验')
    ax4 = fig.add_subplot(gs[1, 1], projection='aitoff')
    ras_arr = np.array(ras)
    decs_arr = np.array(decs)
    lon = np.radians(ras_arr - 180)
    lat = np.radians(decs_arr)
    colors = []
    sizes = []
    for n in n_muts:
        if n >= 10:
            colors.append('#e74c3c')
            sizes.append(60)
        elif n >= 5:
            colors.append('#f39c12')
            sizes.append(40)
        else:
            colors.append('#3498db')
            sizes.append(20)
    sc = ax4.scatter(lon, lat, c=colors, s=sizes, alpha=0.8, edgecolors='black', linewidth=0.5)
    ax4.set_title(f'天球分布 (Rayleigh p = {rayleigh_p:.2e})')
    ax4.grid(alpha=0.3)
    legend_elements = [
        Patch(facecolor='#e74c3c', edgecolor='black', label='突变点 ≥10'),
        Patch(facecolor='#f39c12', edgecolor='black', label='突变点 5-9'),
        Patch(facecolor='#3498db', edgecolor='black', label='突变点 <5')
    ]
    ax4.legend(handles=legend_elements, loc='upper left', fontsize=9)
    fig.suptitle('脉冲星计时残差奇点检测 — 汇总分析', fontsize=18, fontweight='bold', y=0.98)
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ 综合汇总图已保存: {save_path}")

# ============================================================
#  主程序
# ============================================================
def main():
    print("\n" + "="*70)
    print(" 脉冲星计时残差奇点检测 — 完整分析流水线（含随机验证）")
    print("="*70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- 步骤1：坐标提取 ----
    coord_map = {}
    if os.path.exists(COORD_JSON):
        with open(COORD_JSON, 'r') as f:
            coord_map = json.load(f)
        print(f"✅ 加载坐标缓存: {len(coord_map)} 颗")
    else:
        coord_map = extract_coordinates(DATA_DIR, COORD_JSON)
        if not coord_map:
            print("❌ 未提取到任何坐标，请检查数据目录结构（需包含 .par 文件）")
            sys.exit(1)

    # ---- 步骤2：扫描所有 .tim 文件 ----
    tim_files = []
    for root, _, files in os.walk(DATA_DIR):
        for f in files:
            if f.endswith('.tim') and 'excise' not in f:
                tim_files.append(os.path.join(root, f))

    print(f"\n📁 找到 {len(tim_files)} 个 .tim 文件")
    if not tim_files:
        print("❌ 未找到任何 .tim 文件")
        sys.exit(1)

    # ---- 步骤3：逐颗分析（全量） ----
    results = []
    all_mut_times = []
    ras, decs = [], []
    n_muts = []
    vmd_flags = []
    # 同时存储数据用于随机区间采样
    full_mjd_resid = []

    print("\n🔍 开始全量分析脉冲星...")
    for idx, fp in enumerate(tim_files):
        psr = os.path.basename(fp).replace('.tim', '')
        if psr not in coord_map:
            continue
        mjd, resid = parse_tim(fp)
        if mjd is None:
            continue
        # 存储原始数据用于随机区间（只存前N颗）
        if len(full_mjd_resid) < RANDOM_INTERVAL_NSTARS:
            full_mjd_resid.append((mjd, resid))
        times, offsets = extract_freq(mjd, resid)
        if times is None:
            continue
        peaks, _ = deriv_analysis(times, offsets)
        n_mut = len(peaks)
        vres = vmd_analysis(times, offsets) if len(offsets) > 20 else None
        harmonic = vres['harmonic'] if vres else False
        ra = coord_map[psr]['RA']
        dec = coord_map[psr]['DEC']
        results.append({
            'name': psr,
            'RA': ra,
            'DEC': dec,
            'n_mutations': n_mut,
            'vmd_harmonic': harmonic,
            'mutation_times': peaks
        })
        all_mut_times.append(peaks)
        ras.append(ra)
        decs.append(dec)
        n_muts.append(n_mut)
        vmd_flags.append(harmonic)

        if (idx+1) % 50 == 0 or idx == len(tim_files)-1:
            print(f"  进度: {idx+1}/{len(tim_files)} 已扫描, 有效星数: {len(results)}")

    print(f"✅ 成功分析 {len(results)} 颗脉冲星")
    if len(results) < 10:
        print("❌ 有效样本太少，无法进行有意义的统计。")
        sys.exit(1)

    # ---- 步骤4：跨星同步 ----
    sync_res = None
    non_empty = [lst for lst in all_mut_times if lst]
    if len(non_empty) >= 2:
        sync_res = sync_analysis(non_empty, TIME_WINDOW)
        if sync_res:
            print(f"✅ 发现 {sync_res['n']} 个同步簇")
        else:
            print("ℹ️  未发现显著同步簇")

    # ---- 步骤5：天球各向异性 ----
    R, p_val = rayleigh_test(np.array(ras), np.array(decs))
    print(f"\n📊 Rayleigh 均匀性检验: R = {R:.4f}, p = {p_val:.2e}")
    if p_val < 0.05:
        print("   ✅ 存在显著各向异性 (p < 0.05)")
    else:
        print("   ❌ 未检测到显著各向异性 (p ≥ 0.05)")

    # ---- 步骤6：随机区间采样（前50颗） ----
    print(f"\n🔄 开始随机区间采样（前{RANDOM_INTERVAL_NSTARS}颗星，每颗{RANDOM_INTERVAL_SAMPLES}次）...")
    interval_results = []
    for i, (mjd, resid) in enumerate(full_mjd_resid[:RANDOM_INTERVAL_NSTARS]):
        res = random_interval_test(mjd, resid, RANDOM_INTERVAL_SAMPLES)
        interval_results.append(res)
        if (i+1) % 10 == 0:
            print(f"  进度: {i+1}/{len(full_mjd_resid)}")
    print(f"✅ 随机区间采样完成")

    # ---- 步骤7：随机批次采样 ----
    batch_results = batch_analysis(tim_files, coord_map, BATCH_SIZE, N_BATCHES)
    if batch_results:
        print("  批次统计:")
        for br in batch_results:
            print(f"    Batch {br['batch_id']}: 星数={br['n_stars']}, 平均突变={br['avg_mut']:.2f}, VMD率={br['vmd_ratio']:.1f}%")

    # ---- 步骤8：保存详细结果 ----
    json_results = []
    for r in results:
        json_results.append({
            'name': r['name'],
            'RA': r['RA'],
            'DEC': r['DEC'],
            'n_mutations': r['n_mutations'],
            'vmd_harmonic': r['vmd_harmonic']
        })
    with open(os.path.join(OUTPUT_DIR, OUTPUT_JSON), 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"✅ 详细结果已保存: {os.path.join(OUTPUT_DIR, OUTPUT_JSON)}")

    # ---- 步骤9：生成独立汇总图 ----
    plot_mutation_histogram(n_muts, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_mutations.png"))
    plot_vmd_pie(vmd_flags, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_vmd.png"))
    plot_sync_density(sync_res, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_sync.png"))
    plot_sky_map(ras, decs, n_muts, p_val, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_sky.png"))
    # 随机区间结果图
    if interval_results:
        plot_random_interval(interval_results, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_random_interval.png"))
    # 批次结果图
    if batch_results:
        plot_batch_results(batch_results, os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_batch.png"))
    # 综合汇总图
    plot_combined(
        [r['name'] for r in results],
        n_muts,
        vmd_flags,
        all_mut_times,
        ras,
        decs,
        sync_res,
        p_val,
        os.path.join(OUTPUT_DIR, f"{OUTPUT_PREFIX}_combined.png")
    )

    # ---- 最终统计报告 ----
    print("\n" + "="*70)
    print(" 📈 最终统计摘要")
    print("="*70)
    print(f"  分析脉冲星总数: {len(results)}")
    print(f"  突变点总数: {sum(n_muts)}")
    print(f"  平均每星突变点: {np.mean(n_muts):.2f} ± {np.std(n_muts):.2f}")
    print(f"  VMD谐波检测率: {sum(vmd_flags)/len(vmd_flags)*100:.1f}%")
    if sync_res:
        print(f"  同步簇数量: {sync_res['n']}")
        if sync_res['clusters']:
            max_size = max(len(c[1]) for c in sync_res['clusters'])
            print(f"  最大同步簇包含: {max_size} 颗星")
    print(f"  Rayleigh p值: {p_val:.2e}")
    if interval_results:
        valid_means = [r['mean'] for r in interval_results if r is not None]
        if valid_means:
            print(f"  随机区间采样（前{RANDOM_INTERVAL_NSTARS}颗）均值: {np.mean(valid_means):.2f} ± {np.std(valid_means):.2f}")
    if batch_results:
        avg_batch_mut = np.mean([r['avg_mut'] for r in batch_results])
        avg_batch_vmd = np.mean([r['vmd_ratio'] for r in batch_results])
        print(f"  随机批次采样（{N_BATCHES}批）平均突变: {avg_batch_mut:.2f}, VMD率: {avg_batch_vmd:.1f}%")

    print(f"\n✅ 所有输出保存在: {OUTPUT_DIR}/")
    print("   - summary_mutations.png   : 突变点数量分布")
    print("   - summary_vmd.png         : VMD谐波检测饼图")
    print("   - summary_sync.png        : 跨星同步密度图")
    print("   - summary_sky.png         : 天球Aitoff投影图")
    print("   - summary_random_interval.png : 随机区间采样稳定性")
    print("   - summary_batch.png       : 随机批次采样结果")
    print("   - summary_combined.png    : 四合一综合汇总图")
    print("   - psr_analysis_results.json")
    print("="*70)

if __name__ == "__main__":
    main()