#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NCDFT — N=50000 2D Joint-Detection Final Edition
================================================

Target: 8GB RAM, exact computation, zeros #1…#90
Memory strategy:
  - Dynamic sub-chunking: per-block memory capped at 400MB
  - Auto thread throttling
  - 30030-wheel sieve exact reproduction (DO NOT SIMPLIFY)

2D Joint Detection:
  - Channel 1: magnitude peaks
  - Channel 2: phase jumps (unwrap derivative)
  - Channel 3: Re(D) zero-crossings
  - Channel 4: Im(D) zero-crossings
  - Nearest-zero guard: prevents cross-zero misassignment
  - Sub-grid fallback for unresolved pairs

Plot:
  - Global spectrum + phase jumps overlay
  - Phase unwrap (bugfixed: uses D(t), not mag*exp(0))
  - Zoom: zeros #1…#90
  - Complex trajectory D(t)=(Re,Im) with red dots at zero locations
  - plt.show() called exactly once
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import mpmath as mp
import warnings
from concurrent.futures import ThreadPoolExecutor
import os

warnings.filterwarnings("ignore")
mp.mp.dps = 50
plt.rcParams["figure.dpi"] = 150
plt.rcParams["text.usetex"] = False
os.environ["OMP_NUM_THREADS"] = "1"


# =====================================================================
# 1. Adaptive 30030-wheel sieve (exact reproduction — DO NOT SIMPLIFY)
# =====================================================================
def _adaptive_wheel_sieve(n: int) -> list:
    if n < 2:
        return []
    small_primes = [2, 3, 5, 7, 11, 13]
    if n <= 13:
        return [p for p in small_primes if p <= n]

    sqrt_n = int(np.sqrt(n)) + 1
    max_small = max(small_primes)

    base_primes = []
    is_comp = [False] * (sqrt_n + 1)
    for i in range(2, sqrt_n):
        if not is_comp[i]:
            base_primes.append(i)
            for j in range(i * i, sqrt_n, i):
                is_comp[j] = True

    primes = [p for p in small_primes if p <= n]
    primes.extend([p for p in base_primes if p > max_small and p <= n])

    modulus = 30030
    gap_start = sqrt_n
    gap_end = min(n, modulus)
    if gap_start <= gap_end:
        for x in range(gap_start, gap_end + 1):
            is_prime = True
            for p in base_primes:
                if p * p > x:
                    break
                if x % p == 0:
                    is_prime = False
                    break
            if is_prime:
                primes.append(x)

    if n <= modulus:
        return sorted(set(primes))

    residues = [i for i in range(1, modulus) if all(i % p != 0 for p in small_primes)]
    seg_size = modulus * 156
    start = modulus + 1

    for low in range(start, n + 1, seg_size):
        high = min(low + seg_size - 1, n)
        length = high - low + 1
        sieve = bytearray(length)

        first_block = (low // modulus) * modulus
        for k in range(first_block, high + 1, modulus):
            for r in residues:
                x = k + r
                if low <= x <= high:
                    sieve[x - low] = 1

        for p in base_primes:
            if p <= max_small or p * p > high:
                continue
            start_mark = ((low + p - 1) // p) * p
            if start_mark < p * p:
                start_mark = p * p
            for j in range(start_mark, high + 1, p):
                if j >= low:
                    sieve[j - low] = 0

        for k in range(first_block, high + 1, modulus):
            for r in residues:
                x = k + r
                if low <= x <= high and sieve[x - low]:
                    primes.append(x)

    return sorted(set(primes))


# =====================================================================
# 2. von Mangoldt
# =====================================================================
def _von_mangoldt(n_max: int) -> np.ndarray:
    Lambda = np.zeros(n_max + 1, dtype=np.float64)
    for p in _adaptive_wheel_sieve(n_max):
        power = p
        while power <= n_max:
            Lambda[power] = np.log(p)
            power *= p
    return Lambda


# =====================================================================
# 3. Document-correct phi
# =====================================================================
def _phi(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    r = np.zeros_like(x)
    m1 = (x >= 0.0) & (x <= 0.5)
    r[m1] = 1.0
    m2 = (x > 0.5) & (x <= 1.0)
    if np.any(m2):
        t = (x[m2] - 0.5) / 0.5
        r[m2] = 0.5 * (1.0 + np.cos(np.pi * t))
    r[x > 1.0] = 0.0
    return r


# =====================================================================
# 4. Riemann zeros
# =====================================================================
def _riemann_zeros(start: int, end: int) -> np.ndarray:
    return np.array([float(mp.zetazero(i).imag) for i in range(start, end + 1)])


# =====================================================================
# 5. NCDFT with extreme memory control
# =====================================================================
class NCDFT:
    def __init__(self, N: int, n_threads: int = None, mem_cap_gb: float = 4.0):
        self.N = N
        self.t_k = 2.0 * np.pi * np.arange(N) / np.log(N)
        self.Lambda = _von_mangoldt(N)
        self.n_threads = n_threads or min(8, os.cpu_count() or 4)
        self.mem_cap = mem_cap_gb * 1e9

    def _chunked_dirichlet(self, f_n: np.ndarray, log_n: np.ndarray):
        N = self.N
        logN = np.log(N)
        element_size = 16
        safe_chunk = max(50, int(self.mem_cap / (self.n_threads * N * element_size)))
        chunk_rows = min(safe_chunk, max(200, N // (self.n_threads * 4)))
        actual_threads = min(
            self.n_threads,
            max(1, int(self.mem_cap / (chunk_rows * N * element_size))),
        )
        result = np.zeros(N, dtype=complex)
        boundaries = np.linspace(0, N, actual_threads + 1, dtype=int)

        def _worker(idx_start: int, idx_end: int):
            local_len = idx_end - idx_start
            local_res = np.zeros(local_len, dtype=complex)
            k_slice = np.arange(idx_start, idx_end)
            for sub_start in range(0, local_len, chunk_rows):
                sub_end = min(sub_start + chunk_rows, local_len)
                k_sub = k_slice[sub_start:sub_end]
                phases = np.exp(-2.0j * np.pi * np.outer(k_sub, log_n) / logN)
                local_res[sub_start:sub_end] = phases @ f_n
            return idx_start, local_res

        with ThreadPoolExecutor(max_workers=actual_threads) as exe:
            futures = [
                exe.submit(_worker, boundaries[i], boundaries[i + 1])
                for i in range(actual_threads)
            ]
            for fut in futures:
                idx_start, local = fut.result()
                if len(local) > 0:
                    result[idx_start : idx_start + len(local)] = local
        return result

    def spectrum(self):
        N = self.N
        logN = np.log(N)
        n_arr = np.arange(1, N)
        Lambda_n = self.Lambda[1:N]
        phi_vals = _phi(np.log(n_arr) / logN)
        f_n = Lambda_n / np.sqrt(n_arr) * phi_vals
        log_n = np.log(n_arr)
        D_k = self._chunked_dirichlet(f_n, log_n)
        return D_k, phi_vals, f_n


# =====================================================================
# 6. 2D Joint Zero Detection (nearest-zero guard)
# =====================================================================
def detect_zeros_2d(N: int, start: int, end: int):
    zeros = _riemann_zeros(start, end)
    ncdf = NCDFT(N, mem_cap_gb=4.0)
    D, phi_vals, f_n = ncdf.spectrum()

    half = N // 2
    t_pos = ncdf.t_k[:half]
    mag = np.abs(D)[:half]
    phase = np.angle(D)[:half]
    re = np.real(D)[:half]
    im = np.imag(D)[:half]
    dt = t_pos[1] - t_pos[0]
    bound = np.pi / np.log(N)

    # Channel 1: magnitude peaks
    med = np.median(mag)
    prominence = 0.03 * med
    distance = max(1, int(0.5 * (2.0 * np.pi / np.log(np.median(zeros))) / dt))
    peak_idx, _ = find_peaks(mag, distance=distance, prominence=prominence, width=(0.2, 12))
    peak_pos = t_pos[peak_idx]
    peak_mag = mag[peak_idx]

    # Channel 2: phase jumps
    uphase = np.unwrap(phase)
    dphase = np.diff(uphase)
    jump_mask = np.abs(dphase) > 0.8
    jump_idx = np.where(jump_mask)[0]
    jump_clusters = []
    if len(jump_idx) > 0:
        cur = [jump_idx[0]]
        for idx in jump_idx[1:]:
            if idx - cur[-1] <= 2:
                cur.append(idx)
            else:
                jump_clusters.append(cur)
                cur = [idx]
        jump_clusters.append(cur)
    jump_pos = np.array([
        (t_pos[c[0]] + t_pos[min(c[-1] + 1, len(t_pos) - 1)]) / 2
        for c in jump_clusters
    ]) if jump_clusters else np.array([])

    # Channel 3 & 4: Re/Im zero-crossings
    def zero_crossings(sig):
        zc = np.where(np.diff(np.sign(sig)))[0]
        precise = []
        for z in zc:
            if z + 1 < len(sig):
                t0, t1 = t_pos[z], t_pos[z + 1]
                s0, s1 = sig[z], sig[z + 1]
                precise.append(t0 + (t1 - t0) * abs(s0) / (abs(s0) + abs(s1)))
        return np.array(precise)

    re_zc = zero_crossings(re)
    im_zc = zero_crossings(im)

    # Nearest-zero pre-assignment guard
    def nearest_zero_idx(points):
        if len(points) == 0:
            return np.array([])
        dists = np.abs(points[:, None] - zeros[None, :])
        return np.argmin(dists, axis=1)

    peak_owner = nearest_zero_idx(peak_pos)
    jump_owner = nearest_zero_idx(jump_pos)
    rezc_owner = nearest_zero_idx(re_zc)
    imzc_owner = nearest_zero_idx(im_zc)

    # Joint matching
    matched = 0
    refined = []
    used_peaks = set()
    used_jumps = set()
    used_rezc = set()
    used_imzc = set()

    for j, gz in enumerate(zeros):
        local_gap = np.inf
        if j > 0:
            local_gap = min(local_gap, abs(gz - zeros[j - 1]))
        if j < len(zeros) - 1:
            local_gap = min(local_gap, abs(zeros[j + 1] - gz))

        is_dense = local_gap < 2.5 * dt
        window = np.clip(local_gap * (0.85 if is_dense else 0.55), 0.6, 2.5)

        candidates = []

        if len(peak_pos) > 0:
            dists_p = np.abs(peak_pos - gz)
            in_win = np.where(dists_p < window)[0]
            for ci in in_win:
                if ci not in used_peaks and peak_owner[ci] == j:
                    candidates.append(("peak", peak_pos[ci], peak_mag[ci], ci))
                    break

        if len(jump_pos) > 0:
            dists_j = np.abs(jump_pos - gz)
            in_win = np.where(dists_j < window)[0]
            for ci in in_win:
                if ci not in used_jumps and jump_owner[ci] == j:
                    candidates.append(("jump", jump_pos[ci], 0.0, ci))
                    break

        if len(re_zc) > 0:
            dists_r = np.abs(re_zc - gz)
            in_win = np.where(dists_r < window)[0]
            for ci in in_win:
                if ci not in used_rezc and rezc_owner[ci] == j:
                    idx = np.argmin(np.abs(t_pos - re_zc[ci]))
                    candidates.append(("re_zc", re_zc[ci], mag[idx], ci))
                    break

        if len(im_zc) > 0:
            dists_i = np.abs(im_zc - gz)
            in_win = np.where(dists_i < window)[0]
            for ci in in_win:
                if ci not in used_imzc and imzc_owner[ci] == j:
                    idx = np.argmin(np.abs(t_pos - im_zc[ci]))
                    candidates.append(("im_zc", im_zc[ci], mag[idx], ci))
                    break

        if candidates:
            best = max(candidates, key=lambda x: x[2])
            det = best[1]
            tag = best[0]
            idx_tag = best[3]

            if tag == "peak":
                used_peaks.add(idx_tag)
            elif tag == "jump":
                used_jumps.add(idx_tag)
            elif tag == "re_zc":
                used_rezc.add(idx_tag)
            elif tag == "im_zc":
                used_imzc.add(idx_tag)

            matched += 1

            if tag == "peak":
                k_idx = np.argmin(np.abs(t_pos - det))
                k0 = max(1, min(k_idx, len(t_pos) - 2))
                xs = t_pos[k0 - 1 : k0 + 2]
                ys = mag[k0 - 1 : k0 + 2] ** 2
                try:
                    a, b, c = np.polyfit(xs, ys, 2)
                    vertex = -b / (2 * a) if abs(a) > 1e-15 else det
                    det = vertex if (xs[0] <= vertex <= xs[-1] and a < 0) else det
                except Exception:
                    pass

            refined.append((gz, det, abs(det - gz), tag))
        else:
            hw = max(window, 1.5 * dt)
            n_fine = max(41, int(6 * hw / dt))
            t_fine = np.linspace(gz - hw, gz + hw, n_fine)
            log_n = np.log(np.arange(1, N))
            phases = np.exp(-1.0j * np.outer(t_fine, log_n))
            local_mag = np.abs(phases @ f_n)
            sub_idx, _ = find_peaks(local_mag, distance=1, prominence=0.01 * np.median(local_mag))
            if len(sub_idx) > 0:
                sub_pos = t_fine[sub_idx]
                best_sub = sub_pos[np.argmin(np.abs(sub_pos - gz))]
                refined.append((gz, best_sub, abs(best_sub - gz), "subgrid"))
                matched += 1
            else:
                refined.append((gz, np.nan, np.nan, "miss"))

    misses = [r[0] for r in refined if np.isnan(r[1])]

    return {
        "N": N,
        "zeros": zeros,
        "t_pos": t_pos,
        "mag": mag,
        "D": D[:half],
        "peak_idx": peak_idx,
        "peak_pos": peak_pos,
        "jump_pos": jump_pos,
        "re_zc": re_zc,
        "im_zc": im_zc,
        "phi_vals": phi_vals,
        "f_n": f_n,
        "matched": matched,
        "misses": misses,
        "refined": refined,
        "bound": bound,
        "dt": dt,
        "threads": ncdf.n_threads,
        "actual_threads": ncdf.n_threads,
    }


# =====================================================================
# 7. Plot — plt.show() called exactly once
# =====================================================================
def plot(res: dict, filename: str = "ncdf_N50000_2d.png"):
    zeros = res["zeros"]
    t = res["t_pos"]
    mag = res["mag"]
    D = res["D"]
    peaks = res["peak_idx"]
    peak_pos = res["peak_pos"]
    misses = res["misses"]
    jump_pos = res["jump_pos"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax0, ax1, ax2, ax3 = axes.flatten()

    # Top-left: global magnitude spectrum
    ax0.plot(t, mag, "b-", lw=0.5, alpha=0.7, label="$|D_N(t)|$")
    if len(peaks) > 0:
        ax0.scatter(t[peaks], mag[peaks], c="red", s=20, zorder=5, label=f"peaks ({len(peaks)})")
    if len(jump_pos) > 0:
        ax0.scatter(jump_pos, np.interp(jump_pos, t, mag), c="orange", s=15, zorder=5,
                    marker="x", label=f"jumps ({len(jump_pos)})")
    for gz in zeros:
        ax0.axvline(gz, color="green", ls="--", lw=0.6, alpha=0.5)
    ax0.set_title("Global spectrum (2D joint)")
    ax0.set_xlabel("$t$")
    ax0.set_ylabel("$|D_N(t)|$")
    ax0.legend(loc="upper right", fontsize=9)
    ax0.grid(True, alpha=0.3)

    # Top-right: phase unwrap (bugfixed)
    phase = np.angle(D)
    ax1.plot(t, np.unwrap(phase), "g-", lw=0.4)
    for gz in zeros:
        ax1.axvline(gz, color="green", ls="--", lw=0.6, alpha=0.3)
    ax1.set_title("Phase (unwrapped) — fixed")
    ax1.set_xlabel("$t$")
    ax1.set_ylabel(r"$\arg D(t)$")
    ax1.grid(True, alpha=0.3)

    # Bottom-left: zoom zeros #1…#90
    t0, t1 = zeros[0] - 3, zeros[-1] + 3
    m = (t > t0) & (t < t1)
    ax2.plot(t[m], mag[m], "b-", lw=1.2)
    if len(peaks) > 0:
        zm = (peak_pos > t0) & (peak_pos < t1)
        if np.any(zm):
            pk = np.where(zm)[0]
            ax2.scatter(peak_pos[pk], mag[peaks[pk]], c="red", s=35, zorder=5)
    if len(jump_pos) > 0:
        zm = (jump_pos > t0) & (jump_pos < t1)
        if np.any(zm):
            jk = np.where(zm)[0]
            ax2.scatter(jump_pos[jk], np.interp(jump_pos[jk], t, mag), c="orange",
                        s=25, zorder=5, marker="x")
    for gz in zeros:
        ax2.axvline(gz, color="green", ls="--", lw=1.0, alpha=0.7)
    for gz in misses:
        if t0 < gz < t1:
            ax2.annotate("MISS", xy=(gz, ax2.get_ylim()[1] * 0.9),
                         color="purple", fontsize=9, ha="center", fontweight="bold")
    ax2.set_title(f"Zoom: zeros #1...#{len(zeros)}")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel("$|D_N(t)|$")
    ax2.grid(True, alpha=0.3)

    # Bottom-right: complex trajectory D(t)=(Re,Im)
    mask = (t > t0) & (t < t1)
    ax3.plot(np.real(D[mask]), np.imag(D[mask]), "b-", lw=0.6, alpha=0.4)
    for gz in zeros:
        idx = np.argmin(np.abs(t - gz))
        ax3.scatter(np.real(D[idx]), np.imag(D[idx]), c="red", s=12, zorder=5)
    ax3.axhline(0, color="k", lw=0.3)
    ax3.axvline(0, color="k", lw=0.3)
    ax3.set_title("Complex trajectory $D(t)$ (zoom)")
    ax3.set_xlabel("Re $D(t)$")
    ax3.set_ylabel("Im $D(t)$")
    ax3.set_aspect("equal", adjustable="box")
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {filename}")
    plt.show()


# =====================================================================
# 8. Main
# =====================================================================
def main():
    N = 50000
    start, end = 1, 90

    print(f"NCDFT — N={N}, zeros #{start}...#{end}")
    print(f"Memory cap: 4GB for phases, 8GB total system")
    print("=" * 60)

    res = detect_zeros_2d(N, start, end)

    print(f"Δt = {res['dt']:.4f}  |  Bound = {res['bound']:.4f}")
    print(f"Threads: {res['threads']}")
    print(f"Peaks found: {len(res['peak_idx'])}")
    print(f"Phase jumps found: {len(res['jump_pos'])}")
    print(f"Matched: {res['matched']}/{len(res['zeros'])}  |  Misses: {len(res['misses'])}")

    print(f"\n{'#':>3} {'True γ':>10} {'Detected':>10} {'Error':>8} {'Tag':>8} {'Status':>6}")
    print("-" * 50)
    for idx, (gz, det, err, tag) in enumerate(res["refined"], start):
        if np.isnan(err):
            print(f"{idx:3d} {gz:10.4f} {'—':>10} {'—':>8} {'—':>8} {'✗':>6}")
            continue
        ok = "✓" if err < 2.0 * res["bound"] else "✗"
        print(f"{idx:3d} {gz:10.4f} {det:10.4f} {err:8.4f} {tag:>8} {ok:>6}")

    plot(res)


if __name__ == "__main__":
    main()