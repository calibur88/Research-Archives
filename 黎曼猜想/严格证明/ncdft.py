#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NCDFT — Non-Commutative Discrete Fourier Transform
        Riemann Zero Detection via Peak Spectrum
====================================================

Theory alignment:
  • Smooth cutoff:    φ(x),  x = log(n)/log(N) ∈ [0,1]         (§1.1)
  • Phase kernel:     exp(−2πi·k·log n / log N)               (§4.2)
  • Zero signature:   |D_N(t)| attains local MAXIMA at γ_j     (§4.4)
  • Validity range:   γ_max < √N  (asymptotic L² bound)        (Thm 4.2, §6.3)
  • Peak refinement:  concave-down parabolic fit on |D|²       (a < 0)

References:
  Edwards, H. M. (2001). Riemann's Zeta Function. Dover.
  Titchmarsh, E. C. (1986). The Theory of the Riemann Zeta-Function. OUP.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import mpmath as mp
import warnings

warnings.filterwarnings("ignore")

# High-precision Riemann zeros
mp.mp.dps = 50
plt.rcParams["figure.dpi"] = 150


# ----------------------------------------------------------------------
# Arithmetic primitives
# ----------------------------------------------------------------------
def _sieve(limit: int) -> np.ndarray:
    s = np.ones(limit + 1, dtype=bool)
    s[0:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    return np.where(s)[0]


def _von_mangoldt(n_max: int) -> np.ndarray:
    Lambda = np.zeros(n_max + 1)
    for p in _sieve(n_max):
        power = p
        while power <= n_max:
            Lambda[power] = np.log(p)
            power *= p
    return Lambda


# ----------------------------------------------------------------------
# Smooth cutoff φ(x)  —  §1.1
# ----------------------------------------------------------------------
def _phi(x: np.ndarray) -> np.ndarray:
    """
    Document §1.1:
      supp φ ⊂ [0,1],  φ(x)=1 on [0,1/2],
      smooth cos² roll-off on [1/2,1].
    """
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


# ----------------------------------------------------------------------
# Riemann zeros
# ----------------------------------------------------------------------
def _riemann_zeros(start: int, end: int) -> np.ndarray:
    return np.array([float(mp.zetazero(i).imag) for i in range(start, end + 1)])


# ----------------------------------------------------------------------
# NCDFT core operator
# ----------------------------------------------------------------------
class NCDFT:
    def __init__(self, N: int):
        self.N = N
        # Frequency grid: t_k = 2πk / log N   (§1.2)
        self.t_k = 2.0 * np.pi * np.arange(N) / np.log(N)
        self.Lambda = _von_mangoldt(N)

    def spectrum(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns:
          D_k  — Dirichlet polynomial D_N(t_k)  (§1.1, §4.2)
          φ_n  — cutoff weights
          f_n  — input sequence Λ(n)/√n · φ_n
        """
        N = self.N
        logN = np.log(N)
        n_arr = np.arange(1, N)
        Lambda_n = self.Lambda[1:N]

        x_ratio = np.log(n_arr) / logN
        phi_vals = _phi(x_ratio)
        f_n = Lambda_n / np.sqrt(n_arr) * phi_vals

        log_n = np.log(n_arr)
        # Phase kernel: exp(−2πi·k·log n / log N)   (§4.2)
        phases = np.exp(-2.0j * np.pi * np.outer(np.arange(N), log_n) / logN)
        D_k = phases @ f_n

        return D_k, phi_vals, f_n


# ----------------------------------------------------------------------
# Zero detection
# ----------------------------------------------------------------------
def detect(
    N: int,
    start: int,
    end: int,
) -> dict:
    """
    Peak-based zero detection on the |D_N(t)| spectrum.

    Returns a dictionary with detection results, diagnostics,
    and the full spectral data for plotting.
    """
    zeros = _riemann_zeros(start, end)
    ncdf = NCDFT(N)
    D, phi_vals, f_n = ncdf.spectrum()

    half = N // 2
    t_pos = ncdf.t_k[:half]
    mag = np.abs(D)[:half]

    dt = t_pos[1] - t_pos[0]
    bound = np.pi / np.log(N)
    theory_limit = np.sqrt(N)

    # --- Peak detection -------------------------------------------------
    med = np.median(mag)
    prominence = 0.20 * med
    avg_sp = 2.0 * np.pi / np.log(np.median(zeros))
    distance = max(1, int(0.5 * avg_sp / dt))

    peak_idx, _ = find_peaks(mag, distance=distance, prominence=prominence)
    peak_pos = t_pos[peak_idx]

    # --- Match peaks to true zeros ---------------------------------------
    matched = 0
    misses = []
    used = set()

    for j, gz in enumerate(zeros):
        gap = np.inf
        if j > 0:
            gap = min(gap, abs(gz - zeros[j - 1]))
        if j < len(zeros) - 1:
            gap = min(gap, abs(zeros[j + 1] - gz))
        window = np.clip(gap * 0.45, 0.6, 2.5)

        dists = np.abs(peak_pos - gz)
        if np.any(dists < window):
            matched += 1
            used.add(np.argmin(dists))
        else:
            misses.append(gz)

    fp = [peak_pos[k] for k in range(len(peak_pos)) if k not in used]

    # --- Parabolic vertex refinement (concave-down: a < 0) ---------------
    refined = []
    for j, gz in enumerate(zeros):
        gap = np.inf
        if j > 0:
            gap = min(gap, abs(gz - zeros[j - 1]))
        if j < len(zeros) - 1:
            gap = min(gap, abs(zeros[j + 1] - gz))
        window = np.clip(gap * 0.45, 0.6, 2.5)

        dists = np.abs(peak_pos - gz)
        if not np.any(dists < window):
            refined.append((gz, np.nan, np.nan))
            continue

        k = np.argmin(dists)
        best = peak_pos[k]
        k_idx = np.argmin(np.abs(t_pos - best))

        k0 = max(1, min(k_idx, len(t_pos) - 2))
        xs = t_pos[k0 - 1 : k0 + 2]
        ys = mag[k0 - 1 : k0 + 2] ** 2
        try:
            a, b, c = np.polyfit(xs, ys, 2)
            vertex = -b / (2 * a) if abs(a) > 1e-15 else best
            det = vertex if (xs[0] <= vertex <= xs[-1] and a < 0) else best
        except Exception:
            det = best

        refined.append((gz, det, abs(det - gz)))

    return {
        "N": N,
        "zeros": zeros,
        "t_pos": t_pos,
        "mag": mag,
        "phase": np.unwrap(np.angle(D)[:half]),
        "peak_idx": peak_idx,
        "peak_pos": peak_pos,
        "phi_vals": phi_vals,
        "f_n": f_n,
        "n_arr": np.arange(1, N),
        "matched": matched,
        "misses": misses,
        "fp": fp,
        "refined": refined,
        "bound": bound,
        "dt": dt,
        "theory_limit": theory_limit,
    }


# ----------------------------------------------------------------------
# Visualization
# ----------------------------------------------------------------------
def plot(results: dict, filename: str = "ncdf_spectrum.png") -> None:
    zeros = results["zeros"]
    t = results["t_pos"]
    mag = results["mag"]
    phase = results["phase"]
    peaks = results["peak_idx"]
    peak_pos = results["peak_pos"]
    misses = results["misses"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax0, ax1, ax2, ax3 = axes.flatten()

    # --- Global |D| -------------------------------------------------------
    ax0.plot(t, mag, "b-", lw=0.8, alpha=0.7, label="$|D_N(t)|$")
    if len(peaks) > 0:
        ax0.scatter(
            t[peaks],
            mag[peaks],
            c="red",
            s=30,
            zorder=5,
            label=f"peaks ({len(peaks)})",
        )
    for gz in zeros:
        ax0.axvline(gz, color="green", ls="--", lw=0.8, alpha=0.5)
    ax0.set_title("Global spectrum")
    ax0.set_xlabel("$t$")
    ax0.set_ylabel("$|D_N(t)|$")
    ax0.legend(loc="upper right", fontsize=9)
    ax0.grid(True, alpha=0.3)

    # --- Phase ------------------------------------------------------------
    ax1.plot(t, phase, "g-", lw=0.6, alpha=0.8)
    for gz in zeros:
        ax1.axvline(gz, color="green", ls="--", lw=0.8, alpha=0.5)
    ax1.set_title("Phase (unwrapped)")
    ax1.set_xlabel("$t$")
    ax1.set_ylabel("$\\arg D_N(t)$")
    ax1.grid(True, alpha=0.3)

    # --- Zoom on first 8 zeros --------------------------------------------
    t0, t1 = zeros[0] - 3, zeros[min(7, len(zeros) - 1)] + 3
    m = (t > t0) & (t < t1)
    ax2.plot(t[m], mag[m], "b-", lw=1.5)
    if len(peaks) > 0:
        zm = (peak_pos > t0) & (peak_pos < t1)
        if np.any(zm):
            pk = np.where(zm)[0]
            ax2.scatter(peak_pos[pk], mag[peaks[pk]], c="red", s=50, zorder=5)
    for gz in zeros[:8]:
        ax2.axvline(gz, color="green", ls="--", lw=1.2, alpha=0.7)
    for gz in misses:
        if t0 < gz < t1:
            ax2.annotate(
                "MISS",
                xy=(gz, ax2.get_ylim()[1] * 0.9),
                color="purple",
                fontsize=9,
                ha="center",
                fontweight="bold",
            )
    ax2.set_title("Zoom: zeros #1…#8")
    ax2.set_xlabel("$t$")
    ax2.set_ylabel("$|D_N(t)|$")
    ax2.grid(True, alpha=0.3)

    # --- Weight distribution ----------------------------------------------
    n = results["n_arr"]
    phi = results["phi_vals"]
    f_n = results["f_n"]
    ax3_twin = ax3.twinx()
    (l1,) = ax3.semilogy(n, phi, "k-", lw=1.0, label="$\\varphi(n)$")
    (l2,) = ax3_twin.semilogy(n, np.abs(f_n), "c-", lw=0.8, alpha=0.7, label="$|f_n|$")
    ax3.set_xlabel("$n$")
    ax3.set_ylabel("$\\varphi(n)$", color="k")
    ax3_twin.set_ylabel("$|f_n| = \\Lambda(n)/\\sqrt{n}\\cdot\\varphi(n)$", color="c")
    ax3.set_title("Cutoff weights")
    ax3.grid(True, alpha=0.3)
    ax3.legend(
        [l1, l2], [l1.get_label(), l2.get_label()], loc="upper right", fontsize=9
    )

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {filename}")
    plt.show()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main() -> None:
    N = 2000
    start, end = 1, 8  # Theory-compliant: γ₈ ≈ 43.3 < √N ≈ 44.7

    res = detect(N, start, end)

    print("=" * 60)
    print("NCDFT — Peak-Spectrum Zero Detection")
    print("=" * 60)
    print(f"N = {N},  √N = {res['theory_limit']:.1f}")
    print(f"Grid step Δt = {res['dt']:.4f}  |  Error bound = {res['bound']:.4f}")
    print(f"Peaks found: {len(res['peak_idx'])}")
    print(
        f"Matched: {res['matched']}/{len(res['zeros'])}  |  Misses: {len(res['misses'])}  |  FP: {len(res['fp'])}"
    )

    print("\nDetailed results:")
    print(f"{'#':>3} {'True γ':>10} {'Detected':>10} {'Error':>8} {'Status':>6}")
    print("-" * 45)
    for idx, (gz, det, err) in enumerate(res["refined"], start):
        if np.isnan(err):
            print(f"{idx:3d} {gz:10.4f} {'—':>10} {'—':>8} {'✗':>6}")
            continue
        ok = "✓" if err < 2.0 * res["bound"] else "✗"
        print(f"{idx:3d} {gz:10.4f} {det:10.4f} {err:8.4f} {ok:>6}")

    plot(res, "ncdf_spectrum.png")


if __name__ == "__main__":
    main()
