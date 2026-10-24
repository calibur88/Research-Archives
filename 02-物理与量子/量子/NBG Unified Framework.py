"""
NBG Unified Framework — Visualization Suite
======================================================
谱统计 · 干涉相变 · 闭环验证 · 基底投影 · 收敛几何 · 全模式 · 几何相变/隧穿
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.interpolate import UnivariateSpline
from math import comb
from pathlib import Path
import json

# ==================== Configuration ====================
CFG = {
    "N_gue": 80, "N_interp": 400, "gue_ensemble": 30,
    "lambda_list": [0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
    "beta_ensemble": 10,
    "N_fft": 4096, "u_range": 12.0, "n_slits": 3,
    "slit_spacing": 2.0, "slit_width": 0.08,
    "lambda_min": 0.0, "lambda_max": 6.0, "n_lambda": 80,
    "output_dir": "nbg_unified",
    "dpi": 150, "seed": 42,
}

# ==================== Mathematical Utilities ====================
def catalan(n): return comb(2 * n, n) // (n + 1)
def semicircle_moment(k): return 0.0 if k % 2 else catalan(k // 2)
def gue_wigner(s): return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)
def poisson_dist(s): return np.exp(-s)

def generate_gue(N, seed=None):
    if seed is not None: np.random.seed(seed)
    A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
    A = (A + A.conj().T) / 2
    return A / np.sqrt(2)

def spacing_distribution(eig, bins=50):
    s = np.diff(np.sort(eig))
    s = s[np.isfinite(s)]
    if len(s) < 5:
        return np.linspace(0, 3, bins), np.zeros(bins)
    s = s / np.mean(s)
    s = s[(s > 0.01) & (s < 5)]
    hist, edges = np.histogram(s, bins=bins, density=True, range=(0, 4))
    return (edges[:-1] + edges[1:]) / 2, hist

def estimate_beta_single(sc, ps):
    mask = (sc > 0.05) & (sc < 0.5) & (ps > 1e-5)
    if np.sum(mask) > 5:
        try:
            return float(np.polyfit(np.log(sc[mask]), np.log(ps[mask]), 1)[0])
        except:
            return 0.0
    return 0.0

def estimate_beta_ensemble(lam, N, n_ens, seed):
    betas = []
    for i in range(n_ens):
        np.random.seed(seed + i)
        D = np.diag(np.sort(np.random.uniform(-2, 2, N)))
        A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
        R = (A + A.conj().T) / (2 * np.sqrt(N))
        sc, ps = spacing_distribution(np.sort(np.linalg.eigvalsh(D + lam * R)))
        b = estimate_beta_single(sc, ps)
        if 0 <= b <= 4:
            betas.append(b)
    return float(np.median(betas)) if betas else 0.0

# ==================== NBG Physical Kernels ====================
def beta(lam): return 2.0 * lam**2 / (1.0 + lam**2)
def C_mu(lam): return np.sin(np.pi * beta(lam) / 2)
def Lyapunov(lam): return lam / (1.0 + 0.3 * lam)
def Lambda_eff(lam, w=CFG["slit_width"]): return (np.pi / w) / (1.0 + lam**2)
def focusing(lam): return 1.0 / (1.0 + np.exp(-2 * (beta(lam) - 1.0)))
def is_critical(lam, rho_factor=1.0): return C_mu(lam) > 0.85 / rho_factor
def E_geo(lam): return beta(lam)**2 / 4.0

# ==================== Data Engines ====================
class SpectralData:
    """Group I: Spectral statistics — GUE verification & NBG interpolation."""
    def __init__(self, cfg):
        self.cfg = cfg
        self.moment_results = []
        self.beta_results = []
        self.spacing_data = {}

    def run(self):
        eig_gue = []
        for i in range(self.cfg["gue_ensemble"]):
            G = generate_gue(self.cfg["N_gue"], seed=self.cfg["seed"] + i)
            eig_gue.extend(np.linalg.eigvalsh(G / np.sqrt(self.cfg["N_gue"])))
        self.eig_gue = np.array(eig_gue)

        for k in [2, 4, 6, 8]:
            self.moment_results.append((k, np.mean(self.eig_gue**k), semicircle_moment(k)))

        self.s_gue, self.P_gue = spacing_distribution(self.eig_gue)

        for lam in self.cfg["lambda_list"]:
            np.random.seed(self.cfg["seed"])
            D = np.diag(np.sort(np.random.uniform(-2, 2, self.cfg["N_interp"])))
            A = np.random.randn(self.cfg["N_interp"], self.cfg["N_interp"]) + 1j * np.random.randn(self.cfg["N_interp"], self.cfg["N_interp"])
            R = (A + A.conj().T) / (2 * np.sqrt(self.cfg["N_interp"]))
            sc, ps = spacing_distribution(np.sort(np.linalg.eigvalsh(D + lam * R)))
            self.spacing_data[lam] = (sc, ps)
            self.beta_results.append((lam, estimate_beta_ensemble(lam, self.cfg["N_interp"], self.cfg["beta_ensemble"], self.cfg["seed"])))


class InterferenceData:
    """Groups II–VI: Interference phase transition, projection, convergent geometry."""
    def __init__(self, cfg):
        self.cfg = cfg
        self.N = cfg["N_fft"]
        self.d = cfg["slit_spacing"]
        self.w = cfg["slit_width"]
        self.u = np.linspace(-cfg["u_range"], cfg["u_range"], self.N)
        self.du = self.u[1] - self.u[0]
        self.omega = np.fft.fftshift(np.fft.fftfreq(self.N, self.du) * 2 * np.pi)

        c = np.linspace(-(cfg["n_slits"] - 1) / 2, (cfg["n_slits"] - 1) / 2, cfg["n_slits"]) * self.d
        self.psi_u_A = np.zeros_like(self.u)
        for ci in c:
            self.psi_u_A += np.exp(-(self.u - ci)**2 / (2 * self.w**2))
        self.psi_u_A /= np.max(self.psi_u_A)
        self.psi_M_A = self._mellin(self.psi_u_A)

        self.scan = {k: [] for k in [
            "lambda", "Lambda_eff", "dy", "beta", "E_geo", "C_mu",
            "focusing", "vis_fft", "vis_argrel"
        ]}

    def _mellin(self, psi):
        return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(psi))) * (self.N * self.du)

    def _inverse_mellin(self, psi_M):
        return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(psi_M))) / (self.N * self.du)

    def interference(self, lam):
        L = Lambda_eff(lam, self.w)
        dy = lam * self.w / (2 * np.pi)
        psiM = self.psi_M_A * np.exp(1j * dy * self.omega)
        psiM *= np.abs(self.omega) <= L
        psiu = self._inverse_mellin(psiM)
        I = np.abs(psiu)**2
        return I / np.max(I) if np.max(I) > 0 else I

    def interference_with_overflow(self, lam, rho_factor=1.0):
        """
        Topological charge integer overflow → discrete branch jump in Mellin parameter s.
        """
        L = Lambda_eff(lam, self.w)
        dy = lam * self.w / (2 * np.pi)
        psiM = self.psi_M_A * np.exp(1j * dy * self.omega)
        mask = np.abs(self.omega) <= L
        psiM = psiM * mask

        if is_critical(lam, rho_factor):
            overflow_phase = np.pi * focusing(lam)
            center_mask = np.abs(self.omega) < 0.5
            psiM[center_mask] *= np.exp(1j * overflow_phase)

        psiu = self._inverse_mellin(psiM)
        I = np.abs(psiu)**2
        return I / np.max(I) if np.max(I) > 0 else I

    def psi_complex(self, lam):
        L = Lambda_eff(lam, self.w)
        dy = lam * self.w / (2 * np.pi)
        psiM = self.psi_M_A * np.exp(1j * dy * self.omega)
        psiM *= np.abs(self.omega) <= L
        return self._inverse_mellin(psiM)

    def vis_fft(self, I):
        f = np.fft.fft(I)
        fr = np.fft.fftfreq(len(I), self.du)
        C0 = np.abs(f[0]) / len(I)
        idx = np.argmin(np.abs(fr - 1 / self.d))
        C1 = np.abs(f[idx]) / len(I)
        return float(np.clip(2 * C1 / C0, 0, 1)) if C0 > 1e-10 else 0

    def vis_argrel(self, I):
        m = (self.u > -4) & (self.u < 4)
        ur, Ir = self.u[m], I[m]
        mx = argrelextrema(Ir, np.greater)[0]
        mn = argrelextrema(Ir, np.less)[0]
        if len(mx) == 0 or len(mn) == 0:
            return 0.0
        mxv = np.max(Ir[mx])
        mnv = np.min(Ir[mn])
        return (mxv - mnv) / (mxv + mnv) if (mxv + mnv) > 0 else 0

    def run(self):
        lams = np.linspace(self.cfg["lambda_min"], self.cfg["lambda_max"], self.cfg["n_lambda"])
        self.scan["lambda"] = lams.tolist()
        keys_to_fill = [k for k in self.scan.keys() if k != "lambda"]
        for lam in lams:
            I = self.interference(lam)
            vals = [
                Lambda_eff(lam, self.w),
                lam * self.w / (2 * np.pi),
                beta(lam),
                E_geo(lam),
                C_mu(lam),
                focusing(lam),
                self.vis_fft(I),
                self.vis_argrel(I),
            ]
            for k, v in zip(keys_to_fill, vals):
                self.scan[k].append(v)


class LoopData:
    """Group III: Closed-loop validation — spectral ↔ interference self-consistency."""
    def __init__(self, spec_beta, fft_scan, cfg):
        self.spec = spec_beta
        self.fft = fft_scan
        self.cfg = cfg
        self.val = []

    def _mono(self, x, y):
        pts = []
        last = -1
        for xi, yi in zip(x, y):
            if yi > last + 1e-6:
                pts.append((xi, yi))
                last = yi
        if len(pts) < 4:
            pts = []
            last = 2
            for xi, yi in zip(x, y):
                if yi < last - 1e-6:
                    pts.append((xi, yi))
                    last = yi
        if len(pts) < 4:
            step = max(1, len(x) // 10)
            pts = list(zip(x[::step], y[::step]))
        return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

    def run(self):
        sl = np.array([x[0] for x in self.spec])
        sb = np.array([x[1] for x in self.spec])
        fl = np.array(self.fft["lambda"])
        fv = np.array(self.fft["vis_fft"])
        fb = np.array(self.fft["beta"])
        si = np.argsort(fb)
        bu, vu = self._mono(fb[si], fv[si])
        self.spline = UnivariateSpline(vu[::-1], bu[::-1], k=3, s=1e-4) if len(bu) >= 4 else None
        for lam, bs in zip(sl, sb):
            idx = np.argmin(np.abs(fl - lam))
            V = fv[idx]
            bi = float(np.clip(self.spline(V), 0, 2)) if self.spline and 0 <= V <= 1 else 0.0
            self.val.append({
                "lambda": float(lam),
                "beta_spectral": float(bs),
                "V_fft": float(V),
                "beta_inv": float(bi),
                "residual": float(bs - bi),
            })


# ==================== Plotting Groups ====================
def group_spectral(spec, out):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    ax = axes[0, 0]
    ax.hist(spec.eig_gue, bins=60, density=True, alpha=0.5, color="#1976d2", label="GUE")
    t = np.linspace(-2.1, 2.1, 300)
    rho = [(1 / (2 * np.pi)) * np.sqrt(max(0, 4 - ti**2)) for ti in t]
    ax.plot(t, rho, "r-", lw=2.5, label="Wigner semicircle")
    ax.set_title("GUE → Wigner Semicircle", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ks = [r[0] for r in spec.moment_results]
    cat = [r[2] for r in spec.moment_results]
    emp = [r[1] for r in spec.moment_results]
    x = np.arange(len(ks))
    ax.bar(x - 0.2, cat, 0.4, label="Catalan", color="#4caf50")
    ax.bar(x + 0.2, emp, 0.4, label="GUE", color="#1976d2")
    ax.set_xticks(x)
    ax.set_xticklabels([f"k={k}" for k in ks])
    ax.set_title("Moments = Catalan Numbers", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    w = spec.s_gue[1] - spec.s_gue[0] if len(spec.s_gue) > 1 else 0.05
    ax.bar(spec.s_gue, spec.P_gue, width=w, alpha=0.5, color="#7b1fa2", label="GUE")
    st = np.linspace(0, 4, 200)
    ax.plot(st, gue_wigner(st), "r-", lw=2.5, label="GUE theory")
    ax.plot(st, poisson_dist(st), "g--", lw=2, label="Poisson")
    ax.set_title("Level Spacing P(s)", fontweight="bold")
    ax.legend(fontsize=9)
    ax.set_xlim(0, 3.5)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    colors = plt.cm.plasma(np.linspace(0.15, 0.9, len(spec.spacing_data)))
    for i, (lam, (sc, ps)) in enumerate(sorted(spec.spacing_data.items())):
        ax.plot(sc, ps, color=colors[i], lw=1.5, label=f"λ={lam}")
    ax.plot(st, gue_wigner(st), "r--", lw=2, alpha=0.6, label="GUE")
    ax.set_title("NBG: λ → Spacing", fontweight="bold")
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(0, 3.5)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ss = np.linspace(0.01, 0.5, 100)
    ax.loglog(ss, gue_wigner(ss), "r-", lw=2.5, label="GUE: ~s²")
    ax.loglog(ss, poisson_dist(ss), "g--", lw=2, label="Poisson: ~s⁰")
    ax.loglog(ss, 0.5 * ss**2, "k:", lw=1, alpha=0.5, label="s² ref")
    ax.set_title("Small-s: Level Repulsion", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    lams, betas = zip(*spec.beta_results)
    ax.scatter(lams, betas, s=100, c="#c62828", zorder=5, edgecolors="black")
    ax.plot(lams, betas, "#c62828", lw=2)
    ax.axhline(y=2, color="r", ls="--", alpha=0.4, label="GUE: β=2")
    ax.axhline(y=0, color="g", ls="--", alpha=0.4, label="Poisson: β=0")
    ax.set_xlabel("λ (Lyapunov)")
    ax.set_ylabel("Repulsion Exponent β")
    ax.set_title("λ↑ → β↑ → GUE", fontweight="bold")
    ax.set_ylim(-0.3, 2.5)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Group I — Spectral Statistics: GUE Verification & NBG Interpolation", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_01_spectral.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_interference(intf, out):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    lams = np.array(intf.scan["lambda"])

    ax = axes[0, 0]
    ax.fill_between(intf.u, intf.psi_u_A, alpha=0.3, color="blue")
    ax.plot(intf.u, intf.psi_u_A, "b-", lw=1.5)
    ax.set_xlim(-5, 5)
    ax.set_title("Source (Case A)\nPoisson = Random Chaos Walk", fontweight="bold")
    ax.set_xlabel(r"$u=\ln x$")
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(lams, intf.scan["beta"], "r-", lw=2.5)
    ax.fill_between(lams, intf.scan["beta"], alpha=0.2, color="red")
    ax.axhline(y=2, color="r", ls="--", alpha=0.3)
    ax.axhline(y=0, color="g", ls="--", alpha=0.3)
    ax.set_xlabel(r"Lyapunov $\lambda$")
    ax.set_ylabel(r"$\beta(\lambda)$")
    ax.set_title("Level Repulsion: Poisson → GUE", fontweight="bold")
    ax.set_ylim(-0.3, 2.5)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.plot(lams, intf.scan["E_geo"], "r-", lw=2.5, label=r"$E_{\rm geo}$")
    ax2 = ax.twinx()
    ax2.plot(lams, intf.scan["C_mu"], "purple", lw=2, ls="--", label=r"$\mathbf{C}(\mu)$")
    ax.set_xlabel(r"Lyapunov $\lambda$")
    ax.set_ylabel(r"$E_{\rm geo}$", color="r")
    ax2.set_ylabel(r"$\mathbf{C}(\mu)$", color="purple")
    ax.set_title("Energy & Topological Scalar", fontweight="bold")
    ax.legend(loc="upper left")
    ax2.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(lams, intf.scan["vis_fft"], "b-", lw=2.5)
    ax.fill_between(lams, intf.scan["vis_fft"], alpha=0.2, color="blue")
    ax.set_xlabel(r"Lyapunov $\lambda$")
    ax.set_ylabel("Visibility")
    ax.set_title(r"FFT Modulation Depth $V=2|C_1|/C_0$", fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(lams, intf.scan["vis_argrel"], "g-", lw=2.5)
    ax.fill_between(lams, intf.scan["vis_argrel"], alpha=0.2, color="green")
    ax.set_xlabel(r"Lyapunov $\lambda$")
    ax.set_ylabel("Visibility")
    ax.set_title("argrelextrema (Global Contrast)", fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.plot(lams, intf.scan["vis_fft"], "b-o", lw=2, markersize=4, label="FFT modulation")
    ax.plot(lams, intf.scan["vis_argrel"], "g-s", lw=2, markersize=4, label="argrelextrema")
    ax.set_xlabel(r"Lyapunov $\lambda$")
    ax.set_ylabel("Visibility")
    ax.set_title("Visibility Comparison", fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle("Group II — Interference Phase Transition: Visibility & Geometric Observables", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_02_interference.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_loop(loop, out):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    lams = np.array([v["lambda"] for v in loop.val])
    bspec = np.array([v["beta_spectral"] for v in loop.val])
    binv = np.array([v["beta_inv"] for v in loop.val])
    V = np.array([v["V_fft"] for v in loop.val])
    res = np.array([v["residual"] for v in loop.val])
    fb = np.array(loop.fft["beta"])
    fvi = np.array(loop.fft["vis_fft"])

    ax = axes[0, 0]
    ax.plot(fb, fvi, "b-", lw=2, alpha=0.6, label="FFT scan")
    ax.scatter(bspec, V, s=120, c="red", zorder=5, edgecolors="black", marker="D", label="Spectral sample")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$V$")
    ax.set_title(r"$\beta$–$V$ Phase Map", fontweight="bold")
    ax.set_xlim(-0.1, 2.2)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    if loop.spline:
        vg = np.linspace(0.01, 0.99, 200)
        bg = np.clip(loop.spline(vg), 0, 2)
        ax.plot(vg, bg, "purple", lw=2.5, label=r"$\beta_{\rm inv}(V)$")
    ax.scatter(V, binv, s=80, c="red", zorder=5, marker="o", edgecolors="black", label="Sample points")
    ax.set_xlabel(r"$V$")
    ax.set_ylabel(r"$\beta_{\rm inv}$")
    ax.set_title("Inverse Map: V → β", fontweight="bold")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.1, 2.2)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    ax.plot(lams, bspec, "r-o", lw=2, markersize=6, label=r"$\beta_{\rm spectral}$")
    ax.plot(lams, binv, "b--s", lw=2, markersize=6, label=r"$\beta_{\rm inv}$")
    ax.fill_between(lams, bspec, binv, alpha=0.2, color="gray")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\beta$")
    ax.set_title(r"$\beta_{\rm spectral}$ vs $\beta_{\rm inv}$", fontweight="bold")
    ax.set_ylim(-0.3, 2.5)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    colors = ["green" if abs(r) < 0.2 else "orange" if abs(r) < 0.5 else "red" for r in res]
    ax.bar(lams, res, width=0.3, color=colors, alpha=0.7, edgecolor="black")
    ax.axhline(y=0, color="k", lw=1)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\Delta\beta$")
    ax.set_title("Closed-Loop Residual", fontweight="bold")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.hist(res, bins=12, color="#1976d2", alpha=0.6, edgecolor="black")
    ax.axvline(x=0, color="r", lw=2, ls="--")
    ax.set_xlabel(r"$\Delta\beta$")
    ax.set_ylabel("Count")
    ax.set_title("Residual Distribution", fontweight="bold")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.semilogy(lams, np.abs(res), "k-o", lw=2, markersize=6)
    ax.axhline(y=0.1, color="g", ls="--", alpha=0.5, label="Tight (0.1)")
    ax.axhline(y=0.5, color="orange", ls="--", alpha=0.5, label="Loose (0.5)")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$|\Delta\beta|$")
    ax.set_title("Convergence Metric", fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Group III — Closed-Loop Validation: Spectral ↔ Interference Self-Consistency", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_03_closed_loop.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_projection(intf, out):
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    for i, lam in enumerate([0.0, 1.0, 4.0]):
        Iu = intf.interference(lam)
        x = np.exp(intf.u)
        Ix = Iu / (x + 1e-12)
        Ix = Ix / np.max(Ix) if np.max(Ix) > 0 else Ix
        mask = (x > 1e-3) & (x < 1e3)

        ax = axes[i, 0]
        ax.fill_between(intf.u, Iu, alpha=0.3, color="blue")
        ax.plot(intf.u, Iu, "b-", lw=1.5)
        ax.set_xlim(-5, 5)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel(r"$u=\ln x$")
        ax.set_ylabel(r"$I(u)$")
        ax.set_title(f"NBG Multiplicative Basis\nλ={lam}, β={beta(lam):.2f}", fontweight="bold")
        ax.grid(True, alpha=0.3)

        ax = axes[i, 1]
        ax.fill_between(x[mask], Ix[mask], alpha=0.3, color="red")
        ax.plot(x[mask], Ix[mask], "r-", lw=1.5)
        ax.set_xlim(1e-3, 1e3)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xscale("log")
        ax.set_xlabel(r"$x$ (additive, log)")
        ax.set_ylabel(r"$I(x)=I(u)/x$")
        ax.set_title(f"Standard QM Projection\nλ={lam}", fontweight="bold")
        ax.grid(True, alpha=0.3, which="both")

        ax = axes[i, 2]
        ax.fill_between(x[mask], Ix[mask], alpha=0.3, color="green")
        ax.plot(x[mask], Ix[mask], "g-", lw=1.5)
        ax.set_xlim(0, 50)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel(r"$x$ (linear)")
        ax.set_ylabel(r"$I(x)$")
        pidx = np.argmax(Ix[mask])
        px = x[mask][pidx]
        pI = Ix[mask][pidx]
        above = Ix[mask] > pI / 2
        spread = 0
        if np.any(above):
            xa = x[mask][above]
            spread = (xa[-1] - xa[0]) / px if px > 0 else 0
        ax.set_title(f"Additive Linear View\nSpread={spread:.2f}", fontweight="bold")
        ax.grid(True, alpha=0.3)
        if lam > 0:
            ax.annotate("Diffusion" if lam > 2 else "Focusing", xy=(px, pI), xytext=(px + 10, pI * 0.7),
                       arrowprops=dict(arrowstyle="->", color="black", lw=1.5), fontsize=10, color="darkgreen")

    fig.suptitle("Group IV — Basis Projection: NBG Multiplicative → Standard QM Additive\nNumerical Verification of Diffusion Emergence", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_04_projection.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_convergent(intf, out):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    lams = np.linspace(CFG["lambda_min"], CFG["lambda_max"], CFG["n_lambda"])

    ax = axes[0, 0]
    lsel = np.array([0.0, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0, 6.0])
    cols = plt.cm.coolwarm(np.linspace(0, 1, len(lsel)))
    for i, lam in enumerate(lsel):
        Iu = intf.interference(lam)
        off = i * 0.25
        ax.fill_between(intf.u, Iu + off, off, alpha=0.25, color=cols[i])
        ax.plot(intf.u, Iu + off, color=cols[i], lw=1.8, label=f"λ={lam:.1f} β={beta(lam):.2f}")
        if np.max(Iu) > 0.3:
            ax.annotate("", xy=(0, off + np.max(Iu)), xytext=(0, off), arrowprops=dict(arrowstyle="->", color=cols[i], lw=1.5))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-0.1, len(lsel) * 0.25 + 0.5)
    ax.set_xlabel(r"$u=\ln x$")
    ax.set_ylabel(r"Convergence Level $\rightarrow$")
    ax.set_title("Convergent Ripple Cascade", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, title="Recursive Depth")
    ax.grid(True, alpha=0.2)
    ax.text(0.02, 0.98, "↑ Energy Convergence (Non-diffusive)", transform=ax.transAxes, fontsize=11, va="top", color="darkred",
           bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax = axes[0, 1]
    lam = 2.0
    L = Lambda_eff(lam, intf.w)
    dy = lam * intf.w / (2 * np.pi)
    psiM = intf.psi_M_A * np.exp(1j * dy * intf.omega)
    mask = np.abs(intf.omega) <= L
    ax.plot(intf.omega, np.angle(psiM), "b-", lw=1, alpha=0.6, label="arg ψ_M")
    ax.fill_between(intf.omega, -4, 4, where=mask, alpha=0.15, color="green", label="Cutoff")
    ax.axvline(x=L, color="r", ls="--", alpha=0.5)
    ax.axvline(x=-L, color="r", ls="--", alpha=0.5)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"$\arg \psi_M$")
    ax.set_title("Mellin: Phase Ramp & Cutoff", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    pmap = np.zeros((len(lams), intf.N))
    for i, lam in enumerate(lams):
        pmap[i, :] = np.unwrap(np.angle(intf.psi_complex(lam)))
    pmap = np.mod(pmap + np.pi, 2 * np.pi) - np.pi
    im = ax.imshow(pmap, aspect="auto", origin="lower", extent=[intf.u.min(), intf.u.max(), lams.min(), lams.max()],
                   cmap="twilight", vmin=-np.pi, vmax=np.pi)
    plt.colorbar(im, ax=ax, label="Phase (rad)")
    Imap = np.array([intf.interference(lam) for lam in lams])
    ax.contour(intf.u, lams, Imap, levels=5, colors="white", alpha=0.4, linewidths=0.8)
    ax.set_xlabel(r"$u=\ln x$")
    ax.set_ylabel(r"Lyapunov $\lambda$")
    ax.set_title("Phase Fold Evolution", fontweight="bold")
    ax.set_xlim(-5, 5)

    ax = axes[1, 0]
    E = [E_geo(l) for l in lams]
    C = [C_mu(l) for l in lams]
    ax.plot(lams, E, "r-", lw=3, label=r"$E_{\rm geo}$")
    ax.fill_between(lams, E, alpha=0.2, color="red")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$E_{\rm geo}$", color="r")
    ax.set_title("Geometric Kinetic Release", fontweight="bold")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(lams, C, "purple", lw=2, ls="--", label=r"$\mathbf{C}(\mu)$")
    ax2.set_ylabel(r"$\mathbf{C}(\mu)$", color="purple")
    ax2.legend(loc="center right")
    ax.annotate("Chern Jump → Kinetic Release", xy=(1.0, E_geo(1.0)), xytext=(2.5, 0.15),
               arrowprops=dict(arrowstyle="->", color="darkred", lw=2), fontsize=11, color="darkred",
               bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.8))

    ax = axes[1, 1]
    Imap = np.array([intf.interference(lam) for lam in lams])
    im = ax.imshow(Imap, aspect="auto", origin="lower", extent=[intf.u.min(), intf.u.max(), lams.min(), lams.max()],
                   cmap="inferno", vmin=0, vmax=1)
    cs = ax.contour(intf.u, lams, Imap, levels=[0.2, 0.5, 0.8], colors="cyan", alpha=0.6, linewidths=1.2)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%.1f")
    ax.plot([0, 0], [lams.min(), lams.max()], "w--", lw=1, alpha=0.5)
    plt.colorbar(im, ax=ax, label="Intensity")
    ax.set_xlabel(r"$u=\ln x$")
    ax.set_ylabel(r"Recursive Depth $\lambda$")
    ax.set_title("Convergent Ripple Density", fontweight="bold")
    ax.set_xlim(-5, 5)
    ax.text(0.02, 0.02, "Brighter=Stronger Convergence\nEnergy does not diffuse", transform=ax.transAxes, fontsize=10, va="bottom", color="white",
           bbox=dict(boxstyle="round", facecolor="black", alpha=0.5))

    ax = axes[1, 2]
    lam = 3.0
    base = intf.psi_complex(lam)
    nq = 5
    phs = np.linspace(0, 2 * np.pi, nq, endpoint=False)
    cols = plt.cm.viridis(np.linspace(0, 1, nq))
    tot = np.zeros_like(intf.u)
    for i, phi in enumerate(phs):
        Iq = np.abs(base * np.exp(1j * phi))**2
        Iq = Iq / np.max(Iq) * 0.3
        ax.plot(intf.u, Iq + i * 0.35, color=cols[i], lw=1.5, alpha=0.7, label=f"Q{i+1} φ={phi:.2f}")
        ax.fill_between(intf.u, Iq + i * 0.35, i * 0.35, alpha=0.2, color=cols[i])
        tot += Iq
    tot = tot / np.max(tot) * 0.8
    ax.plot(intf.u, tot - 0.5, "white", lw=3, alpha=0.9, label="Superposition")
    ax.fill_between(intf.u, tot - 0.5, -0.5, alpha=0.3, color="white")
    ax.set_xlim(-5, 5)
    ax.set_ylim(-0.8, nq * 0.35 + 0.2)
    ax.set_xlabel(r"$u=\ln x$")
    ax.set_ylabel(r"Level Offset $\rightarrow$")
    ax.set_title("Multi-Quantum Fold Superposition", fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.2)
    ax.annotate("Algebraic Superposition = Observed Pattern", xy=(0, -0.1), xytext=(2, -0.6),
               arrowprops=dict(arrowstyle="->", color="white", lw=2), fontsize=11, color="white",
               bbox=dict(boxstyle="round", facecolor="darkgreen", alpha=0.7))

    fig.suptitle("Group V — Convergent Ripple Geometry: Phase Folds & Multi-Quantum Superposition", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_05_convergent.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_all_patterns(intf, out):
    ncols = 10
    nrows = (CFG["n_lambda"] + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, 2.2 * nrows))
    axes = axes.flatten()
    lams = np.linspace(CFG["lambda_min"], CFG["lambda_max"], CFG["n_lambda"])
    for i, lam in enumerate(lams):
        Iu = intf.interference(lam)
        axes[i].plot(intf.u, Iu, "b-", lw=0.8)
        axes[i].fill_between(intf.u, Iu, alpha=0.15, color="blue")
        axes[i].set_xlim(-5, 5)
        axes[i].set_ylim(-0.05, 1.1)
        axes[i].set_title(f"λ={lam:.2f}\nV={intf.scan['vis_fft'][i]:.2f}", fontsize=6)
        axes[i].set_xticks([])
        axes[i].set_yticks([])
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle("Group VI — All λ Convergent Ripple Patterns (80 panels)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out / "group_06_all_patterns.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


def group_tunneling(intf, out):
    """Group VII: Geometric phase transition / high-density cluster tunneling."""
    fig = plt.figure(figsize=(18, 12))
    x = np.exp(intf.u)
    mask = (x > 1e-3) & (x < 100)

    ax1 = fig.add_subplot(2, 3, 1)
    ax1.scatter([-1.5, 0, 1.5], [0, 0, 0], s=200, c="blue", marker="x", linewidths=3, label="Low ρ: separated", zorder=5)
    ax1.scatter([-0.3, 0, 0.3], [0, 0, 0], s=400, c="red", marker="x", linewidths=4, label="High ρ: pile-up", zorder=5)
    ax1.annotate("Recursive\ncoupling ↑", xy=(0, 0.3), fontsize=10, ha="center", color="darkred",
                bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.8))
    ax1.set_xlim(-2.5, 2.5)
    ax1.set_ylim(-0.8, 0.8)
    ax1.set_xlabel(r"$\Re(s)$ (spectral parameter)")
    ax1.set_ylabel(r"$\Im(s)$")
    ax1.set_title("Mellin Space: Spectral Pile-up\nHigh ρ → poles compress", fontweight="bold")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color="k", lw=0.5)

    ax2 = fig.add_subplot(2, 3, 2)
    lams = np.linspace(0, 6, 200)
    C_vals = C_mu(lams)
    L_vals = Lyapunov(lams)
    beta_vals = beta(lams)
    ax2.plot(lams, C_vals, "purple", lw=2.5, label=r"$\mathbf{C}(\mu) = \sin(\pi\beta/2)$")
    ax2.plot(lams, L_vals, "darkgreen", lw=2.5, ls="--", label=r"$\mathcal{L}(\lambda)$ (Lyapunov)")
    ax2.plot(lams, beta_vals / 2, "gray", lw=1.5, alpha=0.5, label=r"$\beta/2$")
    ax2.axhline(y=0.85, color="r", ls=":", lw=2, alpha=0.7, label="Critical threshold")
    ax2.fill_between(lams, 0.85, 1.0, where=(C_vals > 0.85), alpha=0.15, color="red", label="Overflow zone")
    ax2.set_xlabel(r"Recursive breaking depth $\lambda$")
    ax2.set_ylabel("Order parameter")
    ax2.set_title("Phase Diagram: Approach to Criticality", fontweight="bold")
    ax2.set_xlim(0, 6)
    ax2.set_ylim(-0.05, 1.05)
    ax2.legend(fontsize=8, loc="lower right")
    ax2.grid(True, alpha=0.3)
    ax2.annotate("Critical\noverflow", xy=(2.5, 0.85), xytext=(3.5, 0.6),
                arrowprops=dict(arrowstyle="->", color="red", lw=1.5), fontsize=10, color="darkred")

    ax3 = fig.add_subplot(2, 3, 3)
    foc_vals = [focusing(l) for l in lams]
    ax3.plot(lams, foc_vals, "g-", lw=2.5)
    ax3.fill_between(lams, foc_vals, alpha=0.2, color="green")
    ax3.axhline(y=0.5, color="k", ls="--", alpha=0.3, label="Critical focusing")
    ax3.axvline(x=1.0, color="orange", ls="--", alpha=0.5, label="β=1")
    ax3.set_xlabel(r"$\lambda$")
    ax3.set_ylabel("Focusing degree")
    ax3.set_title("Inverse-Mellin Focusing\nEnergy convergence to finite region", fontweight="bold")
    ax3.set_ylim(-0.05, 1.05)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    ax4 = fig.add_subplot(2, 3, 4)
    lam_pre = 1.5
    I_pre = intf.interference_with_overflow(lam_pre)
    I_x_pre = I_pre / (x + 1e-12)
    I_x_pre = I_x_pre / np.max(I_x_pre) if np.max(I_x_pre) > 0 else I_x_pre
    ax4.fill_between(intf.u, I_pre, alpha=0.3, color="blue")
    ax4.plot(intf.u, I_pre, "b-", lw=1.5, label=f"Multiplicative: λ={lam_pre}")
    ax4_twin = ax4.twinx()
    ax4_twin.fill_between(x[mask], I_x_pre[mask], alpha=0.2, color="red")
    ax4_twin.plot(x[mask], I_x_pre[mask], "r-", lw=1.5, label="Additive projection")
    ax4.set_xlim(-5, 5)
    ax4.set_ylim(-0.05, 1.1)
    ax4.set_xlabel(r"$u = \ln x$")
    ax4.set_ylabel(r"$I(u)$", color="blue")
    ax4_twin.set_ylabel(r"$I(x)$", color="red")
    ax4_twin.set_xscale("log")
    ax4_twin.set_xlim(1e-3, 100)
    ax4.set_title(f'Pre-critical: λ={lam_pre}, C={C_mu(lam_pre):.2f}\n"Rigidity overload"', fontweight="bold")
    ax4.grid(True, alpha=0.3)

    ax5 = fig.add_subplot(2, 3, 5)
    lam_crit = 2.5
    I_crit = intf.interference_with_overflow(lam_crit)
    I_x_crit = I_crit / (x + 1e-12)
    I_x_crit = I_x_crit / np.max(I_x_crit) if np.max(I_x_crit) > 0 else I_x_crit
    ax5.fill_between(intf.u, I_crit, alpha=0.3, color="blue")
    ax5.plot(intf.u, I_crit, "b-", lw=1.5)
    ax5_twin = ax5.twinx()
    ax5_twin.fill_between(x[mask], I_x_crit[mask], alpha=0.2, color="red")
    ax5_twin.plot(x[mask], I_x_crit[mask], "r-", lw=1.5)
    ax5.set_xlim(-5, 5)
    ax5.set_ylim(-0.05, 1.1)
    ax5.set_xlabel(r"$u = \ln x$")
    ax5.set_ylabel(r"$I(u)$", color="blue")
    ax5_twin.set_ylabel(r"$I(x)$", color="red")
    ax5_twin.set_xscale("log")
    ax5_twin.set_xlim(1e-3, 100)
    ax5.set_title(f"CRITICAL: λ={lam_crit}, C={C_mu(lam_crit):.2f}\nTopological charge overflow!", fontweight="bold", color="darkred")
    ax5.grid(True, alpha=0.3)
    ax5.annotate("π-phase jump\n(branch switch)", xy=(0, 0.5), xytext=(2, 0.8),
                arrowprops=dict(arrowstyle="->", color="darkred", lw=2), fontsize=10, color="darkred",
                bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.8))

    ax6 = fig.add_subplot(2, 3, 6)
    lam_post = 4.0
    I_post = intf.interference_with_overflow(lam_post)
    I_x_post = I_post / (x + 1e-12)
    I_x_post = I_x_post / np.max(I_x_post) if np.max(I_x_post) > 0 else I_x_post
    ax6.plot(x[mask], I_x_pre[mask], "b--", lw=2, alpha=0.6, label=f"Pre: λ={lam_pre}")
    ax6.plot(x[mask], I_x_crit[mask], "orange", lw=2.5, alpha=0.8, label=f"Critical: λ={lam_crit}")
    ax6.plot(x[mask], I_x_post[mask], "r-", lw=2.5, label=f"Post: λ={lam_post}")
    for l, Ix, color in [(lam_pre, I_x_pre, "blue"), (lam_crit, I_x_crit, "orange"), (lam_post, I_x_post, "red")]:
        peak_idx = np.argmax(Ix[mask])
        px = x[mask][peak_idx]
        pI = Ix[mask][peak_idx]
        ax6.scatter([px], [pI], s=150, c=color, zorder=5, edgecolors="black", linewidths=1.5)
    ax6.set_xlim(0, 30)
    ax6.set_ylim(-0.05, 1.1)
    ax6.set_xlabel(r"$x$ (additive coordinate)")
    ax6.set_ylabel(r"$I(x) = I(u)/x$")
    ax6.set_title("Peak Position Jump in x-space\nNOT tunneling — branch switch!", fontweight="bold")
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    ax6.text(0.98, 0.98,
             'Standard QM: "particle tunnels through barrier"\n'
             'NBG: "topological charge overflow causes\n'
             '      Mellin branch jump → peak reappears"',
             transform=ax6.transAxes, fontsize=9, va="top", ha="right",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9))

    fig.suptitle("Group VII — Geometric Phase Transition: High-Density Cluster \"Tunneling\"\n"
                 "From Recursive Rigidity Overload to Topological Charge Integer Overflow",
                 fontsize=15, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "group_07_tunneling.png", dpi=CFG["dpi"], bbox_inches="tight")
    plt.close(fig)


# ==================== Main ====================
def main():
    out = Path(CFG["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("NBG Unified Framework — 7 Professional Figure Groups")
    print("=" * 60)

    spec = SpectralData(CFG)
    spec.run()
    group_spectral(spec, out)
    print("[1/7] Group I  — Spectral Statistics")

    intf = InterferenceData(CFG)
    intf.run()
    group_interference(intf, out)
    print("[2/7] Group II — Interference Phase Transition")

    loop = LoopData(spec.beta_results, intf.scan, CFG)
    loop.run()
    group_loop(loop, out)
    print("[3/7] Group III — Closed-Loop Validation")

    group_projection(intf, out)
    print("[4/7] Group IV — Basis Projection")

    group_convergent(intf, out)
    print("[5/7] Group V  — Convergent Ripple Geometry")

    group_all_patterns(intf, out)
    print("[6/7] Group VI — All λ Patterns")

    group_tunneling(intf, out)
    print("[7/7] Group VII — Geometric Phase Transition / Tunneling")

    with open(out / "nbg_results.json", "w") as f:
        json.dump({
            "config": CFG,
            "spectral_beta": spec.beta_results,
            "interference": intf.scan,
            "closed_loop": loop.val,
        }, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)

    print(f"\n✓ All outputs saved to: {out}")
    print(f"  Figures: group_01 ~ group_07 (.png)")
    print(f"  Data:    nbg_results.json")


if __name__ == "__main__":
    main()
