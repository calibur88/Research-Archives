#!/usr/bin/env python3
"""
Riemann Zeros Asymptotic Analysis
Standard implementation based on von Mangoldt formula
"""

import numpy as np
import matplotlib.pyplot as plt
import os

SAVE_PATH = "/storage/emulated/0/研究报告/黎曼猜想几何版/"
DATA_FILE = os.path.join(SAVE_PATH, "riemann_zeros_analysis.txt")
CHART_FILE = os.path.join(SAVE_PATH, "riemann_analysis.png")
REF_FILE = os.path.join(SAVE_PATH, "zeros1.txt")

class ZeroAnalyzer:
    """Standard asymptotic analyzer for Riemann zeta zeros"""
    
    def __init__(self):
        self.tau = 2 * np.pi
        
    def asymptotic_sequence(self, t0, n_max):
        """Generate zeros using asymptotic formula N(t) ~ (t/2π)ln(t/2π)"""
        if n_max <= 0:
            return np.array([])
        
        seq = [float(t0)]
        t = t0
        
        for _ in range(n_max - 1):
            if t <= self.tau:
                break
            # Standard asymptotic spacing: 2π/ln(t/2π)
            dt = self.tau / np.log(t / self.tau)
            t += dt
            seq.append(t)
        
        return np.array(seq)
    
    def load_odlyzko(self, filepath, limit=10000):
        """Load reference zeros from Odlyzko database"""
        if not os.path.exists(filepath):
            return None
        
        zeros = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if len(zeros) >= limit:
                        break
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    try:
                        val = float(line.split()[-1])
                        zeros.append(val)
                    except:
                        continue
            return np.array(zeros) if zeros else None
        except:
            return None

def calc_error_stats(pred, true):
    """Standard error metrics"""
    n = min(len(pred), len(true))
    if n == 0:
        return None
    p, t = pred[:n], true[:n]
    abs_err = np.abs(p - t)
    rel_err = abs_err / t * 100
    return {
        'mae': np.mean(abs_err),
        'max_err': np.max(abs_err),
        'mean_rel': np.mean(rel_err),
        'corr': np.corrcoef(p, t)[0,1] if n > 1 else 0
    }

def generate_figure(indices, pred, true=None):
    """
    Standard 4-panel Riemann zero analysis figure
    """
    if len(pred) == 0:
        return None
    
    if indices is None or len(indices) != len(pred):
        indices = np.arange(1, len(pred)+1)
    
    has_true = (true is not None) and (len(true) > 0)
    if has_true:
        true = np.array(true)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Riemann Zeta Zeros - Asymptotic Analysis', fontsize=16, fontweight='bold', y=0.98)
    
    # Panel 1: Position comparison
    ax1 = axes[0, 0]
    ax1.plot(indices, pred, 'b-', linewidth=2, label='Asymptotic Formula', alpha=0.8)
    if has_true:
        cmp_len = min(len(true), len(pred))
        ax1.scatter(indices[:cmp_len], true[:cmp_len], c='red', s=12, 
                   alpha=0.6, label='Odlyzko Data', zorder=5)
        errors = np.abs(pred[:cmp_len] - true[:cmp_len])
        ax1.fill_between(indices[:cmp_len], pred[:cmp_len], true[:cmp_len], 
                        alpha=0.2, color='gray', label=f'Error (MAE={np.mean(errors):.2f})')
    ax1.set_xlabel('Index n')
    ax1.set_ylabel('Zero position t_n')
    ax1.set_title('Asymptotic Approximation vs Odlyzko Data')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Error analysis (S(t) term)
    ax2 = axes[0, 1]
    if has_true:
        cmp_len = min(len(true), len(pred))
        err = np.abs(pred[:cmp_len] - true[:cmp_len])
        ax2.plot(indices[:cmp_len], err, 'purple', linewidth=1.5, label='|Δt|')
        ax2.axhline(y=np.mean(err), color='red', linestyle='--', 
                   label=f'MAE = {np.mean(err):.3f}')
        ax2.fill_between(indices[:cmp_len], 0, err, alpha=0.3, color='purple')
        ax2.set_ylabel('Absolute Error')
        ax2.set_title(f'Approximation Error (r = {np.corrcoef(pred[:cmp_len], true[:cmp_len])[0,1]:.6f})')
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, 'No reference data available', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=12)
    ax2.set_xlabel('Index n')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Level spacing distribution
    ax3 = axes[1, 0]
    if len(pred) > 1:
        p_gaps = np.diff(pred)
        ax3.plot(indices[1:], p_gaps, 'b-', linewidth=2, label='Computed gaps', alpha=0.7)
        if has_true and len(true) > 1:
            cmp_len = min(len(true), len(pred))
            t_gaps = np.diff(true[:cmp_len])
            ax3.scatter(indices[1:cmp_len], t_gaps, c='red', s=8, alpha=0.6, label='Observed gaps')
        # Theoretical: 2π/ln(t/2π)
        theory_gaps = 2*np.pi / np.log(pred[:-1]/(2*np.pi))
        ax3.plot(indices[1:], theory_gaps, 'g--', linewidth=1.5, 
                label='von Mangoldt: 2π/ln(t/2π)')
        ax3.set_xlabel('Index n')
        ax3.set_ylabel('Gap Δt')
        ax3.set_title('Consecutive Zero Spacing')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # Panel 4: Nearest-neighbor statistics (GUE)
    ax4 = axes[1, 1]
    src = true if has_true else pred
    if len(src) > 3:
        gaps = np.diff(src)
        ratios = []
        for i in range(len(gaps)-1):
            if gaps[i] > 0 and gaps[i+1] > 0:
                r = min(gaps[i], gaps[i+1]) / max(gaps[i], gaps[i+1])
                ratios.append(r)
        
        if ratios:
            ax4.hist(ratios, bins=20, density=True, alpha=0.7, color='purple',
                    label=f'Measured (n={len(ratios)}, μ={np.mean(ratios):.3f})')
            ax4.axvline(0.536, color='red', linestyle='--', linewidth=2, label='GUE 0.536')
            ax4.axvline(0.386, color='gray', linestyle=':', linewidth=2, label='Poisson 0.386')
            ax4.set_xlabel('Spacing ratio min/max')
            ax4.set_ylabel('Density')
            ax4.set_title('Nearest-Neighbor Spacing Statistics')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
            ax4.set_xlim(0, 1)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig

def save_analysis_data(indices, pred, true, stats, filepath):
    """Save data in standard format"""
    with open(filepath, 'w', encoding='utf-8') as f:
        # Minimal header - standard math only
        f.write("# Riemann Zeta Zeros - Asymptotic Analysis\n")
        f.write("# Formula: t_{n+1} = t_n + 2π/ln(t_n/2π)\n")
        f.write("# Reference: Odlyzko database\n\n")
        
        total = len(pred)
        ref_n = len(true) if true is not None else 0
        
        f.write(f"Total_points: {total}\n")
        f.write(f"Reference_points: {ref_n}\n")
        if stats:
            f.write(f"Mean_Absolute_Error: {stats['mae']:.8f}\n")
            f.write(f"Correlation: {stats['corr']:.8f}\n")
        f.write("\n")
        
        # Data columns
        header = f"{'n':>6} {'t_n_computed':>18}"
        if ref_n > 0:
            header += f" {'t_n_odlyzko':>18} {'error':>12}"
        header += f" {'delta_t':>12}\n"
        f.write(header)
        f.write("-" * 85 + "\n")
        
        for i in range(total):
            delta = pred[i] - pred[i-1] if i > 0 else 0.0
            line = f"{int(indices[i]):6d} {pred[i]:18.10f}"
            if ref_n > 0:
                if i < ref_n:
                    err = pred[i] - true[i]
                    line += f" {true[i]:18.10f} {err:12.6f}"
                else:
                    line += f" {'-':>18} {'-':>12}"
            line += f" {delta:12.6f}\n"
            f.write(line)

def main():
    print("="*60)
    print("Riemann Zeta Zero Analysis")
    print("Asymptotic Formula: N(t) ~ (t/2π)ln(t/2π)")
    print("="*60)
    
    if not os.path.exists(SAVE_PATH):
        os.makedirs(SAVE_PATH)
    
    analyzer = ZeroAnalyzer()
    
    # Generate 500 zeros
    N = 500
    t0 = 14.134725  # First zero
    
    print(f"\nGenerating {N} zeros...")
    pred = analyzer.asymptotic_sequence(t0, N)
    indices = np.arange(1, len(pred)+1)
    print(f"Computed: {len(pred)} zeros")
    
    # Load reference
    print("Loading Odlyzko reference...")
    ref = analyzer.load_odlyzko(REF_FILE, limit=N)
    if ref is not None:
        print(f"Loaded: {len(ref)} reference zeros")
    else:
        print("No reference file found")
    
    # Stats
    stats = None
    if ref is not None and len(ref) > 0:
        stats = calc_error_stats(pred, ref)
        print(f"\nError Statistics:")
        print(f"  MAE: {stats['mae']:.4f}")
        print(f"  Correlation: {stats['corr']:.6f}")
    
    # Generate figure
    print("\nGenerating analysis figure...")
    fig = generate_figure(indices, pred, ref)
    if fig:
        fig.savefig(CHART_FILE, dpi=300, bbox_inches='tight')
        print(f"Figure: {CHART_FILE}")
    
    # Save data
    save_analysis_data(indices, pred, ref, stats, DATA_FILE)
    print(f"Data: {DATA_FILE}")
    
    print(f"\nOutput complete. {len(pred)} zeros analyzed.")

if __name__ == "__main__":
    main()
