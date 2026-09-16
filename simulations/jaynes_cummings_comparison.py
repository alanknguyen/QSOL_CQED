# jaynes_cummings_comparison.py
# Population inversion W(t) = sum_n P(n) cos(2 g sqrt(n+1) t), g = 1
# for coherent (Poisson) and thermal (Bose-Einstein) field states.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: jaynes_cummings_coherent.png, jaynes_cummings_thermal.png,
#            jaynes_cummings_comparison.png

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

PAPER_FIG = Path(__file__).resolve().parents[1] / 'paper' / 'figures'
from scipy.stats import poisson


def inversion_coherent(n_avg, t_max=100, num_points=1000):
    """W(t) for coherent state: P(n) = exp(-n_avg) * n_avg^n / n!"""
    t = np.linspace(0, t_max, num_points)
    # NOTE: an integer power of an integer array overflows int64 above
    # n ~ 16 and silently corrupts the weights; use the float pmf instead.
    n_max = int(n_avg * 5) + 20
    ns = np.arange(n_max)
    p_n = poisson.pmf(ns, n_avg)

    W = np.zeros_like(t)
    for n in range(n_max):
        W += p_n[n] * np.cos(2 * np.sqrt(n + 1) * t)
    return t, W


def inversion_thermal(n_avg, t_max=100, num_points=1000):
    """W(t) for thermal state: P(n) = n_avg^n / (1 + n_avg)^(n+1)"""
    t = np.linspace(0, t_max, num_points)
    n_max = int(n_avg * 10) + 20
    ns = np.arange(n_max, dtype=float)
    # log-space to avoid float overflow of n_avg**n at large n
    p_n = (np.exp(ns * np.log(n_avg) - (ns + 1) * np.log1p(n_avg))
           if n_avg > 0 else np.eye(1, n_max)[0])

    W = np.zeros_like(t)
    for n in range(n_max):
        W += p_n[n] * np.cos(2 * np.sqrt(n + 1) * t)
    return t, W


def plot_single(n_avg_values, state_type, filename=None):
    """Inversion vs time for one field state type."""
    num = len(n_avg_values)
    fig = plt.figure(figsize=(10, 2.5 * num))
    gs = GridSpec(num, 1, figure=fig, hspace=0.3)

    calc = inversion_coherent if state_type == 'coherent' else inversion_thermal
    label = "Coherent State" if state_type == 'coherent' else "Thermal State"

    for i, n_avg in enumerate(n_avg_values):
        ax = fig.add_subplot(gs[i, 0])
        t, W = calc(n_avg)
        ax.plot(t, W, 'k-', linewidth=1.5)
        ax.set_xlim(0, 100)
        ax.set_ylim(-1, 1)
        ax.set_ylabel("Inversion", fontsize=12)
        if i == num - 1:
            ax.set_xlabel(r"$gt$", fontsize=12)
        ax.text(0.98, 0.9, rf"$\langle n \rangle = {n_avg}$",
                ha='right', transform=ax.transAxes, fontsize=12)

    fig.suptitle(f"Population Inversion ({label})", fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    if filename:
        fig.savefig(filename, dpi=300, bbox_inches='tight')


def plot_comparison(n_avg_values, filename=None):
    """Coherent vs thermal side by side."""
    num = len(n_avg_values)
    fig = plt.figure(figsize=(15, 2.5 * num))
    gs = GridSpec(num, 2, figure=fig, hspace=0.3, wspace=0.15)

    for i, n_avg in enumerate(n_avg_values):
        # Coherent
        ax1 = fig.add_subplot(gs[i, 0])
        t, W = inversion_coherent(n_avg)
        ax1.plot(t, W, 'k-', linewidth=1.5)
        ax1.set_xlim(0, 100)
        ax1.set_ylim(-1, 1)
        ax1.set_ylabel("Inversion", fontsize=12)
        if i == 0:
            ax1.set_title("Coherent State", fontsize=12)
        if i == num - 1:
            ax1.set_xlabel(r"$gt$", fontsize=12)
        ax1.text(0.98, 0.9, rf"$\langle n \rangle = {n_avg}$",
                 ha='right', transform=ax1.transAxes, fontsize=12)

        # Thermal
        ax2 = fig.add_subplot(gs[i, 1])
        t, W = inversion_thermal(n_avg)
        ax2.plot(t, W, 'k-', linewidth=1.5)
        ax2.set_xlim(0, 100)
        ax2.set_ylim(-1, 1)
        if i == 0:
            ax2.set_title("Thermal State", fontsize=12)
        if i == num - 1:
            ax2.set_xlabel(r"$gt$", fontsize=12)
        ax2.text(0.98, 0.9, rf"$\langle n \rangle = {n_avg}$",
                 ha='right', transform=ax2.transAxes, fontsize=12)

    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    n_values = [4, 9, 14, 19, 24]
    PAPER_FIG.mkdir(parents=True, exist_ok=True)
    plot_single(n_values, 'coherent', PAPER_FIG / "jaynes_cummings_coherent.png")
    plot_single(n_values, 'thermal', PAPER_FIG / "jaynes_cummings_thermal.png")
    plot_comparison(n_values, PAPER_FIG / "jaynes_cummings_comparison.png")
    print("Saved: jaynes_cummings_*.png")
