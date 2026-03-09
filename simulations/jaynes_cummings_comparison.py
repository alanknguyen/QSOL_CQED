# jaynes_cummings_comparison.py
# Population inversion W(t) = sum_n P(n) cos(2 sqrt(n+1) t)
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
from scipy.special import factorial


def inversion_coherent(n_avg, t_max=100, num_points=1000):
    """W(t) for coherent state: P(n) = exp(-n_avg) * n_avg^n / n!"""
    t = np.linspace(0, t_max, num_points)
    n_max = int(n_avg * 5)
    ns = np.arange(n_max)
    p_n = np.exp(-n_avg) * np.power(n_avg, ns) / factorial(ns)

    W = np.zeros_like(t)
    for n in range(n_max):
        W += p_n[n] * np.cos(2 * np.sqrt(n + 1) * t)
    return t, W


def inversion_thermal(n_avg, t_max=100, num_points=1000):
    """W(t) for thermal state: P(n) = n_avg^n / (1 + n_avg)^(n+1)"""
    t = np.linspace(0, t_max, num_points)
    n_max = int(n_avg * 10)

    p_n = np.zeros(n_max)
    if n_avg > 0:
        for n in range(n_max):
            p_n[n] = (n_avg**n) / ((1 + n_avg)**(n + 1))

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
            ax.set_xlabel(r"Time (units of $(2C\langle n \rangle^{1/2}/\hbar)^{-1}$)", fontsize=12)
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
            ax1.set_xlabel(r"Time (units of $(2C\langle n \rangle^{1/2}/\hbar)^{-1}$)", fontsize=12)
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
            ax2.set_xlabel(r"Time (units of $(2C\langle n \rangle^{1/2}/\hbar)^{-1}$)", fontsize=12)
        ax2.text(0.98, 0.9, rf"$\langle n \rangle = {n_avg}$",
                 ha='right', transform=ax2.transAxes, fontsize=12)

    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    n_values = [4, 9, 14, 19, 24]
    plot_single(n_values, 'coherent', "jaynes_cummings_coherent.png")
    plot_single(n_values, 'thermal', "jaynes_cummings_thermal.png")
    plot_comparison(n_values, "jaynes_cummings_comparison.png")
