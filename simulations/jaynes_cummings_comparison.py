# jaynes_cummings_comparison.py
# Atomic inversion  <sigma_z>(t) = sum_n p(n) cos(2 g sqrt(n+1) t),  g = 1,
# for coherent (Poisson) and thermal (Bose-Einstein) photon statistics.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585, taught by Prof. Luca Dal Negro
#
# Generates (paper/figures/): jaynes_cummings_comparison (.png, .pdf),
#     jaynes_cummings_coherent.png, jaynes_cummings_thermal.png
#
# The sum is evaluated on a time grid fine enough to resolve the fastest Rabi
# frequency that carries weight (at least 20 points per period); a coarser grid
# aliases the oscillations into a spurious envelope.

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

from plotstyle import use_style, save, panel_label, FULL_W, COLORS, PAPER_FIG

use_style()
T_MAX = 100.0
TAIL = 1e-10          # photon-number weight neglected in the sum


def photon_weights(n_avg, kind):
    """p(n) truncated where the neglected weight is below TAIL."""
    if kind == 'coherent':        # p(n) = exp(-n_avg) n_avg^n / n!
        n_max = int(poisson.isf(TAIL, n_avg)) + 2
        return poisson.pmf(np.arange(n_max), n_avg)
    # thermal: p(n) = n_avg^n / (1 + n_avg)^(n+1); tail beyond n is (n_avg/(1+n_avg))^n
    q = n_avg / (1.0 + n_avg)
    n_max = int(np.log(TAIL) / np.log(q)) + 2
    return (1.0 - q) * q ** np.arange(n_max)


def inversion(n_avg, kind):
    p = photon_weights(n_avg, kind)
    omega_max = 2.0 * np.sqrt(len(p))
    t = np.linspace(0, T_MAX, int(T_MAX * omega_max / (2 * np.pi) * 20) + 1)
    W = np.zeros_like(t)
    for n, pn in enumerate(p):
        W += pn * np.cos(2.0 * np.sqrt(n + 1) * t)
    return t, W


def draw_column(axes, n_values, kind, color):
    for ax, n_avg in zip(axes, n_values):
        t, W = inversion(n_avg, kind)
        ax.plot(t, W, color=color, lw=0.4, rasterized=True)
        ax.set_xlim(0, T_MAX)
        ax.set_ylim(-1.05, 1.38)          # headroom, so the label never covers a spike
        ax.set_yticks([-1, 0, 1])
        ax.text(0.985, 0.97, rf'$\bar{{n}} = {n_avg}$', ha='right', va='top',
                transform=ax.transAxes, fontsize=8,
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
        if kind == 'coherent':     # the k-th revival envelope peaks at k * 2 pi sqrt(n_avg + 1)
            t_r = 2 * np.pi * np.sqrt(n_avg + 1)
            ks = np.arange(1, int(T_MAX / t_r) + 1)
            ax.plot(ks * t_r, np.full(len(ks), -0.9), marker='^', ms=3.5, ls='none',
                    mfc='k', mec='k')
    axes[-1].set_xlabel(r'$gt$')


def plot_comparison(n_values, stem):
    fig, axes = plt.subplots(len(n_values), 2, figsize=(FULL_W, 1.08 * len(n_values) + 0.4),
                             sharex=True, sharey=True)
    draw_column(axes[:, 0], n_values, 'coherent', COLORS['blue'])
    draw_column(axes[:, 1], n_values, 'thermal', COLORS['red'])
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$\langle \sigma_z \rangle$')
    panel_label(axes[0, 0], 'a', 'coherent field')
    panel_label(axes[0, 1], 'b', 'thermal field')
    save(fig, stem)


def plot_single(n_values, kind, stem):
    fig, axes = plt.subplots(len(n_values), 1, figsize=(FULL_W, 1.0 * len(n_values) + 0.4),
                             sharex=True, sharey=True)
    draw_column(axes, n_values, kind, COLORS['blue'] if kind == 'coherent' else COLORS['red'])
    for ax in axes:
        ax.set_ylabel(r'$\langle \sigma_z \rangle$')
    fig.savefig(stem.with_suffix('.png'), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    n_values = [4, 9, 14, 19, 24]
    plot_comparison(n_values, PAPER_FIG / 'jaynes_cummings_comparison')
    plot_single(n_values, 'coherent', PAPER_FIG / 'jaynes_cummings_coherent')
    plot_single(n_values, 'thermal', PAPER_FIG / 'jaynes_cummings_thermal')
    print("Saved: jaynes_cummings_coherent.png, jaynes_cummings_thermal.png")
