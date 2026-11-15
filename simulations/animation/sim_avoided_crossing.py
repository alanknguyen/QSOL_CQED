#!/usr/bin/env python3
"""
Dressed-state avoided crossing in the Jaynes-Cummings model.

In the frame rotating at the cavity frequency the manifold {|e,n>, |g,n+1>} has
energies  E_+-(n) = +- (1/2) Omega_n(Delta),  Omega_n = sqrt(Delta^2 + 4 g^2 (n+1)),
with Delta = omega_a - omega_c.  The bare levels +-Delta/2 cross at Delta = 0;
the coupling opens a gap 2 g sqrt(n+1) there.

The closed form is checked against a numerical diagonalisation of the QuTiP
Hamiltonian (markers in the figure).

Generates: figures/fig_avoided_crossing (.png, .pdf), animations/anim_avoided_crossing.gif

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, save, panel_label, fig_to_frame, save_gif, FULL_W,
                       COLORS, FIG_DIR, ANIM_DIR, PAPER_FIG, SUPP_DIR)

use_style()

g = 1.0
N_cav = 40
n_values = [0, 1, 5, 10]
colors_n = [COLORS['blue'], COLORS['orange'], COLORS['green'], COLORS['red']]
styles_n = ['-', (0, (6, 2)), (0, (3, 1, 1, 1)), (0, (1, 1.2))]
Delta = np.linspace(-10, 10, 401)


def omega_n(delta, n):
    """Generalised n-photon Rabi frequency."""
    return np.sqrt(delta ** 2 + 4 * g ** 2 * (n + 1))


# -- numerical check: diagonalise the 2x2 block of H = (Delta/2) sigma_z + g (a^dag sigma- + a sigma+)
a = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
sz = qt.tensor(qt.qeye(N_cav), qt.sigmaz())
H_int = (g * (a.dag() * sm + a * sm.dag())).full()
Sz = sz.full()


def numeric_pair(delta, n):
    """Eigenvalues of the manifold {|e,n>, |g,n+1>} (basis index = 2*n_photon + atom)."""
    idx = [2 * n, 2 * (n + 1) + 1]
    block = (0.5 * delta * Sz + H_int)[np.ix_(idx, idx)]
    return np.linalg.eigvalsh(block)


Delta_chk = np.linspace(-10, 10, 9)
max_err = max(abs(numeric_pair(d, n)[1] - 0.5 * omega_n(d, n)) for d in Delta_chk for n in range(21))
n_sweep = np.arange(0, 21)
split_num = np.array([np.diff(numeric_pair(0.0, n))[0] for n in n_sweep])
print(f"closed form vs numerical diagonalisation: max |error| = {max_err:.1e}")
assert max_err < 1e-12


def draw_levels(ax, alpha=1.0, lw=1.2):
    ax.plot(Delta, 0.5 * Delta, color='0.55', lw=0.6, ls='--')
    ax.plot(Delta, -0.5 * Delta, color='0.55', lw=0.6, ls='--')
    for n, c, ls in zip(n_values, colors_n, styles_n):
        E = 0.5 * omega_n(Delta, n)
        ax.plot(Delta, E, color=c, lw=lw, ls=ls, alpha=alpha, label=rf'$n = {n}$')
        ax.plot(Delta, -E, color=c, lw=lw, ls=ls, alpha=alpha)
    ax.set_xlabel(r'$\Delta/g$')
    ax.set_ylabel(r'$E_\pm(n)/\hbar g$')
    ax.set_xlim(-10, 10); ax.set_ylim(-7, 7)
    ax.set_xticks([-10, -5, 0, 5, 10])


# ---------------------------------------------------------------------------
# Static figure
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FULL_W, 2.7), gridspec_kw={'width_ratios': [1.45, 1]})
draw_levels(ax1)
for n, c in zip(n_values, colors_n):        # numerical eigenvalues
    ev = np.array([numeric_pair(d, n) for d in Delta_chk])
    ax1.plot(Delta_chk, ev[:, 1], 'o', ms=2.6, mfc='white', mec=c, mew=0.7)
    ax1.plot(Delta_chk, ev[:, 0], 'o', ms=2.6, mfc='white', mec=c, mew=0.7)
ax1.annotate('', xy=(0, 1.0), xytext=(0, -1.0),
             arrowprops=dict(arrowstyle='<->', lw=0.7, shrinkA=0, shrinkB=0))
ax1.text(0.35, 0.0, r'$2g$', fontsize=8, va='center')
ax1.text(-9.5, -4.2, 'bare', fontsize=7.5, color='0.4', rotation=24)
ax1.legend(loc='upper center', ncol=4, fontsize=7.5, columnspacing=1.0, handlelength=2.0)
panel_label(ax1, 'a')

ax2.plot(n_sweep, 2 * g * np.sqrt(n_sweep + 1), color='k', lw=0.9, label=r'$2g\sqrt{n+1}$')
ax2.plot(n_sweep, split_num, 'o', ms=3, mfc='white', mec='k', mew=0.7, label='numerical')
for n, c in zip(n_values, colors_n):
    ax2.plot(n, 2 * g * np.sqrt(n + 1), 'o', ms=4.5, color=c, zorder=5)
ax2.set_xlabel(r'photon number $n$')
ax2.set_ylabel(r'gap at $\Delta = 0$, $\Omega_n(0)/g$')
ax2.set_xticks([0, 5, 10, 15, 20])
ax2.legend(loc='lower right')
panel_label(ax2, 'b')
save(fig, FIG_DIR / 'fig_avoided_crossing', copy_to=(PAPER_FIG, SUPP_DIR))


# ---------------------------------------------------------------------------
# Animation: sweep Delta from -10 g to +10 g
# ---------------------------------------------------------------------------
Delta_sweep = np.linspace(-10, 10, 81)
frames, durations = [], []
print("Rendering avoided-crossing animation ...")
for D in Delta_sweep:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.3), gridspec_kw={'width_ratios': [1.5, 1]})
    draw_levels(ax1, alpha=0.45)
    ax1.axvline(D, color='0.4', lw=0.7)
    for n, c in zip(n_values, colors_n):
        E = 0.5 * omega_n(D, n)
        ax1.plot([D, D], [E, -E], 'o', color=c, ms=5, mec='k', mew=0.5, zorder=10)
    ax1.legend(loc='upper center', ncol=4, fontsize=7.5, columnspacing=1.0, handlelength=2.0)
    ax1.set_title(rf'dressed-state energies, $\Delta/g = {D:+.1f}$', fontsize=9)

    # right: the four pairs of levels at the current detuning
    for k, (n, c) in enumerate(zip(n_values, colors_n)):
        Om = omega_n(D, n)
        ax2.plot([k - 0.3, k + 0.3], [Om / 2, Om / 2], color=c, lw=2)
        ax2.plot([k - 0.3, k + 0.3], [-Om / 2, -Om / 2], color=c, lw=2)
        ax2.annotate('', xy=(k, Om / 2), xytext=(k, -Om / 2),
                     arrowprops=dict(arrowstyle='<->', color=c, lw=0.8, shrinkA=1, shrinkB=1))
        ax2.text(k, Om / 2 + 0.35, rf'${Om:.2f}$', ha='center', fontsize=7.5, color=c)
    ax2.set_xticks(range(len(n_values))); ax2.set_xticklabels([rf'$n={n}$' for n in n_values])
    ax2.set_xlim(-0.6, len(n_values) - 0.4); ax2.set_ylim(-7.6, 8.6)
    ax2.set_ylabel(r'$E_\pm(n)/\hbar g$')
    ax2.set_title(r'splitting $\Omega_n(\Delta)/g$', fontsize=9)
    frames.append(fig_to_frame(fig, dpi=120))
    durations.append(400 if abs(D) < 0.2 else (140 if abs(D) < 2 else 80))
durations[0] = durations[-1] = 800
save_gif(frames, ANIM_DIR / 'anim_avoided_crossing.gif', durations)
