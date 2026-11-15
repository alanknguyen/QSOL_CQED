#!/usr/bin/env python3
"""
Animation for "Rotation and shear in the Jaynes-Cummings model" (paper/ejp/jc_beyond_linear).

The half revival is a revival of the photon-number parity, and whether the cat state is
"odd" or "even" is a matter of the Rabi phase.  The mean photon number n_bar is swept
from 8.5 to 11.5 and the field is looked at each time at t = t_r/2 = pi sqrt(n_bar)/g:
  - left: the photon-number distribution P(n) (even n blue, odd n orange).  It is a comb
    of even numbers at n_bar = 9 and 11 and of odd numbers at n_bar = 10;
  - right: the parity <(-1)^n> at t_r/2 against n_bar, exact (solid) and the closed form
    -cos(pi n_bar - eta) / (1 + pi^2/4)^(1/4), eta = arctan(pi/2)/2 = 0.502 (dashed).
    The amplitude 0.733 and the phase eta are what the shear does to the ideal value
    -cos(pi n_bar).

Everything is evaluated from the exact solution, P(n,t) = p_n cos^2(theta_n) +
p_{n-1} sin^2(theta_{n-1}) with theta_n = g t sqrt(n+1).

Generates: animations/anim_parity_comb.gif
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import use_style, fig_to_frame, save_gif, COLORS, ANIM_DIR

use_style()
g = 1.0
PI = np.pi
ETA = 0.5 * np.arctan(PI / 2)
AMP = (1 + PI ** 2 / 4) ** -0.25
N_SHOW = 27
n = np.arange(200)


def distribution(nb):
    """Photon-number distribution of the field at t_r/2."""
    p = poisson.pmf(n, nb)
    th = PI * np.sqrt(nb) * np.sqrt(n + 1.0)
    return p * np.cos(th) ** 2 + np.concatenate([[0.0], (p * np.sin(th) ** 2)[:-1]])


def parity(nb):
    pn = distribution(nb)
    return pn[0::2].sum() - pn[1::2].sum()


def closed_form(nb):
    return -np.cos(PI * nb - ETA) * AMP


nb_frames = np.round(np.arange(8.5, 11.5001, 0.025), 3)        # integers 9, 10, 11 are frames
nb_line = np.linspace(8.4, 11.6, 641)
exact_line = np.array([parity(x) for x in nb_line])
print("parity at t_r/2, exact / closed form:")
for x in (9, 10, 10.5, 11):
    print(f"   n_bar = {x:4.1f}: {parity(x):+.4f} / {closed_form(x):+.4f}")
assert max(abs(parity(x) - closed_form(x)) for x in (9, 10, 11)) < 3e-3

even, odd = n[:N_SHOW] % 2 == 0, n[:N_SHOW] % 2 == 1
frames, durations = [], []
for nb in nb_frames:
    pn = distribution(nb)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.8, 3.2), gridspec_kw={'width_ratios': [1.25, 1]})
    ax1.bar(n[:N_SHOW][even], pn[:N_SHOW][even], width=0.8, color=COLORS['blue'], lw=0, label='even $n$')
    ax1.bar(n[:N_SHOW][odd], pn[:N_SHOW][odd], width=0.8, color=COLORS['orange'], lw=0, label='odd $n$')
    ax1.plot(n[:N_SHOW], poisson.pmf(n[:N_SHOW], nb), color='k', lw=0.6, ls='--', marker='.', ms=2, label='Poisson')
    ax1.set_xlim(-0.8, N_SHOW - 0.2); ax1.set_ylim(0, 0.27)
    ax1.set_xlabel(r'$n$'); ax1.set_ylabel(r'$P(n)$ at $t_r/2$')
    ax1.legend(loc='upper right', fontsize=7.5, handlelength=1.2)
    kind = 'even comb' if parity(nb) > 0.45 else ('odd comb' if parity(nb) < -0.45 else 'no comb')
    ax1.set_title(rf'$\bar{{n}} = {nb:.2f}$: {kind}', fontsize=9, loc='left')

    ax2.axhline(0, color='0.6', lw=0.4)
    for sgn in (+1, -1):
        ax2.axhline(sgn * AMP, color='0.5', lw=0.5, ls=':')
    ax2.text(8.45, AMP + 0.04, r'$\pm(1+\pi^2/4)^{-1/4}$', fontsize=7, color='0.35')
    ax2.plot(nb_line, exact_line, color=COLORS['purple'], lw=1.1, label='exact')
    ax2.plot(nb_line, closed_form(nb_line), color='k', lw=0.8, ls='--', label='closed form')
    ax2.plot(nb, parity(nb), 'o', color=COLORS['red'], ms=5.5, mec='k', mew=0.5, zorder=5)
    ax2.set_xlim(8.4, 11.6); ax2.set_ylim(-1.0, 1.12)
    ax2.set_xticks([9, 10, 11])
    ax2.set_xlabel(r'mean photon number $\bar{n}$'); ax2.set_ylabel(r'parity $\langle(-1)^{\hat n}\rangle$ at $t_r/2$')
    ax2.legend(loc='lower left', fontsize=7.5, handlelength=1.6, ncol=2, columnspacing=1.0)
    ax2.set_title(rf'$\langle(-1)^{{\hat n}}\rangle = {parity(nb):+.2f}$', fontsize=9)
    frames.append(fig_to_frame(fig, dpi=115))
    durations.append(1700 if abs(nb - round(nb)) < 1e-9 else 80)
durations[0] = durations[-1] = 900
save_gif(frames, ANIM_DIR / 'anim_parity_comb.gif', durations)
