#!/usr/bin/env python3
"""
Photon-number distribution P(n, t) of the cavity field during JC dynamics.

For |e> x |alpha> on resonance the distribution is known in closed form,

    P(n, t) = p(n) cos^2(g sqrt(n+1) t) + p(n-1) sin^2(g sqrt(n) t),

with p(n) the initial Poisson weights; the script checks the numerics against it.

The two components of the half-revival cat have the SAME amplitude |alpha| and
opposite phases, so their Poisson envelopes coincide and P(n) never becomes
bimodal.  What appears at t_r/2 is their interference: a parity comb.  For
n_bar = 10 odd photon numbers dominate (parity <(-1)^n> = -0.64).  Which
parity wins depends on n_bar: the comb is odd for even integer n_bar and even
for odd integer n_bar (parity +0.64 for n_bar = 9 or 11).  At the first
revival the comb has washed out into a broadened single hump.

Generates: figures/fig_photon_number_evolution (.png, .pdf), animations/anim_photon_number.gif

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson as poisson_dist
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, save, panel_label, time_guides, fig_to_frame, save_gif,
                       FULL_W, COLORS, FIG_DIR, ANIM_DIR, PAPER_FIG, SUPP_DIR)

use_style()

g = 1.0
n_bar = 10
N_cav = 40
t_r = 2 * np.pi * np.sqrt(n_bar) / g
t_c = 1.0 / g                     # operational collapse time (paper Sec. III)
t_max = 1.3 * t_r
N_SHOW = 26

a = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
sz = qt.tensor(qt.qeye(N_cav), qt.sigmaz())
H_JC = g * (a.dag() * sm + a * sm.dag())
psi0 = qt.tensor(qt.coherent(N_cav, np.sqrt(n_bar)), qt.basis(2, 0))

snapshots = [
    ('a', 0.0,          r'$t = 0$'),
    ('b', 0.5 * t_c,    r'$t = 0.5\,t_c$'),
    ('c', 2.0 * t_c,    r'$t = 2\,t_c$'),
    ('d', 0.5 * t_r,    r'$t = t_r/2$'),
    ('e', 0.75 * t_r,   r'$t = 0.75\,t_r$'),
    ('f', t_r,          r'$t = t_r$'),
]
tlist = np.array(sorted(set(np.linspace(0, t_max, 2601)) | {s[1] for s in snapshots}))
print("Solving JC dynamics ...")
res = qt.mesolve(H_JC, psi0, tlist, e_ops=[sz], options={'store_states': True})
sigma_z = np.real(res.expect[0])
n_axis = np.arange(N_cav)
p0 = poisson_dist.pmf(n_axis, n_bar)


def pn_numeric(t):
    st = res.states[int(np.argmin(np.abs(tlist - t)))]
    return np.real(st.ptrace(0).diag())


def pn_exact(t):
    pm1 = np.concatenate([[0.0], p0[:-1]])
    return p0 * np.cos(g * np.sqrt(n_axis + 1) * t) ** 2 + pm1 * np.sin(g * np.sqrt(n_axis) * t) ** 2


err = max(np.abs(pn_numeric(t) - pn_exact(t)).max() for _, t, _ in snapshots)
print(f"numerical P(n,t) vs closed form at the snapshot times: max |error| = {err:.1e}")
assert err < 1e-5
pn_half = pn_numeric(0.5 * t_r)
parity = pn_half[0::2].sum() - pn_half[1::2].sum()
print(f"t_r/2: P(odd) = {pn_half[1::2].sum():.3f}, parity = {parity:+.3f}")


def bars(ax, pn, color=COLORS['blue']):
    ax.bar(n_axis[:N_SHOW], pn[:N_SHOW], width=0.8, color=color, lw=0)
    ax.plot(n_axis[:N_SHOW], p0[:N_SHOW], color='k', lw=0.6, ls='--', marker='.', ms=2,
            label='initial Poisson')
    ax.set_xlim(-0.8, N_SHOW - 0.2)


# ---------------------------------------------------------------------------
# Static figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(FULL_W, 3.5), sharex=True, sharey=True)
for ax, (letter, t, label) in zip(axes.flat, snapshots):
    bars(ax, pn_numeric(t))
    panel_label(ax, letter, label)
    ax.set_ylim(0, 0.26)
axes[0, 0].legend(loc='upper right', fontsize=7)
axes[1, 0].text(0.97, 0.93, rf'$\langle(-1)^{{\hat n}}\rangle = {parity:+.2f}$',
                transform=axes[1, 0].transAxes, ha='right', va='top', fontsize=7.5)
for ax in axes[-1]:
    ax.set_xlabel(r'$n$')
for ax in axes[:, 0]:
    ax.set_ylabel(r'$P(n)$')
save(fig, FIG_DIR / 'fig_photon_number_evolution', copy_to=(PAPER_FIG, SUPP_DIR))


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------
N_frames = 131                       # spacing t_r/100: t_r/2 and t_r are frames 50 and 100 exactly
frame_times = np.linspace(0, t_max, N_frames)
assert abs(frame_times[50] - t_r / 2) < 1e-9 and abs(frame_times[100] - t_r) < 1e-9
odd = n_axis[:N_SHOW] % 2 == 1
frames, durations = [], []
print("Rendering photon-number animation ...")
for t_now in frame_times:
    if t_now < 2.5 * t_c:
        phase = 'Rabi oscillations collapse'
    elif abs(t_now - t_r / 2) < 1.5:
        phase = 'cat state: parity comb'
    elif t_now < t_r / 2:
        phase = 'collapsed (two components separating)'
    elif abs(t_now - t_r) < 2.5:
        phase = 'revival'
    else:
        phase = 'components recombining' if t_now < t_r else 'second collapse'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.1), gridspec_kw={'width_ratios': [1.35, 1]})
    pn = pn_exact(t_now)
    bars(ax1, pn)
    ax1.bar(n_axis[:N_SHOW][odd], pn[:N_SHOW][odd], width=0.8, color=COLORS['orange'], lw=0)   # odd n on top
    ax1.set_ylim(0, 0.26)
    ax1.text(0.03, 0.95, rf'$\langle(-1)^{{\hat n}}\rangle = {pn[0::2].sum() - pn[1::2].sum():+.2f}$',
             transform=ax1.transAxes, ha='left', va='top', fontsize=8.5)
    ax1.text(0.03, 0.85, 'even $n$', color=COLORS['blue'], transform=ax1.transAxes, ha='left', va='top', fontsize=7.5)
    ax1.text(0.03, 0.78, 'odd $n$', color=COLORS['orange'], transform=ax1.transAxes, ha='left', va='top', fontsize=7.5)
    ax1.set_xlabel(r'$n$'); ax1.set_ylabel(r'$P(n,t)$')
    ax1.legend(loc='upper right', fontsize=7.5)
    ax1.set_title(rf'$gt = {t_now * g:5.1f}$: {phase}', fontsize=9, loc='left')

    ax2.plot(tlist * g, sigma_z, color=COLORS['blue'], lw=0.6)
    ax2.axhline(0, color='0.6', lw=0.4)
    time_guides(ax2, t_r * g, labels=True)
    ax2.plot(t_now * g, np.interp(t_now, tlist, sigma_z), 'o', color=COLORS['red'], ms=5,
             mec='k', mew=0.5, zorder=5)
    ax2.set_xlim(0, t_max * g); ax2.set_ylim(-1.1, 1.1)
    ax2.set_xlabel(r'$gt$'); ax2.set_ylabel(r'$\langle\sigma_z\rangle$')
    frames.append(fig_to_frame(fig, dpi=115))
    step = t_max / (N_frames - 1)
    durations.append(1500 if abs(t_now - t_r / 2) < step / 2 or abs(t_now - t_r) < step / 2
                     else (800 if t_now == 0 else 90))
save_gif(frames, ANIM_DIR / 'anim_photon_number.gif', durations)
