#!/usr/bin/env python3
"""
Bloch-vector trajectory of the reduced atomic state during JC dynamics.

The reduced atomic density matrix maps to a point r = (<sigma_x>, <sigma_y>,
<sigma_z>) inside the Bloch sphere.  For a real coherent amplitude and exact
resonance <sigma_x> = 0 at all times, so the whole trajectory lies in the
y-z plane and is drawn there.

During the collapse the Bloch vector spirals in toward the centre (atom
strongly entangled with the field).  Near the HALF revival time t = t_r/2 it
returns close to the surface, |r| = 0.96, on the sigma_y axis: the atom and
field approximately disentangle and the field is left in a nearly pure
Schrodinger-cat state (Gea-Banacloche, PRL 65, 3385 (1990)).  Around the full
revival t_r the atom is mixed again (|r| ~ 0.5, oscillating).

For a pure atom-field state the atomic purity is (1 + |r|^2)/2 and the
entanglement entropy is S = h((1 + |r|)/2), so the distance from the surface
is a direct picture of the entanglement.

Generates: figures/fig_bloch_sphere_trajectory (.png, .pdf), animations/anim_bloch_sphere.gif

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, save, panel_label, time_guides, fig_to_frame, save_gif,
                       FULL_W, COLORS, FIG_DIR, ANIM_DIR, PAPER_FIG, SUPP_DIR)

use_style()

g = 1.0
n_bar = 10
N_cav = 40
t_r = 2 * np.pi * np.sqrt(n_bar) / g
t_max = 1.3 * t_r
tlist = np.linspace(0, t_max, 2601)

a = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
H_JC = g * (a.dag() * sm + a * sm.dag())
psi0 = qt.tensor(qt.coherent(N_cav, np.sqrt(n_bar)), qt.basis(2, 0))
paulis = [qt.tensor(qt.qeye(N_cav), s) for s in (qt.sigmax(), qt.sigmay(), qt.sigmaz())]

print("Solving JC dynamics ...")
res = qt.mesolve(H_JC, psi0, tlist, e_ops=paulis)
rx, ry, rz = (np.real(e) for e in res.expect)
r = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
p = np.clip((1 + r) / 2, 1e-15, 1 - 1e-15)
entropy = -p * np.log2(p) - (1 - p) * np.log2(1 - p)      # S = h((1 + |r|)/2), pure global state
i_half = int(np.argmin(np.abs(tlist - t_r / 2)))
i_rev = int(np.argmin(np.abs(tlist - t_r)))
print(f"max |<sigma_x>| = {np.abs(rx).max():.1e} (trajectory is planar)")
print(f"|r|(t_r/2) = {r[i_half]:.3f}, S(t_r/2) = {entropy[i_half]:.3f} bit; "
      f"|r|(t_r) = {r[i_rev]:.3f}; min |r| = {r.min():.3f}")
assert np.abs(rx).max() < 1e-6


def disc(ax):
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(th), np.sin(th), color='k', lw=0.7)
    ax.axhline(0, color='0.8', lw=0.4); ax.axvline(0, color='0.8', lw=0.4)
    ax.set_aspect('equal')
    ax.set_xlim(-1.12, 1.12); ax.set_ylim(-1.12, 1.12)
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
    ax.set_xlabel(r'$\langle\sigma_y\rangle$'); ax.set_ylabel(r'$\langle\sigma_z\rangle$')


def trajectory(ax, upto=None, cmap='viridis'):
    n = len(tlist) if upto is None else upto + 1
    pts = np.column_stack([ry[:n], rz[:n]]).reshape(-1, 1, 2)
    lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1), cmap=cmap,
                        norm=plt.Normalize(0, t_max * g), lw=0.7, rasterized=True)
    lc.set_array(tlist[:n - 1] * g)
    ax.add_collection(lc)
    return lc


# ---------------------------------------------------------------------------
# Static figure
# ---------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(FULL_W, 2.3),
                                    gridspec_kw={'width_ratios': [1.3, 1, 1]})
disc(ax1)
lc = trajectory(ax1)
cb = fig.colorbar(lc, ax=ax1, shrink=0.85, aspect=22, pad=0.03, ticks=[0, 10, 20])
cb.set_label(r'$gt$', labelpad=1)
box = dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85)
for i, mk, lab, off, ha in ((0, 'o', r'$t = 0$', (0.10, -0.02), 'left'),
                            (i_half, '*', r'$t_r/2$', (-0.06, -0.22), 'right'),
                            (i_rev, '^', r'$t_r$', (-0.10, 0.10), 'right')):
    ax1.plot(ry[i], rz[i], mk, ms=7.5 if mk == '*' else 4.5, mfc='white', mec='k', mew=0.8, zorder=6)
    ax1.text(ry[i] + off[0], rz[i] + off[1], lab, fontsize=7.5, ha=ha, va='center', bbox=box, zorder=7)
panel_label(ax1, 'a', x=-0.02, y=1.12)

ax2.plot(tlist * g, r, color=COLORS['blue'], lw=0.9)
time_guides(ax2, t_r * g, labels=True)
ax2.set_xlabel(r'$gt$'); ax2.set_ylabel(r'$|\mathbf{r}|$')
ax2.set_xlim(0, t_max * g); ax2.set_ylim(0, 1.05)
panel_label(ax2, 'b', x=-0.02, y=1.12)

ax3.plot(tlist * g, entropy, color=COLORS['red'], lw=0.9)
ax3.axhline(1, color='0.6', lw=0.4)
time_guides(ax3, t_r * g, labels=True)
ax3.set_xlabel(r'$gt$'); ax3.set_ylabel(r'$S(\rho_{\rm atom})$ (bit)')
ax3.set_xlim(0, t_max * g); ax3.set_ylim(-0.03, 1.08)
panel_label(ax3, 'c', x=-0.02, y=1.12)
save(fig, FIG_DIR / 'fig_bloch_sphere_trajectory', copy_to=(PAPER_FIG, SUPP_DIR))


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------
N_frames = 131                       # every 20th point of tlist: spacing t_r/100, so t_r/2 and t_r are frames
frame_idx = np.unique(np.round(np.linspace(0, len(tlist) - 1, N_frames)).astype(int))
assert i_half in frame_idx and i_rev in frame_idx
frames, durations = [], []
print("Rendering Bloch-vector animation ...")
for idx in frame_idx:
    t_now = tlist[idx]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8.0, 3.0),
                                        gridspec_kw={'width_ratios': [1.05, 1, 1]})
    disc(ax1)
    ax1.plot(ry, rz, color='0.85', lw=0.4)
    if idx > 0:
        trajectory(ax1, upto=idx)
    ax1.plot([0, ry[idx]], [0, rz[idx]], color=COLORS['red'], lw=1.0)
    ax1.plot(ry[idx], rz[idx], 'o', color=COLORS['red'], ms=5, mec='k', mew=0.5, zorder=6)
    ax1.set_title(rf'Bloch vector (y–z plane), $gt = {t_now * g:.1f}$', fontsize=9)

    for ax, y, c, lab in ((ax2, r, COLORS['blue'], r'$|\mathbf{r}|$'),
                          (ax3, entropy, COLORS['red'], r'$S(\rho_{\rm atom})$ (bit)')):
        ax.plot(tlist * g, y, color='0.82', lw=0.6)
        ax.plot(tlist[:idx + 1] * g, y[:idx + 1], color=c, lw=1.0)
        ax.plot(t_now * g, y[idx], 'o', color=c, ms=4.5, mec='k', mew=0.5, zorder=5)
        time_guides(ax, t_r * g, labels=True)
        ax.set_xlim(0, t_max * g); ax.set_ylim(-0.03, 1.08)
        ax.set_xlabel(r'$gt$'); ax.set_ylabel(lab)
    frames.append(fig_to_frame(fig, dpi=115))
    step = t_max / N_frames
    durations.append(1500 if abs(t_now - t_r / 2) < step / 2 or abs(t_now - t_r) < step / 2
                     else (800 if idx == 0 else 90))
save_gif(frames, ANIM_DIR / 'anim_bloch_sphere.gif', durations)
