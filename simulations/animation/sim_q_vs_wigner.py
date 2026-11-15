#!/usr/bin/env python3
"""
Husimi Q function versus Wigner function of the cavity field during JC dynamics.

In the convention a = (x + ip)/sqrt(2) used throughout,
    Q(x, p) = <alpha| rho |alpha> / (2 pi),   alpha = (x + ip)/sqrt(2),
is normalised to one in (x, p), is never negative, and equals the Wigner
function smoothed by a vacuum-sized Gaussian.

At t = 0 both are Gaussians.  At the cat-state time t = t_r/2 Q shows the two
components but NO interference fringes, while W shows fringes with negative
values: the smoothing in Q averages over a fringe period (0.70) that is
smaller than the smoothing width, so the coherence between the components is
invisible in Q.  Wigner negativity is the more sensitive witness.

Generates: figures/fig_q_vs_wigner (.png, .pdf), animations/anim_q_vs_wigner.gif

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, save, panel_label, wigner_image, fig_to_frame, save_gif,
                       FULL_W, W_MAX, COLORS, FIG_DIR, ANIM_DIR, PAPER_FIG, SUPP_DIR)
from utils import wigner_negativity

use_style()

g = 1.0
n_bar = 10
N_cav = 40
t_r = 2 * np.pi * np.sqrt(n_bar) / g
t_c = 1.0 / g
t_max = 1.3 * t_r
Q_MAX = 1.0 / (2 * np.pi)          # peak of Q for a coherent state

a = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
H_JC = g * (a.dag() * sm + a * sm.dag())
psi0 = qt.tensor(qt.coherent(N_cav, np.sqrt(n_bar)), qt.basis(2, 0))

snapshots = [('a', 'e', 0.0, r'$t = 0$'), ('b', 'f', 2 * t_c, r'$t = 2\,t_c$'),
             ('c', 'g', 0.5 * t_r, r'$t = t_r/2$'), ('d', 'h', t_r, r'$t = t_r$')]
N_frames = 92                        # spacing t_r/70: t_r/2 and t_r are frames 35 and 70 exactly
frame_times = np.linspace(0, t_max, N_frames)
assert abs(frame_times[35] - t_r / 2) < 1e-9 and abs(frame_times[70] - t_r) < 1e-9
tlist = np.array(sorted(set(frame_times) | {s[2] for s in snapshots}))
print("Solving JC dynamics ...")
states = qt.mesolve(H_JC, psi0, tlist).states
field_at = lambda t: states[int(np.argmin(np.abs(tlist - t)))].ptrace(0)


def q_image(ax, Q, xvec):
    dx = xvec[1] - xvec[0]
    ext = [xvec[0] - dx / 2, xvec[-1] + dx / 2] * 2
    im = ax.imshow(np.clip(Q, 0, None), origin='lower', extent=ext, cmap='magma', vmin=0,
                   vmax=Q_MAX, interpolation='bilinear', rasterized=True)
    ax.set_aspect('equal')
    return im


# ---------------------------------------------------------------------------
# Static figure: W (top) and Q (bottom) at four times
# ---------------------------------------------------------------------------
xvec = np.linspace(-7, 7, 281)
fig, axes = plt.subplots(2, 4, figsize=(FULL_W, 3.3), sharex=True, sharey=True)
for col, (lw_, lq_, t, label) in enumerate(snapshots):
    rho = field_at(t)
    W = qt.wigner(rho, xvec, xvec)
    Q = qt.qfunc(rho, xvec, xvec)
    dx = xvec[1] - xvec[0]
    delta = wigner_negativity(W, xvec)
    print(f"  {label:18s} delta = {delta:6.3f}, min W = {W.min():+.3f}, min Q = {Q.min():+.1e}, "
          f"int Q = {Q.sum() * dx * dx:.4f}")
    im_w = wigner_image(axes[0, col], W, xvec)
    im_q = q_image(axes[1, col], Q, xvec)
    panel_label(axes[0, col], lw_, label)
    panel_label(axes[1, col], lq_)
    axes[0, col].text(0.05, 0.05, rf'$\delta = {max(delta, 0):.2f}$', transform=axes[0, col].transAxes,
                      fontsize=7.5, va='bottom')
    axes[1, col].set_xlabel(r'$x$')
    for ax in axes[:, col]:
        ax.set_xticks([-5, 0, 5]); ax.set_yticks([-5, 0, 5])
axes[0, 0].set_ylabel(r'$p$'); axes[1, 0].set_ylabel(r'$p$')
cbw = fig.colorbar(im_w, ax=axes[0], shrink=0.9, aspect=16, pad=0.02, ticks=[-0.3, 0, 0.3])
cbw.set_label(r'$W(x,p)$')
cbq = fig.colorbar(im_q, ax=axes[1], shrink=0.9, aspect=16, pad=0.02, ticks=[0, 0.05, 0.1, 0.15])
cbq.set_label(r'$Q(x,p)$')
save(fig, FIG_DIR / 'fig_q_vs_wigner', copy_to=(PAPER_FIG, SUPP_DIR))


# ---------------------------------------------------------------------------
# Animation: W, Q and their cuts along p = 0
# ---------------------------------------------------------------------------
xa = np.linspace(-7, 7, 141)
x_line = np.linspace(-7, 7, 701)
frames, durations = [], []
print("Rendering Q-vs-Wigner animation ...")
for t_now in frame_times:
    rho = field_at(t_now)
    W = qt.wigner(rho, xa, xa)
    Q = qt.qfunc(rho, xa, xa)
    W_cut = qt.wigner(rho, x_line, np.array([0.0]))[0]
    Q_cut = qt.qfunc(rho, x_line, np.array([0.0]))[0]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8.4, 3.0), gridspec_kw={'width_ratios': [1, 1, 1.1]})
    wigner_image(ax1, W, xa)
    delta_now = wigner_negativity(qt.wigner(rho, xvec, xvec), xvec)      # fine grid: the frame grid is too coarse for delta
    ax1.set_title(rf'Wigner $W$, $\delta = {max(delta_now, 0):.2f}$', fontsize=9)
    q_image(ax2, Q, xa)
    ax2.set_title(r'Husimi $Q \geq 0$', fontsize=9)
    for ax in (ax1, ax2):
        ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$p$')
        ax.set_xticks([-5, 0, 5]); ax.set_yticks([-5, 0, 5])
    ax3.axhline(0, color='0.6', lw=0.4)
    ax3.fill_between(x_line, W_cut, 0, where=W_cut < 0, color=COLORS['red'], alpha=0.3, lw=0)
    ax3.plot(x_line, W_cut, color=COLORS['blue'], lw=1.0, label=r'$W(x, 0)$')
    ax3.plot(x_line, Q_cut, color=COLORS['orange'], lw=1.2, ls=(0, (4, 2)), label=r'$Q(x, 0)$')
    ax3.set_xlim(-7, 7); ax3.set_ylim(-0.3, 0.34)
    ax3.set_xlabel(r'$x$')
    ax3.legend(loc='upper left', fontsize=7.5)
    ax3.set_title(rf'cut along $p = 0$, $gt = {t_now * g:.1f}$', fontsize=9)
    frames.append(fig_to_frame(fig, dpi=110))
    step = t_max / (N_frames - 1)
    durations.append(1800 if abs(t_now - t_r / 2) < step / 2 else
                     (1200 if abs(t_now - t_r) < step / 2 else (800 if t_now == 0 else 100)))
save_gif(frames, ANIM_DIR / 'anim_q_vs_wigner.gif', durations)
