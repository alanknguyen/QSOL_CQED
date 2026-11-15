# anim_wigner_evolution.py
# Animated Wigner function W(x,p) of the cavity field during JC dynamics:
# the coherent state splits into two counter-rotating components, which form a
# Schrodinger cat at t_r/2 and meet again at the revival t_r.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585 and EC 777, taught by Prof. Luca Dal Negro
#
# Generates: animations/anim_wigner_evolution.gif

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, tensor, destroy, sigmaz, sigmam, qeye, mesolve, wigner, coherent

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, wigner_image, time_guides, fig_to_frame, save_gif,
                       COLORS, ANIM_DIR)
from utils import wigner_negativity

use_style()

N_cav = 35
g = 1.0
n_bar = 10.0
t_revival = 2 * np.pi * np.sqrt(n_bar) / g
t_max = 1.2 * t_revival

a = tensor(destroy(N_cav), qeye(2))
sm = tensor(qeye(N_cav), sigmam())
sz = tensor(qeye(N_cav), sigmaz())
H_JC = g * (a.dag() * sm + a * sm.dag())
psi0 = tensor(coherent(N_cav, np.sqrt(n_bar)), basis(2, 0))

# The inversion is solved on a fine grid (the Rabi period is ~0.95/g); the
# frames are a subset of it.  Using the frame times alone aliases the curve.
# 121 frames over 1.2 t_r: the frame spacing is t_r/100, so t_r/2 and t_r are
# frames 50 and 100 exactly and the numbers shown at the pauses are those of
# the paper (delta = 0.85 at t_r/2).
n_frames = 121
STRIDE = 20
tlist = np.linspace(0, t_max, (n_frames - 1) * STRIDE + 1)
print(f"Solving JC dynamics on {len(tlist)} time points, {n_frames} frames")
result = mesolve(H_JC, psi0, tlist, e_ops=[sz], options={'store_states': True})
inversion = np.array(result.expect[0])

xvec = np.linspace(-7, 7, 281)          # 14 points per fringe, as in the figures
frames, durations = [], []
for k in range(n_frames):
    i = k * STRIDE
    W = wigner(result.states[i].ptrace(0), xvec, xvec)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.3), gridspec_kw={'width_ratios': [1.15, 1]})
    im = wigner_image(ax1, W, xvec)
    ax1.set_xlabel(r'$x$'); ax1.set_ylabel(r'$p$')
    ax1.set_xticks([-5, 0, 5]); ax1.set_yticks([-5, 0, 5])
    fig.colorbar(im, ax=ax1, shrink=0.9, aspect=22, pad=0.03, ticks=[-0.3, -0.15, 0, 0.15, 0.3],
                 label=r'$W(x,p)$')
    ax1.set_title(rf'cavity field, $\delta = {max(wigner_negativity(W, xvec), 0):.2f}$', fontsize=9)

    ax2.plot(tlist * g, inversion, color='0.8', lw=0.6)
    ax2.plot(tlist[:i + 1] * g, inversion[:i + 1], color=COLORS['blue'], lw=0.8)
    ax2.plot(tlist[i] * g, inversion[i], 'o', color=COLORS['red'], ms=5, mec='k', mew=0.5, zorder=5)
    ax2.axhline(0, color='0.6', lw=0.4)
    time_guides(ax2, t_revival * g, labels=True)
    ax2.set_xlabel(r'$gt$'); ax2.set_ylabel(r'$\langle \sigma_z \rangle$')
    ax2.set_xlim(0, t_max * g); ax2.set_ylim(-1.1, 1.1)
    ax2.set_title(rf'atomic inversion, $t/t_r = {tlist[i] / t_revival:.2f}$', fontsize=9, pad=16)
    frames.append(fig_to_frame(fig, dpi=115))
    durations.append(80)
    if (k + 1) % 30 == 0:
        print(f"  {k + 1}/{n_frames} frames")

# pause on the initial state, the cat state (t_r/2) and the revival (t_r)
frame_t = tlist[::STRIDE]
durations[0] = 1000
assert abs(frame_t[50] - 0.5 * t_revival) < 1e-9 and abs(frame_t[100] - t_revival) < 1e-9
durations[int(np.argmin(np.abs(frame_t - 0.5 * t_revival)))] = 1800
durations[int(np.argmin(np.abs(frame_t - t_revival)))] = 1200
save_gif(frames, ANIM_DIR / 'anim_wigner_evolution.gif', durations)
