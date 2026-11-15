# anim_entanglement.py
# Animated Wigner function, inversion, entanglement entropy and field purity.
# The entropy RISES during the collapse (the reduced field is then a two-branch
# mixture, purity ~ 0.5) and reaches its MINIMUM at t_r/2, where the atom and
# field disentangle and the field is a nearly pure cat (purity 0.96).  Toward
# the revival t_r the entropy grows again.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585 and EC 777, taught by Prof. Luca Dal Negro
#
# Generates: animations/anim_entanglement.gif

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qutip import (basis, tensor, destroy, sigmaz, sigmam, qeye, mesolve, wigner,
                   entropy_vn, coherent)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, wigner_image, time_guides, fig_to_frame, save_gif,
                       COLORS, ANIM_DIR)

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

n_frames = 121                   # frame spacing t_r/100: t_r/2 and t_r are frames 50 and 100 exactly
STRIDE = 20                      # fine time grid for the curves; every STRIDE-th point is a frame
tlist = np.linspace(0, t_max, (n_frames - 1) * STRIDE + 1)
print(f"Solving JC dynamics on {len(tlist)} time points, {n_frames} frames")
result = mesolve(H_JC, psi0, tlist, e_ops=[sz], options={'store_states': True})
inversion = np.array(result.expect[0])
entropy = np.array([entropy_vn(st.ptrace(1), 2) for st in result.states])
purity = np.array([(st.ptrace(0) ** 2).tr().real for st in result.states])

xvec = np.linspace(-7, 7, 161)
frames, durations = [], []
for k in range(n_frames):
    i = k * STRIDE
    W = wigner(result.states[i].ptrace(0), xvec, xvec)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.2, 3.1),
                                        gridspec_kw={'width_ratios': [1.05, 1, 1]})
    im = wigner_image(ax1, W, xvec)
    ax1.set_xlabel(r'$x$'); ax1.set_ylabel(r'$p$')
    ax1.set_xticks([-5, 0, 5]); ax1.set_yticks([-5, 0, 5])
    fig.colorbar(im, ax=ax1, shrink=0.85, aspect=22, pad=0.03, ticks=[-0.3, 0, 0.3])
    ax1.set_title(rf'$W(x,p)$, $t/t_r = {tlist[i] / t_revival:.2f}$', fontsize=9)

    def trace(ax, y, color, ls='-'):
        ax.plot(tlist * g, y, color='0.82', lw=0.5, ls=ls)
        ax.plot(tlist[:i + 1] * g, y[:i + 1], color=color, lw=0.9, ls=ls)
        ax.plot(tlist[i] * g, y[i], 'o', color=color, ms=4.5, mec='k', mew=0.5, zorder=5)

    trace(ax2, inversion, COLORS['blue'])
    ax2.axhline(0, color='0.6', lw=0.4)
    ax2.set_ylim(-1.1, 1.1); ax2.set_ylabel(r'$\langle \sigma_z \rangle$')
    trace(ax3, entropy, COLORS['red'])
    trace(ax3, purity, COLORS['green'], ls=(0, (4, 2)))
    ax3.plot([], [], color=COLORS['red'], lw=0.9, label=rf'$S = {entropy[i]:.2f}$ bit')
    ax3.plot([], [], color=COLORS['green'], lw=0.9, ls=(0, (4, 2)), label=rf'$\mathcal{{P}} = {purity[i]:.2f}$')
    ax3.legend(loc='lower right', fontsize=7.5, frameon=True, framealpha=0.9, edgecolor='none')
    ax3.set_ylim(-0.03, 1.08)
    ax3.set_ylabel(r'entropy $S$, field purity $\mathcal{P}$')
    for ax in (ax2, ax3):
        time_guides(ax, t_revival * g, labels=True)
        ax.set_xlim(0, t_max * g); ax.set_xlabel(r'$gt$')
    frames.append(fig_to_frame(fig, dpi=105))
    durations.append(80)
    if (k + 1) % 30 == 0:
        print(f"  {k + 1}/{n_frames} frames")

frame_t = tlist[::STRIDE]
durations[0] = 1000
assert abs(frame_t[50] - 0.5 * t_revival) < 1e-9 and abs(frame_t[100] - t_revival) < 1e-9
durations[int(np.argmin(np.abs(frame_t - 0.5 * t_revival)))] = 1800
durations[int(np.argmin(np.abs(frame_t - t_revival)))] = 1200
save_gif(frames, ANIM_DIR / 'anim_entanglement.gif', durations)
