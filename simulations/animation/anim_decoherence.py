# anim_decoherence.py
# Destruction of the Schrodinger cat state by cavity loss: the reduced field
# at t = t_r/2 as kappa/g is swept from 0 to 0.15.
#
# The fringe visibility follows  exp(-pi kappa n_bar^{3/2} / g):  the fringe
# decoherence rate kappa |Delta alpha|^2 / 2 = 2 kappa n_bar sin^2(g t / 2 sqrt(n_bar))
# integrated over the time the two components take to separate.  That estimate
# is drawn as a dashed line next to the simulated Wigner negativity.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585 and EC 777, taught by Prof. Luca Dal Negro
#
# Generates: animations/anim_decoherence.gif

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, tensor, destroy, sigmam, qeye, mesolve, wigner, coherent

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import use_style, wigner_image, fig_to_frame, save_gif, COLORS, ANIM_DIR
from utils import wigner_negativity

use_style()

N_cav = 35
g = 1.0
n_bar = 10.0
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

a = tensor(destroy(N_cav), qeye(2))
sm = tensor(qeye(N_cav), sigmam())
H_JC = g * (a.dag() * sm + a * sm.dag())
psi0 = tensor(coherent(N_cav, np.sqrt(n_bar)), basis(2, 0))

n_frames = 60
# Quadratic spacing: the fringes are gone by kappa/g ~ 0.03, so a uniform grid would spend
# five frames on everything that happens and fifty on the aftermath.
kappa_values = 0.15 * np.linspace(0, 1, n_frames) ** 2
xvec = np.linspace(-7, 7, 281)

print(f"Solving {n_frames} dissipative runs ...")
W_frames, deltas, purities = [], [], []
for kap in kappa_values:
    c_ops = [np.sqrt(kap) * a] if kap > 0 else []
    rho_f = mesolve(H_JC, psi0, [0, 0.5 * t_revival], c_ops=c_ops).states[-1].ptrace(0)
    W = wigner(rho_f, xvec, xvec)
    W_frames.append(W)
    deltas.append(max(wigner_negativity(W, xvec), 1e-6))
    purities.append((rho_f * rho_f).tr().real)
deltas = np.array(deltas)
model = deltas[0] * np.exp(-np.pi * kappa_values * n_bar ** 1.5 / g)

frames = []
for i, kap in enumerate(kappa_values):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.3), gridspec_kw={'width_ratios': [1.15, 1]})
    im = wigner_image(ax1, W_frames[i], xvec)
    ax1.set_xlabel(r'$x$'); ax1.set_ylabel(r'$p$')
    ax1.set_xticks([-5, 0, 5]); ax1.set_yticks([-5, 0, 5])
    fig.colorbar(im, ax=ax1, shrink=0.9, aspect=22, pad=0.03, ticks=[-0.3, -0.15, 0, 0.15, 0.3],
                 label=r'$W(x,p)$')
    ax1.set_title(rf'$t = t_r/2$, $\kappa/g = {kap:.4f}$, $\mathcal{{P}} = {purities[i]:.2f}$', fontsize=9)

    ax2.semilogy(kappa_values, model, color='0.4', lw=0.8, ls='--',
                 label=r'$\delta_0\,e^{-\pi\kappa\bar{n}^{3/2}/g}$')
    ax2.semilogy(kappa_values[:i + 1], deltas[:i + 1], color=COLORS['red'], lw=1.1, label='simulation')
    ax2.semilogy(kap, deltas[i], 'o', color=COLORS['red'], ms=5, mec='k', mew=0.5, zorder=5)
    ax2.set_xlim(0, 0.15); ax2.set_ylim(3e-4, 1.5)
    ax2.set_xlabel(r'$\kappa/g$'); ax2.set_ylabel(r'Wigner negativity $\delta$')
    ax2.legend(loc='upper right', fontsize=7.5)
    ax2.set_title(rf'$\delta = {deltas[i]:.3f}$', fontsize=9)
    frames.append(fig_to_frame(fig, dpi=115))

durations = [130] * n_frames
durations[0] = 1800      # pause on the pure cat
durations[-1] = 1800     # pause on the fully decohered mixture
save_gif(frames, ANIM_DIR / 'anim_decoherence.gif', durations)
