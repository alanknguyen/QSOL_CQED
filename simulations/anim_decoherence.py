# anim_decoherence.py
# Animated destruction of the Schrodinger cat state by cavity loss.
# Sweeps kappa/g from 0 to 0.15 at fixed t = t_r/2.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: anim_decoherence.gif

import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
from qutip import (
    basis, tensor, destroy, sigmam, sigmap,
    mesolve, wigner, ptrace, coherent, Qobj
)

# -- Parameters --
N_cav = 35
g = 1.0
alpha = np.sqrt(10.0)
n_bar = np.abs(alpha)**2
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# -- Operators --
a  = tensor(destroy(N_cav), Qobj(np.eye(2)))
sm = tensor(Qobj(np.eye(N_cav)), sigmam())
sp = tensor(Qobj(np.eye(N_cav)), sigmap())

H_JC = g * (a.dag() * sm + a * sp)
psi0 = tensor(coherent(N_cav, alpha), basis(2, 0))

# -- Kappa sweep --
n_frames = 60
kappa_values = np.linspace(0, 0.15, n_frames)
tlist_short = np.linspace(0, 0.5 * t_revival, 400)
xvec = np.linspace(-7, 7, 180)

frame_dir = '/home/claude/deco_frames'
os.makedirs(frame_dir, exist_ok=True)

global_wlim = 0.25
negativities = []

print(f"Rendering {n_frames} decoherence frames...")
for i, kap in enumerate(kappa_values):
    c_ops = [np.sqrt(kap) * a] if kap > 0 else []
    res = mesolve(H_JC, psi0, tlist_short, c_ops, [],
                  options={'store_states': True})
    rho_f = ptrace(res.states[-1], 0)
    W = wigner(rho_f, xvec, xvec)

    dx = xvec[1] - xvec[0]
    neg_vol = np.sum(np.abs(W)) * dx**2 - 1
    purity = (rho_f * rho_f).tr().real
    negativities.append(neg_vol)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5),
                                     gridspec_kw={'width_ratios': [1.3, 1]})

    # Left: Wigner
    im = ax1.contourf(xvec, xvec, W, levels=80, cmap='RdBu_r',
                       vmin=-global_wlim, vmax=global_wlim, extend='both')
    ax1.set_xlabel(r'$x$', fontsize=14)
    ax1.set_ylabel(r'$p$', fontsize=14)
    ax1.set_aspect('equal')
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label=r'$W(x,p)$')
    ax1.set_title(rf'$\kappa/g = {kap:.3f}$    $\delta = {neg_vol:.3f}$    '
                  rf'$\mathcal{{P}} = {purity:.3f}$', fontsize=11)

    # Right: negativity vs kappa (running trace)
    ax2.plot(kappa_values[:i+1], negativities, 'r-', linewidth=2)
    ax2.plot(kap, neg_vol, 'ro', markersize=10, zorder=5)
    ax2.set_xlim(0, 0.15)
    ax2.set_ylim(-0.02, 0.5)
    ax2.set_xlabel(r'$\kappa / g$', fontsize=14)
    ax2.set_ylabel(r'Wigner negativity $\delta$', fontsize=14)
    ax2.set_title('Non-classicality vs. cavity loss', fontsize=12)
    ax2.axhline(0, color='gray', linewidth=0.5, linestyle='--')

    fig.suptitle(
        rf'Decoherence of Schrodinger Cat State at $t = t_r/2$  '
        rf'($\bar{{n}} = {n_bar:.0f}$)',
        fontsize=13, y=0.98
    )

    fig.savefig(f'{frame_dir}/frame_{i:04d}.png', dpi=100, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

    if (i + 1) % 15 == 0:
        print(f"  {i+1}/{n_frames} frames done")

# -- Stitch GIF --
print("Assembling GIF...")
frames = []
for i in range(n_frames):
    frames.append(imageio.imread(f'{frame_dir}/frame_{i:04d}.png'))

durations = [0.12] * n_frames
durations[0] = 1.5   # pause on pure cat state
durations[-1] = 1.5  # pause on fully decohered

imageio.mimsave('/home/claude/anim_decoherence.gif', frames,
                duration=durations, loop=0)
print("Saved: anim_decoherence.gif")

import shutil
shutil.rmtree(frame_dir)
