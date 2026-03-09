# anim_wigner_evolution.py
# Animated Wigner function W(x,p) of the cavity field during JC dynamics.
# Shows coherent state splitting into a Schrodinger cat state.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: anim_wigner_evolution.gif

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import imageio
import os
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap,
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
sz = tensor(Qobj(np.eye(N_cav)), sigmaz())

H_JC = g * (a.dag() * sm + a * sp)
psi0 = tensor(coherent(N_cav, alpha), basis(2, 0))

# -- Solve for full time range --
n_frames = 120
tlist = np.linspace(0, 1.2 * t_revival, n_frames)

print(f"Solving JC dynamics: {n_frames} frames, t_max = {tlist[-1]:.1f}")
result = mesolve(H_JC, psi0, tlist, [], [sz], options={'store_states': True})
inversion = np.array(result.expect[0])

# -- Generate frames --
xvec = np.linspace(-7, 7, 180)
frame_dir = '/home/claude/wigner_frames'
os.makedirs(frame_dir, exist_ok=True)

# Fixed colorscale across all frames
global_wlim = 0.30

print("Rendering frames...")
for i in range(n_frames):
    rho_field = ptrace(result.states[i], 0)
    W = wigner(rho_field, xvec, xvec)

    fig = plt.figure(figsize=(14, 5.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.3)

    # Left: Wigner function
    ax1 = fig.add_subplot(gs[0])
    im = ax1.contourf(xvec, xvec, W, levels=80, cmap='RdBu_r',
                       vmin=-global_wlim, vmax=global_wlim,
                       extend='both')
    ax1.set_xlabel(r'$x$', fontsize=14)
    ax1.set_ylabel(r'$p$', fontsize=14)
    ax1.set_aspect('equal')
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label=r'$W(x,p)$')

    # Wigner negativity
    dx = xvec[1] - xvec[0]
    neg_vol = np.sum(np.abs(W)) * dx**2 - 1
    ax1.set_title(rf'Cavity Field Wigner Function    $\delta = {neg_vol:.3f}$',
                  fontsize=12)

    # Right: inversion time trace with moving marker
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(tlist * g, inversion, 'b-', linewidth=0.8, alpha=0.4)
    ax2.plot(tlist[:i+1] * g, inversion[:i+1], 'b-', linewidth=1.5)
    ax2.plot(tlist[i] * g, inversion[i], 'ro', markersize=10, zorder=5)
    ax2.set_xlabel(r'$gt$', fontsize=14)
    ax2.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=14)
    ax2.set_title('Atomic Inversion', fontsize=12)
    ax2.set_xlim(0, 1.2 * t_revival * g)
    ax2.set_ylim(-1.1, 1.1)
    ax2.axhline(0, color='gray', linewidth=0.5, linestyle='--')

    # Time label
    t_ratio = tlist[i] / t_revival
    fig.suptitle(
        rf'Jaynes-Cummings Dynamics: $\bar{{n}} = {n_bar:.0f}$,  $t/t_r = {t_ratio:.2f}$',
        fontsize=14, y=0.98
    )

    fig.savefig(f'{frame_dir}/frame_{i:04d}.png', dpi=100, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{n_frames} frames done")

# -- Stitch into GIF --
print("Assembling GIF...")
frames = []
for i in range(n_frames):
    frames.append(imageio.imread(f'{frame_dir}/frame_{i:04d}.png'))

# Pause on key moments: t=0, cat state (t_r/2), revival (t_r)
cat_idx = np.argmin(np.abs(tlist - 0.5 * t_revival))
rev_idx = np.argmin(np.abs(tlist - t_revival))

durations = [0.08] * n_frames
durations[0] = 1.0          # pause on initial state
durations[cat_idx] = 1.5    # pause on cat state
durations[rev_idx] = 1.0    # pause on revival

imageio.mimsave('/home/claude/anim_wigner_evolution.gif', frames,
                duration=durations, loop=0)
print("Saved: anim_wigner_evolution.gif")

# Cleanup frames
import shutil
shutil.rmtree(frame_dir)
