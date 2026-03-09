# anim_entanglement.py
# Animated entropy S(rho_atom) + Wigner W(x,p) side by side.
# Shows entropy rising during collapse (cat forms) and dipping at revival.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: anim_entanglement.gif

import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap,
    mesolve, wigner, ptrace, entropy_vn, coherent, Qobj, ket2dm
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

# -- Solve --
n_frames = 120
tlist = np.linspace(0, 1.2 * t_revival, n_frames)

print(f"Solving JC dynamics: {n_frames} frames")
result = mesolve(H_JC, psi0, tlist, [], [sz], options={'store_states': True})
inversion = np.array(result.expect[0])

# Precompute entropy and purity
entropy = np.zeros(n_frames)
purity = np.zeros(n_frames)
for i in range(n_frames):
    rho_at = ptrace(result.states[i], 1)
    entropy[i] = entropy_vn(rho_at, 2)
    rho_f = ptrace(result.states[i], 0)
    purity[i] = (rho_f * rho_f).tr().real

# -- Render frames --
xvec = np.linspace(-7, 7, 150)
frame_dir = '/home/claude/ent_frames'
os.makedirs(frame_dir, exist_ok=True)
global_wlim = 0.30

print("Rendering frames...")
for i in range(n_frames):
    rho_field = ptrace(result.states[i], 0)
    W = wigner(rho_field, xvec, xvec)

    fig = plt.figure(figsize=(16, 5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

    # Left: Wigner
    ax1 = fig.add_subplot(gs[0])
    im = ax1.contourf(xvec, xvec, W, levels=80, cmap='RdBu_r',
                       vmin=-global_wlim, vmax=global_wlim, extend='both')
    ax1.set_xlabel(r'$x$', fontsize=12)
    ax1.set_ylabel(r'$p$', fontsize=12)
    ax1.set_aspect('equal')
    ax1.set_title(r'$W(x,p)$ of cavity field', fontsize=11)
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

    # Middle: inversion + entropy
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(tlist * g, inversion, color='#2166ac', linewidth=0.6, alpha=0.3)
    ax2.plot(tlist[:i+1] * g, inversion[:i+1], color='#2166ac', linewidth=1.2)
    ax2.plot(tlist[i] * g, inversion[i], 'o', color='#2166ac', markersize=8, zorder=5)
    ax2.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=12, color='#2166ac')
    ax2.set_xlabel(r'$gt$', fontsize=12)
    ax2.set_ylim(-1.1, 1.1)
    ax2.set_xlim(0, 1.2 * t_revival * g)
    ax2.set_title('Inversion', fontsize=11)
    ax2.axhline(0, color='gray', linewidth=0.5, linestyle='--')

    # Right: entropy + purity
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(tlist * g, entropy, 'r-', linewidth=0.6, alpha=0.3)
    ax3.plot(tlist[:i+1] * g, entropy[:i+1], 'r-', linewidth=1.5)
    ax3.plot(tlist[i] * g, entropy[i], 'ro', markersize=8, zorder=5)
    ax3.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=12, color='red')
    ax3.set_xlabel(r'$gt$', fontsize=12)
    ax3.set_ylim(-0.05, 1.15)
    ax3.set_xlim(0, 1.2 * t_revival * g)
    ax3.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
    ax3.set_title('Entanglement entropy', fontsize=11)

    # Purity on twin axis
    ax3b = ax3.twinx()
    ax3b.plot(tlist * g, purity, color='#1b7837', linewidth=0.6, alpha=0.3)
    ax3b.plot(tlist[:i+1] * g, purity[:i+1], color='#1b7837', linewidth=1.2)
    ax3b.plot(tlist[i] * g, purity[i], 'o', color='#1b7837', markersize=6, zorder=5)
    ax3b.set_ylabel(r'$\mathrm{Tr}(\rho_f^2)$', fontsize=11, color='#1b7837')
    ax3b.set_ylim(-0.05, 1.1)

    t_ratio = tlist[i] / t_revival
    fig.suptitle(
        rf'Jaynes-Cummings: $\bar{{n}} = {n_bar:.0f}$,  $t/t_r = {t_ratio:.2f}$,  '
        rf'$S = {entropy[i]:.2f}$ bits', fontsize=13, y=0.98
    )

    fig.savefig(f'{frame_dir}/frame_{i:04d}.png', dpi=90, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{n_frames} frames done")

# -- Stitch GIF --
print("Assembling GIF...")
frames = []
for i in range(n_frames):
    frames.append(imageio.imread(f'{frame_dir}/frame_{i:04d}.png'))

cat_idx = np.argmin(np.abs(tlist - 0.5 * t_revival))
rev_idx = np.argmin(np.abs(tlist - t_revival))

durations = [0.08] * n_frames
durations[0] = 1.0
durations[cat_idx] = 1.5
durations[rev_idx] = 1.0

imageio.mimsave('/home/claude/anim_entanglement.gif', frames,
                duration=durations, loop=0)
print("Saved: anim_entanglement.gif")

import shutil
shutil.rmtree(frame_dir)
