#!/usr/bin/env python3
"""
Animation: Cat-state survival as κ/g increases.

For fixed n̄ = 10, sweeps κ/g from 0 to 0.12, showing:
  - Left: Wigner function W(x,p) at t = t_r/2
  - Right top: δ(κ/g) building up as a curve
  - Right bottom: Cat fidelity F(κ/g) building up

Produces: anim_phase_diagram.gif
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    destroy, sigmap, sigmam, tensor, qeye,
    basis, coherent, mesolve, wigner, ket2dm
)
import imageio
import os, shutil

fig_dir = '/home/claude/animations'
os.makedirs(fig_dir, exist_ok=True)
frame_dir = '/home/claude/_frames_phase'
os.makedirs(frame_dir, exist_ok=True)

# ── Parameters ───────────────────────────────────────────────
g = 1.0
n_bar = 10
alpha = np.sqrt(n_bar)
N_cav = 50
t_r = 2 * np.pi * np.sqrt(n_bar) / g
t_cat = t_r / 2.0

xvec = np.linspace(-7, 7, 150)

kappa_values = np.linspace(0.0, 0.12, 40)

# ── Precompute all frames ────────────────────────────────────
print("Precomputing Wigner functions...")
wigner_frames = []
delta_vals = []
fid_vals = []

for i, kappa_g in enumerate(kappa_values):
    kappa = kappa_g * g

    a = tensor(destroy(N_cav), qeye(2))
    sm = tensor(qeye(N_cav), sigmam())
    sp = tensor(qeye(N_cav), sigmap())

    H = g * (a.dag() * sm + a * sp)
    c_ops = [np.sqrt(kappa) * a] if kappa > 0 else []

    psi0 = tensor(coherent(N_cav, alpha), basis(2, 0))
    result = mesolve(H, psi0, [0, t_cat], c_ops, [])
    rho_final = result.states[-1]
    if rho_final.isket:
        rho_final = ket2dm(rho_final)

    rho_field = rho_final.ptrace(0)

    W = wigner(rho_field, xvec, xvec)
    wigner_frames.append(W)

    dx = xvec[1] - xvec[0]
    delta = max(np.sum(np.abs(W)) * dx * dx - 1.0, 0.0)
    delta_vals.append(delta)

    purity = np.real((rho_field * rho_field).tr())
    fid_vals.append(purity)

    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{len(kappa_values)} (κ/g = {kappa_g:.3f}, δ = {delta:.3f}, P = {purity:.3f})")

delta_vals = np.array(delta_vals)
fid_vals = np.array(fid_vals)

# ── Render frames ────────────────────────────────────────────
print("Rendering frames...")
wmax = max(np.max(np.abs(W)) for W in wigner_frames) * 0.9

for i in range(len(kappa_values)):
    fig = plt.figure(figsize=(14, 5.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], hspace=0.35, wspace=0.3)

    # Left: Wigner function
    ax_w = fig.add_subplot(gs[:, 0])
    im = ax_w.pcolormesh(xvec, xvec, wigner_frames[i], cmap='RdBu_r',
                         vmin=-wmax, vmax=wmax, shading='auto')
    ax_w.set_xlabel('$x$', fontsize=12)
    ax_w.set_ylabel('$p$', fontsize=12)
    ax_w.set_aspect('equal')
    ax_w.set_title(rf'$W(x,p)$ at $t = t_r/2$    '
                   rf'$\kappa/g = {kappa_values[i]:.3f}$',
                   fontsize=13)
    plt.colorbar(im, ax=ax_w, shrink=0.8)

    # Right top: δ vs κ/g
    ax_d = fig.add_subplot(gs[0, 1])
    ax_d.plot(kappa_values[:i+1], delta_vals[:i+1], 'b-o', markersize=4, linewidth=2)
    ax_d.plot(kappa_values[i], delta_vals[i], 'ro', markersize=10, zorder=5)
    ax_d.set_xlim(-0.005, 0.125)
    ax_d.set_ylim(-0.02, max(delta_vals) * 1.1)
    ax_d.axhline(0.05, color='grey', ls=':', lw=1)
    ax_d.set_ylabel(r'$\delta$', fontsize=13)
    ax_d.set_title(r'Wigner Negativity', fontsize=12)
    ax_d.text(0.08, 0.06, r'$\delta = 0.05$', fontsize=9, color='grey')

    # Right bottom: Fidelity vs κ/g
    ax_f = fig.add_subplot(gs[1, 1])
    ax_f.plot(kappa_values[:i+1], fid_vals[:i+1], 'g-o', markersize=4, linewidth=2)
    ax_f.plot(kappa_values[i], fid_vals[i], 'ro', markersize=10, zorder=5)
    ax_f.set_xlim(-0.005, 0.125)
    ax_f.set_ylim(-0.02, 1.05)
    ax_f.set_xlabel(r'$\kappa / g$', fontsize=13)
    ax_f.set_ylabel(r'$\mathrm{Tr}[\rho^2]$', fontsize=13)
    ax_f.set_title(r'Field-State Purity', fontsize=12)

    fig.suptitle(rf'Cat-State Decoherence — $\bar{{n}} = {n_bar}$',
                 fontsize=14, y=1.0)

    fig.savefig(f'{frame_dir}/frame_{i:04d}.png', dpi=120, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)

# ── Assemble GIF ─────────────────────────────────────────────
print("Assembling GIF...")
frames = []
for i in range(len(kappa_values)):
    frames.append(imageio.v2.imread(f'{frame_dir}/frame_{i:04d}.png'))

durations = [0.15] * len(kappa_values)
durations[0] = 1.5   # pause on κ=0 (ideal cat)
durations[-1] = 1.5   # pause on final

outpath = f'{fig_dir}/anim_phase_diagram.gif'
imageio.mimsave(outpath, frames, duration=durations, loop=0)
print(f"Saved: {outpath}")
shutil.rmtree(frame_dir)
