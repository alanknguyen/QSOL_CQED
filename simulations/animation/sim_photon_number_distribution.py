#!/usr/bin/env python3
"""
Photon number distribution P(n,t) during Jaynes-Cummings dynamics.

At t=0 the distribution is Poissonian (coherent state, n_bar=10).
During collapse the distribution splits into two peaks — the two
coherent components of the Schrödinger cat state separating in
photon-number space.  At revival it reforms into a single peak.

This is the Fock-space complement to the Wigner function evolution
(Fig. 2 in the paper): the photon number distribution literally
splits, which is WHY the Wigner function develops two lobes.

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import qutip as qt
from PIL import Image
import io
import os

# ── Parameters ──────────────────────────────────────────────────────
g = 1.0                      # vacuum Rabi coupling (sets units)
n_bar = 10                    # mean photon number
alpha = np.sqrt(n_bar)        # coherent state amplitude
N_cav = 40                    # Fock space truncation
t_r = 2 * np.pi * np.sqrt(n_bar) / g   # revival time
t_c = np.pi / g               # collapse time

N_frames = 100                # animation frames
t_max = 1.3 * t_r
tlist = np.linspace(0, t_max, 800)

# ── Build JC Hamiltonian (resonant) ─────────────────────────────────
a  = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
sp = qt.tensor(qt.qeye(N_cav), qt.sigmap())
sz = qt.tensor(qt.qeye(N_cav), qt.sigmaz())

H_JC = g * (a.dag() * sm + a * sp)

# ── Initial state: |e> ⊗ |α> ───────────────────────────────────────
psi0 = qt.tensor(qt.coherent(N_cav, alpha), qt.basis(2, 0))

# ── Time evolution (unitary) ────────────────────────────────────────
print("Running JC time evolution...")
result = qt.mesolve(H_JC, psi0, tlist, e_ops=[])
# Compute sigma_z from states
sigma_z = np.array([qt.expect(sz, result.states[i]) for i in range(len(tlist))])

# ── Extract P(n,t) and inversion at snapshot times ──────────────────
frame_times = np.linspace(0, t_max, N_frames)
frame_indices = [np.argmin(np.abs(tlist - ft)) for ft in frame_times]

# Fock state projectors for the field (traced over atom)
def get_pn(state, N):
    """Get photon number distribution from joint state."""
    rho_field = state.ptrace(0)
    return np.real(np.array([rho_field[n, n] for n in range(N)]))

print("Computing photon number distributions...")
pn_snapshots = []
for idx in frame_indices:
    psi_t = result.states[idx]
    pn = get_pn(psi_t, N_cav)
    pn_snapshots.append(pn)

# ── Static figure: 6 key snapshots ──────────────────────────────────
snapshot_labels = [
    (0,           r'$t = 0$' + '\n(coherent state)'),
    (t_c * 0.5,   r'$t = 0.5\,t_c$' + '\n(early Rabi)'),
    (t_c * 2,     r'$t = 2\,t_c$' + '\n(collapse onset)'),
    (t_r / 2,     r'$t = t_r/2$' + '\n(cat state)'),
    (t_r * 0.75,  r'$t = 0.75\,t_r$' + '\n(pre-revival)'),
    (t_r,         r'$t = t_r$' + '\n(first revival)'),
]

fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
fig.suptitle(
    r'Photon Number Distribution $P(n,t)$ — Jaynes–Cummings Model'
    + f'\n' + r'($\bar{{n}} = {0}$, $\Delta = 0$)'.format(n_bar),
    fontsize=14, fontweight='bold'
)

n_axis = np.arange(N_cav)
from math import factorial as _fact
poisson = np.array([np.exp(-n_bar) * n_bar**k / _fact(k) for k in range(N_cav)])

for i, (t_snap, label) in enumerate(snapshot_labels):
    ax = axes.flat[i]
    idx = np.argmin(np.abs(tlist - t_snap))
    psi_t = result.states[idx]
    pn = get_pn(psi_t, N_cav)

    ax.bar(n_axis[:30], pn[:30], color='steelblue', alpha=0.8, width=0.8,
           edgecolor='navy', linewidth=0.3)
    if i == 0:
        ax.plot(n_axis[:30], poisson[:30], 'r--', lw=1.5, label='Poisson')
        ax.legend(fontsize=8)
    ax.set_title(label, fontsize=10)
    ax.set_xlim(-0.5, 28)
    ax.set_xlabel(r'$n$', fontsize=10)
    if i % 3 == 0:
        ax.set_ylabel(r'$P(n)$', fontsize=10)
    ax.tick_params(labelsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig('/home/claude/figures/fig_photon_number_evolution.png', dpi=200,
            bbox_inches='tight')
fig.savefig('/home/claude/figures/fig_photon_number_evolution.pdf',
            bbox_inches='tight')
print("Static figure saved.")
plt.close()

# ── GIF animation ───────────────────────────────────────────────────
print("Generating animation frames...")
frames = []
cmap = plt.cm.coolwarm
norm = Normalize(vmin=0, vmax=t_max)

for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]
    pn = pn_snapshots[fi]

    fig = plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.4, 1], wspace=0.3)

    # ─ Left: P(n,t) bar chart
    ax1 = fig.add_subplot(gs[0])
    colors = ['crimson' if (t_r/2 - 1.0) < t_now < (t_r/2 + 1.0) else
              'steelblue' for _ in range(30)]
    ax1.bar(n_axis[:30], pn[:30], color='steelblue', alpha=0.85, width=0.8,
            edgecolor='navy', linewidth=0.3)
    # Overlay initial Poisson
    ax1.plot(n_axis[:30], poisson[:30], 'r--', lw=1.2, alpha=0.4,
             label=r'$P_{\rm Poisson}$')
    ax1.set_xlim(-0.5, 28)
    ax1.set_ylim(0, 0.22)
    ax1.set_xlabel(r'Photon number $n$', fontsize=11)
    ax1.set_ylabel(r'$P(n, t)$', fontsize=11)
    ax1.set_title(r'Photon Number Distribution', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9)

    # Time annotation
    phase = ''
    if t_now < t_c:
        phase = 'Rabi oscillations'
    elif t_now < t_r * 0.35:
        phase = 'Collapse → splitting'
    elif abs(t_now - t_r/2) < 1.5:
        phase = 'Cat state (two peaks)'
    elif t_now < t_r * 0.85:
        phase = 'Pre-revival'
    elif abs(t_now - t_r) < 1.5:
        phase = 'Revival (single peak)'
    else:
        phase = 'Post-revival'

    ax1.text(0.02, 0.95, f'$gt = {t_now:.1f}$\n{phase}',
             transform=ax1.transAxes, fontsize=10, va='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

    # ─ Right: inversion with time marker
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(tlist, sigma_z, 'k-', lw=0.6, alpha=0.5)
    ax2.axhline(0, color='gray', ls=':', lw=0.5)
    ax2.axvline(t_c, color='blue', ls='--', lw=0.8, alpha=0.4, label=r'$t_c$')
    ax2.axvline(t_r, color='purple', ls='--', lw=0.8, alpha=0.4, label=r'$t_r$')
    ax2.axvline(t_r/2, color='red', ls='--', lw=0.8, alpha=0.4,
                label=r'$t_r/2$')
    ax2.plot(t_now, np.interp(t_now, tlist, sigma_z), 'ro', ms=8, zorder=5)
    ax2.set_xlim(0, t_max)
    ax2.set_ylim(-1.1, 1.1)
    ax2.set_xlabel(r'$gt$', fontsize=11)
    ax2.set_ylabel(r'$\langle\sigma_z\rangle$', fontsize=11)
    ax2.set_title('Atomic Inversion', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=8, loc='lower right')

    fig.suptitle(
        r'Photon Number Distribution $P(n,t)$ — Jaynes–Cummings'
        + f' ($\\bar{{n}} = {n_bar}$)',
        fontsize=13, fontweight='bold', y=0.98
    )

    # Render to PIL
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
    buf.seek(0)
    frames.append(Image.open(buf).copy())
    buf.close()
    plt.close(fig)

    if (fi + 1) % 20 == 0:
        print(f"  Frame {fi+1}/{N_frames}")

# Save GIF
gif_path = '/home/claude/animations/anim_photon_number.gif'
# Add pauses at key times
durations = []
for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]
    if abs(t_now - t_r/2) < t_max / N_frames * 1.5:
        durations.append(300)   # pause at cat state
    elif abs(t_now - t_r) < t_max / N_frames * 1.5:
        durations.append(300)   # pause at revival
    elif t_now < 0.1:
        durations.append(200)   # pause at start
    else:
        durations.append(80)

frames[0].save(
    gif_path, save_all=True, append_images=frames[1:],
    duration=durations, loop=0
)
print(f"Animation saved: {gif_path}")
print("Done.")
