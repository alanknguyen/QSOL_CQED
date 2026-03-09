#!/usr/bin/env python3
"""
Q-function vs Wigner function side-by-side during JC dynamics.

The Husimi Q-function Q(α) = <α|ρ|α>/π is always non-negative (it's
the Wigner function convolved with a Gaussian of width 1/2).

At t=0 both look identical (Gaussian).  At the cat-state time t = t_r/2,
Q shows two blobs but NO interference fringes, while W shows the fringes
with negative regions.  This visually demonstrates why Wigner negativity
is a stronger non-classicality test, which is discussed in Section V
but never shown in the paper.

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import qutip as qt
from PIL import Image
import io

# ── Parameters ──────────────────────────────────────────────────────
g = 1.0
n_bar = 10
alpha = np.sqrt(n_bar)
N_cav = 40
t_r = 2 * np.pi * np.sqrt(n_bar) / g
t_c = np.pi / g
t_max = 1.3 * t_r
N_time = 600
tlist = np.linspace(0, t_max, N_time)

# Phase space grid
N_grid = 150  # for Wigner
xvec = np.linspace(-7, 7, N_grid)
N_grid_q = 150  # for Q
xvec_q = np.linspace(-7, 7, N_grid_q)

# ── JC Hamiltonian ──────────────────────────────────────────────────
a  = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
sp = qt.tensor(qt.qeye(N_cav), qt.sigmap())
sz = qt.tensor(qt.qeye(N_cav), qt.sigmaz())

H_JC = g * (a.dag() * sm + a * sp)
psi0 = qt.tensor(qt.coherent(N_cav, alpha), qt.basis(2, 0))

# ── Solve ───────────────────────────────────────────────────────────
print("Running JC dynamics...")
result = qt.mesolve(H_JC, psi0, tlist, e_ops=[])
sigma_z = np.array([qt.expect(sz, result.states[i]) for i in range(N_time)])

# ── Wigner negativity volume ────────────────────────────────────────
def wigner_negativity(W, dx):
    """δ = ∫|W| dx dp - 1"""
    return np.sum(np.abs(W)) * dx * dx - 1.0

# ── Static figure: 6 snapshots side-by-side ─────────────────────────
snapshot_times = [0, t_c * 0.5, t_c * 2, t_r / 2, t_r * 0.75, t_r]
snapshot_labels = [
    r'$t = 0$', r'$t = 0.5\,t_c$', r'$t = 2\,t_c$',
    r'$t = t_r/2$ (cat)', r'$t = 0.75\,t_r$', r'$t = t_r$ (revival)'
]

fig, axes = plt.subplots(2, 6, figsize=(20, 7))
fig.suptitle(
    r'Wigner $W(x,p)$ vs Husimi $Q(x,p)$ — Jaynes–Cummings'
    + f' ($\\bar{{n}} = {n_bar}$)',
    fontsize=14, fontweight='bold'
)

for col, (t_snap, label) in enumerate(zip(snapshot_times, snapshot_labels)):
    idx = np.argmin(np.abs(tlist - t_snap))
    rho_field = result.states[idx].ptrace(0)

    # Wigner
    W = qt.wigner(rho_field, xvec, xvec)
    dx = xvec[1] - xvec[0]
    delta = wigner_negativity(W, dx)

    # Q function
    Q = qt.qfunc(rho_field, xvec_q, xvec_q)

    # Top row: Wigner
    ax_w = axes[0, col]
    wlim = max(abs(W.min()), abs(W.max()))
    ax_w.contourf(xvec, xvec, W, levels=60, cmap='RdBu_r',
                  vmin=-wlim, vmax=wlim)
    ax_w.set_aspect('equal')
    ax_w.set_title(label, fontsize=9)
    if col == 0:
        ax_w.set_ylabel(r'Wigner $W(x,p)$' + '\n' + r'$p$', fontsize=10)
    ax_w.text(0.02, 0.98, f'$\\delta = {delta:.3f}$', transform=ax_w.transAxes,
              fontsize=7, va='top', color='black',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    ax_w.tick_params(labelsize=6)

    # Bottom row: Q
    ax_q = axes[1, col]
    ax_q.contourf(xvec_q, xvec_q, Q, levels=60, cmap='inferno')
    ax_q.set_aspect('equal')
    if col == 0:
        ax_q.set_ylabel(r'Husimi $Q(\alpha)$' + '\n' + r'$p$', fontsize=10)
    ax_q.set_xlabel(r'$x$', fontsize=9)
    ax_q.text(0.02, 0.98, r'$Q \geq 0$ always', transform=ax_q.transAxes,
              fontsize=7, va='top', color='white',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5))
    ax_q.tick_params(labelsize=6)

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig('/home/claude/figures/fig_q_vs_wigner.png', dpi=200,
            bbox_inches='tight')
fig.savefig('/home/claude/figures/fig_q_vs_wigner.pdf', bbox_inches='tight')
print("Static figure saved.")
plt.close()

# ── GIF animation ───────────────────────────────────────────────────
N_frames = 80
frame_times = np.linspace(0, t_max, N_frames)
frame_indices = [np.argmin(np.abs(tlist - ft)) for ft in frame_times]

# Use smaller grid for animation speed
N_anim = 100
xvec_a = np.linspace(-7, 7, N_anim)
dx_a = xvec_a[1] - xvec_a[0]

print("Generating Q vs Wigner animation...")
frames = []

for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]
    rho_field = result.states[idx].ptrace(0)

    W = qt.wigner(rho_field, xvec_a, xvec_a)
    Q = qt.qfunc(rho_field, xvec_a, xvec_a)
    delta = wigner_negativity(W, dx_a)

    fig = plt.figure(figsize=(14, 5.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.8], wspace=0.25)

    # ─ Wigner function
    ax1 = fig.add_subplot(gs[0])
    wlim = max(abs(W.min()), abs(W.max()), 0.01)
    ax1.contourf(xvec_a, xvec_a, W, levels=50, cmap='RdBu_r',
                 vmin=-wlim, vmax=wlim)
    ax1.set_aspect('equal')
    ax1.set_xlabel(r'$x$', fontsize=11)
    ax1.set_ylabel(r'$p$', fontsize=11)
    ax1.set_title(f'Wigner $W(x,p)$\n$\\delta = {delta:.3f}$',
                  fontsize=11, fontweight='bold')

    # Highlight negative regions
    if delta > 0.01:
        ax1.contour(xvec_a, xvec_a, W, levels=[0], colors='black',
                    linewidths=0.5, linestyles='--')

    # ─ Q function
    ax2 = fig.add_subplot(gs[1])
    ax2.contourf(xvec_a, xvec_a, Q, levels=50, cmap='inferno')
    ax2.set_aspect('equal')
    ax2.set_xlabel(r'$x$', fontsize=11)
    ax2.set_ylabel(r'$p$', fontsize=11)
    ax2.set_title(r'Husimi $Q(\alpha) \geq 0$' + '\n(no fringes)',
                  fontsize=11, fontweight='bold')

    # ─ Right panel: cross-sections + inversion
    ax3 = fig.add_subplot(gs[2])
    mid = N_anim // 2

    # W and Q cross-sections at p = 0
    ax3.plot(xvec_a, W[mid, :], 'b-', lw=1.5, label=r'$W(x, 0)$')
    ax3.plot(xvec_a, Q[mid, :], 'r-', lw=1.5, alpha=0.7,
             label=r'$Q(x, 0)$')
    ax3.axhline(0, color='gray', ls=':', lw=0.5)
    ax3.fill_between(xvec_a, W[mid, :], 0,
                     where=(W[mid, :] < 0), color='blue', alpha=0.2)
    ax3.set_xlabel(r'$x$', fontsize=10)
    ax3.set_ylabel('Value', fontsize=10)
    ax3.set_title('Cross-section at $p=0$', fontsize=10, fontweight='bold')
    ax3.legend(fontsize=8)
    ax3.set_ylim(-0.3, 0.4)
    ax3.set_xlim(-7, 7)

    # Phase annotation
    if t_now < t_c:
        phase = 'Rabi oscillations'
    elif t_now < t_r * 0.35:
        phase = 'Collapse'
    elif abs(t_now - t_r/2) < 2.0:
        phase = 'Cat state — fringes in W only!'
    elif t_now < t_r * 0.85:
        phase = 'Pre-revival'
    elif abs(t_now - t_r) < 2.0:
        phase = 'Revival'
    else:
        phase = 'Post-revival'

    fig.suptitle(
        f'$gt = {t_now:.1f}$ — {phase}\n'
        + r'Wigner negativity is a stronger non-classicality witness than $Q \geq 0$',
        fontsize=12, fontweight='bold', y=1.01
    )

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=110, bbox_inches='tight')
    buf.seek(0)
    frames.append(Image.open(buf).copy())
    buf.close()
    plt.close(fig)

    if (fi + 1) % 20 == 0:
        print(f"  Frame {fi+1}/{N_frames}")

# Durations
durations = []
for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]
    if abs(t_now - t_r/2) < t_max / N_frames * 1.5:
        durations.append(350)
    elif abs(t_now - t_r) < t_max / N_frames * 1.5:
        durations.append(250)
    elif t_now < 0.1:
        durations.append(200)
    else:
        durations.append(80)

gif_path = '/home/claude/animations/anim_q_vs_wigner.gif'
frames[0].save(gif_path, save_all=True, append_images=frames[1:],
               duration=durations, loop=0)
print(f"Animation saved: {gif_path}")
print("Done.")
