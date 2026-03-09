#!/usr/bin/env python3
"""
Bloch sphere animation of the reduced atomic state during JC dynamics.

The reduced atomic density matrix is 2×2, mapping to a point (r_x, r_y, r_z)
inside the Bloch sphere with r_i = Tr[ρ_atom σ_i].

During Rabi oscillations the Bloch vector traces circles near the surface
(nearly pure state).  During collapse it spirals inward to the center
(maximally mixed, S = 1 bit).  At revival it pops back toward the surface.

The distance from the center is |r| = sqrt(1 - 2(1-P)) where P is the
purity.  Being INSIDE the sphere is a geometric signature of entanglement
with the field.

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
N_time = 800
tlist = np.linspace(0, t_max, N_time)

# ── JC Hamiltonian ──────────────────────────────────────────────────
a  = qt.tensor(qt.destroy(N_cav), qt.qeye(2))
sm = qt.tensor(qt.qeye(N_cav), qt.sigmam())
sp = qt.tensor(qt.qeye(N_cav), qt.sigmap())
sz = qt.tensor(qt.qeye(N_cav), qt.sigmaz())
sx = qt.tensor(qt.qeye(N_cav), qt.sigmax())
sy = qt.tensor(qt.qeye(N_cav), qt.sigmay())

H_JC = g * (a.dag() * sm + a * sp)
psi0 = qt.tensor(qt.coherent(N_cav, alpha), qt.basis(2, 0))

# ── Solve ───────────────────────────────────────────────────────────
print("Running JC dynamics...")
result = qt.mesolve(H_JC, psi0, tlist, e_ops=[])
rx = np.array([qt.expect(sx, result.states[i]) for i in range(N_time)])
ry = np.array([qt.expect(sy, result.states[i]) for i in range(N_time)])
rz = np.array([qt.expect(sz, result.states[i]) for i in range(N_time)])
bloch_r = np.sqrt(rx**2 + ry**2 + rz**2)

# Von Neumann entropy
print("Computing entanglement entropy...")
entropy = np.zeros(N_time)
for i in range(N_time):
    rho_atom = result.states[i].ptrace(1)
    entropy[i] = qt.entropy_vn(rho_atom, 2)

# ── Static figure: trajectory + entropy ─────────────────────────────
fig = plt.figure(figsize=(14, 6))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.35)

# 3D Bloch trajectory
ax1 = fig.add_subplot(gs[0], projection='3d')

# Draw Bloch sphere wireframe
u = np.linspace(0, 2*np.pi, 50)
v = np.linspace(0, np.pi, 30)
xs = np.outer(np.cos(u), np.sin(v))
ys = np.outer(np.sin(u), np.sin(v))
zs = np.outer(np.ones_like(u), np.cos(v))
ax1.plot_surface(xs, ys, zs, alpha=0.05, color='lightblue')

# Axes
ax1.plot([-1.2, 1.2], [0, 0], [0, 0], 'k-', lw=0.3, alpha=0.3)
ax1.plot([0, 0], [-1.2, 1.2], [0, 0], 'k-', lw=0.3, alpha=0.3)
ax1.plot([0, 0], [0, 0], [-1.2, 1.2], 'k-', lw=0.3, alpha=0.3)

# Color by time
colors = plt.cm.viridis(np.linspace(0, 1, N_time))
for i in range(N_time - 1):
    ax1.plot(rx[i:i+2], ry[i:i+2], rz[i:i+2],
             color=colors[i], lw=0.8, alpha=0.7)

# Mark key points
ax1.scatter(*[rx[0]], *[ry[0]], *[rz[0]], c='red', s=50, zorder=10,
            label=r'$t=0$ (surface)')
i_cat = np.argmin(np.abs(tlist - t_r/2))
ax1.scatter(*[rx[i_cat]], *[ry[i_cat]], *[rz[i_cat]], c='gold', s=50,
            zorder=10, marker='*', label=r'$t=t_r/2$ (center)')
i_rev = np.argmin(np.abs(tlist - t_r))
ax1.scatter(*[rx[i_rev]], *[ry[i_rev]], *[rz[i_rev]], c='lime', s=50,
            zorder=10, marker='^', label=r'$t=t_r$ (surface)')

ax1.set_xlabel(r'$\langle\sigma_x\rangle$', fontsize=9)
ax1.set_ylabel(r'$\langle\sigma_y\rangle$', fontsize=9)
ax1.set_zlabel(r'$\langle\sigma_z\rangle$', fontsize=9)
ax1.set_title('Bloch Sphere Trajectory', fontsize=11, fontweight='bold')
ax1.legend(fontsize=7, loc='upper left')
ax1.set_xlim(-1.1, 1.1)
ax1.set_ylim(-1.1, 1.1)
ax1.set_zlim(-1.1, 1.1)

# Bloch vector magnitude
ax2 = fig.add_subplot(gs[1])
ax2.plot(tlist, bloch_r, 'steelblue', lw=1.2)
ax2.axvline(t_c, color='blue', ls='--', lw=0.8, alpha=0.5, label=r'$t_c$')
ax2.axvline(t_r, color='purple', ls='--', lw=0.8, alpha=0.5, label=r'$t_r$')
ax2.axvline(t_r/2, color='red', ls='--', lw=0.8, alpha=0.5, label=r'$t_r/2$')
ax2.set_xlabel(r'$gt$', fontsize=11)
ax2.set_ylabel(r'$|\mathbf{r}|$', fontsize=11)
ax2.set_title('Bloch Vector Length', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=8)
ax2.text(0.5, 0.15, r'$|\mathbf{r}| = 0$: maximally mixed' + '\n'
         + r'(max entanglement)',
         transform=ax2.transAxes, fontsize=8, ha='center',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Entropy
ax3 = fig.add_subplot(gs[2])
ax3.plot(tlist, entropy, 'crimson', lw=1.2)
ax3.axvline(t_c, color='blue', ls='--', lw=0.8, alpha=0.5)
ax3.axvline(t_r, color='purple', ls='--', lw=0.8, alpha=0.5)
ax3.axvline(t_r/2, color='red', ls='--', lw=0.8, alpha=0.5)
ax3.axhline(1, color='gray', ls=':', lw=0.5)
ax3.set_xlabel(r'$gt$', fontsize=11)
ax3.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=11)
ax3.set_title('Entanglement Entropy', fontsize=11, fontweight='bold')
ax3.set_ylim(-0.05, 1.1)

fig.suptitle(
    r'Bloch Sphere Dynamics of the Reduced Atomic State — JC Model'
    + f' ($\\bar{{n}} = {n_bar}$)',
    fontsize=13, fontweight='bold', y=1.02
)
plt.tight_layout()
fig.savefig('/home/claude/figures/fig_bloch_sphere_trajectory.png', dpi=200,
            bbox_inches='tight')
fig.savefig('/home/claude/figures/fig_bloch_sphere_trajectory.pdf',
            bbox_inches='tight')
print("Static figure saved.")
plt.close()

# ── GIF animation ───────────────────────────────────────────────────
N_frames = 100
frame_times = np.linspace(0, t_max, N_frames)
frame_indices = [np.argmin(np.abs(tlist - ft)) for ft in frame_times]

print("Generating Bloch sphere animation...")
frames = []

for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]

    fig = plt.figure(figsize=(13, 5.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.3, 1, 1], wspace=0.3)

    # ─ 3D Bloch sphere
    ax1 = fig.add_subplot(gs[0], projection='3d')

    # Wireframe sphere
    ax1.plot_surface(xs, ys, zs, alpha=0.04, color='lightblue')
    # Great circles
    theta_c = np.linspace(0, 2*np.pi, 100)
    ax1.plot(np.cos(theta_c), np.sin(theta_c), np.zeros_like(theta_c),
             'k-', lw=0.3, alpha=0.15)
    ax1.plot(np.cos(theta_c), np.zeros_like(theta_c), np.sin(theta_c),
             'k-', lw=0.3, alpha=0.15)
    ax1.plot(np.zeros_like(theta_c), np.cos(theta_c), np.sin(theta_c),
             'k-', lw=0.3, alpha=0.15)

    # Trajectory up to now (trail)
    trail_len = min(idx, 60)  # show last ~60 steps
    start = max(0, idx - trail_len)
    trail_alpha = np.linspace(0.1, 0.9, idx - start + 1) if idx > start else [0.9]
    for j in range(start, idx):
        a_val = trail_alpha[j - start]
        ax1.plot(rx[j:j+2], ry[j:j+2], rz[j:j+2],
                 color=(0.2, 0.4, 0.8, a_val), lw=1.2)

    # Current point
    ax1.scatter([rx[idx]], [ry[idx]], [rz[idx]], c='red', s=80, zorder=10,
                edgecolors='darkred', linewidth=0.8)
    # Line from center to point (shows |r|)
    ax1.plot([0, rx[idx]], [0, ry[idx]], [0, rz[idx]],
             'r-', lw=1.0, alpha=0.5)

    # Axis labels
    for coord, label, pos in [
        ([1.3, 0, 0], r'$x$', 'x'),
        ([0, 1.3, 0], r'$y$', 'y'),
        ([0, 0, 1.3], r'$|e\rangle$', 'z')
    ]:
        ax1.text(*coord, label, fontsize=9, ha='center')
    ax1.text(0, 0, -1.3, r'$|g\rangle$', fontsize=9, ha='center')

    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-1.2, 1.2)
    ax1.set_zlim(-1.2, 1.2)
    ax1.set_axis_off()
    ax1.view_init(elev=20, azim=30 + t_now * 2)  # slow rotation
    ax1.set_title(f'Bloch Sphere\n$gt = {t_now:.1f}$',
                  fontsize=11, fontweight='bold')

    # ─ Bloch magnitude
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(tlist[:idx+1], bloch_r[:idx+1], 'steelblue', lw=1.2)
    ax2.plot(tlist[idx:], bloch_r[idx:], color='steelblue', lw=0.5, alpha=0.15)
    ax2.plot(t_now, bloch_r[idx], 'ro', ms=7, zorder=5)
    ax2.axvline(t_r, color='purple', ls='--', lw=0.7, alpha=0.4)
    ax2.axvline(t_r/2, color='red', ls='--', lw=0.7, alpha=0.4)
    ax2.set_xlim(0, t_max)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel(r'$gt$', fontsize=10)
    ax2.set_ylabel(r'$|\mathbf{r}|$ (purity)', fontsize=10)
    ax2.set_title('Bloch Vector Length', fontsize=11, fontweight='bold')

    # Annotation for regime
    r_now = bloch_r[idx]
    if r_now > 0.7:
        regime = 'Nearly pure\n(weak entanglement)'
    elif r_now > 0.3:
        regime = 'Partially mixed\n(growing entanglement)'
    else:
        regime = 'Maximally mixed\n(max entanglement)'
    ax2.text(0.98, 0.95, regime, transform=ax2.transAxes, fontsize=8,
             va='top', ha='right',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # ─ Entropy
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(tlist[:idx+1], entropy[:idx+1], 'crimson', lw=1.2)
    ax3.plot(tlist[idx:], entropy[idx:], color='crimson', lw=0.5, alpha=0.15)
    ax3.plot(t_now, entropy[idx], 'ro', ms=7, zorder=5)
    ax3.axhline(1, color='gray', ls=':', lw=0.5)
    ax3.axvline(t_r, color='purple', ls='--', lw=0.7, alpha=0.4)
    ax3.axvline(t_r/2, color='red', ls='--', lw=0.7, alpha=0.4)
    ax3.set_xlim(0, t_max)
    ax3.set_ylim(-0.05, 1.1)
    ax3.set_xlabel(r'$gt$', fontsize=10)
    ax3.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=10)
    ax3.set_title('Entanglement Entropy', fontsize=11, fontweight='bold')

    fig.suptitle(
        r'Bloch Sphere: Atomic State Inside $\Leftrightarrow$ '
        + r'Entangled with Field',
        fontsize=12, fontweight='bold', y=0.99
    )

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=110, bbox_inches='tight')
    buf.seek(0)
    frames.append(Image.open(buf).copy())
    buf.close()
    plt.close(fig)

    if (fi + 1) % 20 == 0:
        print(f"  Frame {fi+1}/{N_frames}")

# Save GIF with pauses
durations = []
for fi, idx in enumerate(frame_indices):
    t_now = tlist[idx]
    if abs(t_now - t_r/2) < t_max / N_frames * 1.5:
        durations.append(300)
    elif abs(t_now - t_r) < t_max / N_frames * 1.5:
        durations.append(300)
    elif t_now < 0.1:
        durations.append(200)
    else:
        durations.append(80)

gif_path = '/home/claude/animations/anim_bloch_sphere.gif'
frames[0].save(gif_path, save_all=True, append_images=frames[1:],
               duration=durations, loop=0)
print(f"Animation saved: {gif_path}")
print("Done.")
