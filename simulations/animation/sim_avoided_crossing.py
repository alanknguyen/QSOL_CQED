#!/usr/bin/env python3
"""
Dressed-state avoided crossing in the Jaynes–Cummings model.

Animates the dressed-state energy eigenvalues as the detuning Δ sweeps
from -10g to +10g.  The bare-state energies cross at Δ = 0, but the
JC coupling opens an avoided crossing of width 2g√(n+1).

Overlays n = 0, 1, 5, 10 simultaneously so the √(n+1) scaling of the
splitting is visible.

This connects the analytical theory (Sec. II E–F on dressed states and
limiting cases) to a visual that makes the resonance/dispersive crossover
intuitive.

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
from PIL import Image
import io

# ── Parameters ──────────────────────────────────────────────────────
g = 1.0
Omega0 = 2 * g  # single-photon Rabi frequency
Delta_range = np.linspace(-10, 10, 500)
n_values = [0, 1, 5, 10]
colors_n = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
labels_n = [r'$n=0$', r'$n=1$', r'$n=5$', r'$n=10$']

# ── Eigenvalue computation ──────────────────────────────────────────
def dressed_energies(Delta, n, omega_c=1.0):
    """
    E±(n) = ℏωc(n + 1/2) ± (ℏ/2)Ωn(Δ)
    where Ωn(Δ) = sqrt(Δ² + Ω0²(n+1))

    We subtract the common energy ℏωc(n + 1/2) and plot in units of g.
    """
    Omega_n = np.sqrt(Delta**2 + Omega0**2 * (n + 1))
    E_plus  = +0.5 * Omega_n
    E_minus = -0.5 * Omega_n
    return E_plus, E_minus

def bare_energies(Delta):
    """Bare state energies: ±Δ/2 (in the rotating frame)."""
    return +0.5 * Delta, -0.5 * Delta

# ── Static figure: all n overlaid ───────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1.5, 1]})

ax = axes[0]
# Bare states (dashed)
E_bare_plus, E_bare_minus = bare_energies(Delta_range)
ax.plot(Delta_range, E_bare_plus, 'k--', lw=1.0, alpha=0.3, label='Bare states')
ax.plot(Delta_range, E_bare_minus, 'k--', lw=1.0, alpha=0.3)

for i, n in enumerate(n_values):
    Ep, Em = dressed_energies(Delta_range, n)
    ax.plot(Delta_range, Ep, color=colors_n[i], lw=2.0, label=labels_n[i])
    ax.plot(Delta_range, Em, color=colors_n[i], lw=2.0)

    # Annotate splitting at Δ=0
    split = Omega0 * np.sqrt(n + 1)
    ax.annotate('', xy=(10.5 - i*0.8, split/2), xytext=(10.5 - i*0.8, -split/2),
                arrowprops=dict(arrowstyle='<->', color=colors_n[i], lw=1.2))
    ax.text(10.8 - i*0.8, 0, f'{split:.1f}g',
            color=colors_n[i], fontsize=8, va='center', ha='left')

ax.set_xlabel(r'Detuning $\Delta/g$', fontsize=12)
ax.set_ylabel(r'Energy $/ \hbar g$', fontsize=12)
ax.set_title('Dressed-State Avoided Crossing', fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper left')
ax.set_xlim(-10, 13)
ax.set_ylim(-8, 8)
ax.axvline(0, color='gray', ls=':', lw=0.5)
ax.axhline(0, color='gray', ls=':', lw=0.5)

# Right panel: splitting vs n at Δ=0
ax2 = axes[1]
n_sweep = np.arange(0, 21)
splitting = Omega0 * np.sqrt(n_sweep + 1)
ax2.plot(n_sweep, splitting, 'ko-', ms=5, lw=1.5)
ax2.plot(n_sweep, Omega0 * np.sqrt(n_sweep + 1), 'r--', lw=1.5, alpha=0.5,
         label=r'$2g\sqrt{n+1}$')
for i, n in enumerate(n_values):
    ax2.plot(n, Omega0 * np.sqrt(n + 1), 'o', color=colors_n[i], ms=10,
             zorder=5)
ax2.set_xlabel(r'Photon number $n$', fontsize=12)
ax2.set_ylabel(r'Splitting $\Omega_n(0) / g$', fontsize=12)
ax2.set_title(r'$\sqrt{n+1}$ Scaling at $\Delta = 0$', fontsize=13,
              fontweight='bold')
ax2.legend(fontsize=10)

fig.suptitle(
    r'Jaynes–Cummings Dressed States: Avoided Crossing and $\sqrt{n+1}$ Scaling',
    fontsize=14, fontweight='bold', y=1.02
)
plt.tight_layout()
fig.savefig('/home/claude/figures/fig_avoided_crossing.png', dpi=200,
            bbox_inches='tight')
fig.savefig('/home/claude/figures/fig_avoided_crossing.pdf',
            bbox_inches='tight')
print("Static figure saved.")
plt.close()

# ── GIF animation: sweep Δ from -10g to +10g ───────────────────────
N_frames = 80
Delta_sweep = np.linspace(-10, 10, N_frames)
frames = []

print("Generating avoided crossing animation...")
for fi, Delta_now in enumerate(Delta_sweep):
    fig = plt.figure(figsize=(13, 5.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1], wspace=0.3)

    ax1 = fig.add_subplot(gs[0])

    # Full energy curves
    for i, n in enumerate(n_values):
        Ep, Em = dressed_energies(Delta_range, n)
        ax1.plot(Delta_range, Ep, color=colors_n[i], lw=1.5, alpha=0.3)
        ax1.plot(Delta_range, Em, color=colors_n[i], lw=1.5, alpha=0.3)

    # Bare states
    ax1.plot(Delta_range, E_bare_plus, 'k--', lw=0.8, alpha=0.2)
    ax1.plot(Delta_range, E_bare_minus, 'k--', lw=0.8, alpha=0.2)

    # Current detuning line
    ax1.axvline(Delta_now, color='gray', ls='-', lw=1.5, alpha=0.3)

    # Markers at current Δ
    for i, n in enumerate(n_values):
        Ep_now = 0.5 * np.sqrt(Delta_now**2 + Omega0**2 * (n + 1))
        Em_now = -Ep_now
        ax1.plot(Delta_now, Ep_now, 'o', color=colors_n[i], ms=8,
                 zorder=10, markeredgecolor='black', markeredgewidth=0.5)
        ax1.plot(Delta_now, Em_now, 'o', color=colors_n[i], ms=8,
                 zorder=10, markeredgecolor='black', markeredgewidth=0.5)

    ax1.set_xlabel(r'Detuning $\Delta/g$', fontsize=11)
    ax1.set_ylabel(r'Energy $/ \hbar g$', fontsize=11)
    ax1.set_title('Dressed-State Energies', fontsize=12, fontweight='bold')
    ax1.set_xlim(-10, 10)
    ax1.set_ylim(-8, 8)
    ax1.axhline(0, color='gray', ls=':', lw=0.3)

    # Legend
    for i, n in enumerate(n_values):
        ax1.plot([], [], 'o-', color=colors_n[i], label=labels_n[i], ms=5)
    ax1.legend(fontsize=9, loc='upper left')

    # Regime annotation
    if abs(Delta_now) < 1.5:
        regime = 'Resonant regime\nMax splitting'
        regime_color = 'lightyellow'
    elif abs(Delta_now) < 5:
        regime = 'Intermediate regime'
        regime_color = 'lightyellow'
    else:
        regime = 'Dispersive regime\nAC Stark shifts'
        regime_color = 'lightcyan'

    ax1.text(0.98, 0.05,
             f'$\\Delta/g = {Delta_now:.1f}$\n{regime}',
             transform=ax1.transAxes, fontsize=9, va='bottom', ha='right',
             bbox=dict(boxstyle='round', facecolor=regime_color, alpha=0.8))

    # Right panel: splitting diagram at current Δ
    ax2 = fig.add_subplot(gs[1])

    for i, n in enumerate(n_values):
        Omega_n = np.sqrt(Delta_now**2 + Omega0**2 * (n + 1))
        # Draw the two levels as horizontal bars
        y_center = i * 3
        ax2.barh(y_center + Omega_n/2, 0.8, height=0.15, left=-0.4,
                 color=colors_n[i], edgecolor='black', linewidth=0.5)
        ax2.barh(y_center - Omega_n/2, 0.8, height=0.15, left=-0.4,
                 color=colors_n[i], edgecolor='black', linewidth=0.5)
        # Arrow showing splitting
        ax2.annotate('', xy=(1.0, y_center + Omega_n/2),
                     xytext=(1.0, y_center - Omega_n/2),
                     arrowprops=dict(arrowstyle='<->', color=colors_n[i], lw=1.5))
        ax2.text(1.4, y_center, f'$\\Omega_{{{n}}} = {Omega_n:.2f}g$',
                 fontsize=9, va='center', color=colors_n[i])
        ax2.text(-1.0, y_center, labels_n[i], fontsize=10, va='center',
                 fontweight='bold', color=colors_n[i])

    ax2.set_xlim(-1.5, 4)
    ax2.set_ylim(-2, 11)
    ax2.set_title(f'Level Splittings at $\\Delta = {Delta_now:.1f}g$',
                  fontsize=11, fontweight='bold')
    ax2.set_axis_off()

    fig.suptitle(
        r'Jaynes–Cummings Avoided Crossing — Detuning Sweep',
        fontsize=13, fontweight='bold', y=0.99
    )

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=110, bbox_inches='tight')
    buf.seek(0)
    frames.append(Image.open(buf).copy())
    buf.close()
    plt.close(fig)

    if (fi + 1) % 20 == 0:
        print(f"  Frame {fi+1}/{N_frames}")

# Durations: pause near resonance
durations = []
for fi, D in enumerate(Delta_sweep):
    if abs(D) < 0.5:
        durations.append(200)
    elif abs(D) < 2:
        durations.append(120)
    else:
        durations.append(70)

gif_path = '/home/claude/animations/anim_avoided_crossing.gif'
frames[0].save(gif_path, save_all=True, append_images=frames[1:],
               duration=durations, loop=0)
print(f"Animation saved: {gif_path}")
print("Done.")
