#!/usr/bin/env python3
"""
Animation: Photon blockade transition.

Sweeps g/κ from 0.1 to 12, showing:
  - Left: JC dressed-state energy ladder (units of g, constant layout)
  - Right: Live cavity transmission spectrum ⟨n⟩(Δ) that morphs from
    a single Lorentzian into a vacuum Rabi doublet as g crosses κ

With g^(2)(0) and ⟨n⟩ shown as numerical readouts on the figure.

Produces: anim_g2_blockade.gif
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    destroy, sigmaz, sigmap, sigmam, tensor, qeye,
    steadystate, expect
)
import imageio
import os, shutil

fig_dir = '/home/claude/animations'
os.makedirs(fig_dir, exist_ok=True)
frame_dir = '/home/claude/_frames_g2'
if os.path.exists(frame_dir):
    shutil.rmtree(frame_dir)
os.makedirs(frame_dir, exist_ok=True)

N_cav = 15  # small is fine for weak drive

def build_driven_jc(g_val, kappa, gamma, epsilon, delta=0.0):
    a  = tensor(destroy(N_cav), qeye(2))
    sm = tensor(qeye(N_cav), sigmam())
    sp = tensor(qeye(N_cav), sigmap())
    sz = tensor(qeye(N_cav), sigmaz())
    H = (delta * a.dag() * a
         + delta / 2.0 * sz
         + g_val * (a.dag() * sm + a * sp)
         + epsilon * (a.dag() + a))
    c_ops = [np.sqrt(kappa) * a]
    if gamma > 0:
        c_ops.append(np.sqrt(gamma) * sm)
    return H, c_ops, a

kappa = 1.0; gamma = 0.1; epsilon = 0.3
g_over_kappa = np.linspace(0.1, 12, 50)

# Detuning grid for transmission spectrum
delta_grid = np.linspace(-18, 18, 100)  # in units of κ (=1)

# ── Precompute everything ────────────────────────────────────
print("Precomputing transmission spectra + g^(2)(0)...")
spectra = []
g2_vals = []
n_on_res = []

for i, gk in enumerate(g_over_kappa):
    g_val = gk * kappa

    # On-resonance g^(2)(0)
    H0, c0, a0 = build_driven_jc(g_val, kappa, gamma, epsilon, delta=0.0)
    rho0 = steadystate(H0, c0)
    n0 = np.real(expect(a0.dag() * a0, rho0))
    g2_0 = np.real(expect(a0.dag() * a0.dag() * a0 * a0, rho0)) / n0**2 if n0 > 1e-15 else 2.0
    g2_vals.append(g2_0)
    n_on_res.append(n0)

    # Transmission spectrum: ⟨n⟩ vs Δ
    spec = np.zeros(len(delta_grid))
    for j, delta in enumerate(delta_grid):
        H, c_ops, a = build_driven_jc(g_val, kappa, gamma, epsilon, delta=delta)
        rho_ss = steadystate(H, c_ops)
        spec[j] = np.real(expect(a.dag() * a, rho_ss))
    spectra.append(spec)

    if (i + 1) % 5 == 0:
        print(f"  {i+1}/{len(g_over_kappa)} (g/κ = {gk:.1f})")

g2_vals = np.array(g2_vals)
n_on_res = np.array(n_on_res)

# Global y-max for spectrum
spec_ymax = max(np.max(s) for s in spectra) * 1.15

# ── Render frames ────────────────────────────────────────────
print("Rendering frames...")

manifold_offsets = [0, 5, 11]
manifold_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
manifold_labels = ['$n=0$', '$n=1$', '$n=2$']

for i in range(len(g_over_kappa)):
    gk = g_over_kappa[i]
    g_val = gk * kappa
    kappa_over_g = kappa / g_val

    fig, (ax_e, ax_s) = plt.subplots(1, 2, figsize=(15, 6.5),
                                      gridspec_kw={'width_ratios': [1, 1.2]})

    # ════════════════════════════════════════════════════════════
    # LEFT: JC energy ladder
    # ════════════════════════════════════════════════════════════
    for n_man in range(3):
        y0 = manifold_offsets[n_man]
        split = np.sqrt(n_man + 1)
        color = manifold_colors[n_man]

        # Bare level
        ax_e.plot([0.5, 1.5], [y0, y0], 'k--', lw=1, alpha=0.3)

        # Dressed states with linewidth bands
        y_plus = y0 + split
        y_minus = y0 - split
        band_hw = min(kappa_over_g * 0.5, split * 0.4)

        ax_e.fill_between([2.5, 4.5], y_plus - band_hw, y_plus + band_hw,
                         color=color, alpha=0.25)
        ax_e.fill_between([2.5, 4.5], y_minus - band_hw, y_minus + band_hw,
                         color=color, alpha=0.25)
        ax_e.plot([2.5, 4.5], [y_plus, y_plus], color=color, lw=2.5)
        ax_e.plot([2.5, 4.5], [y_minus, y_minus], color=color, lw=2.5)

        # Splitting annotation
        ax_e.annotate('', xy=(5.0, y_plus), xytext=(5.0, y_minus),
                     arrowprops=dict(arrowstyle='<->', color=color, lw=1.2))
        ax_e.text(5.2, y0, rf'$2\sqrt{{{n_man+1}}}$',
                 fontsize=9, color=color, va='center')

        ax_e.text(-0.3, y0, manifold_labels[n_man], fontsize=11,
                 ha='right', va='center', fontweight='bold')

    # Drive arrow
    y_target = manifold_offsets[0] - 1.0
    ax_e.annotate('', xy=(3.5, y_target + 0.15), xytext=(3.5, -3.5),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2.5))
    ax_e.text(3.8, -2.5, r'$\omega_L$', fontsize=12, color='red')

    # Blocked 2nd photon
    if gk > 2.0:
        y_from = manifold_offsets[0] + 1.0
        y_to = manifold_offsets[1] - np.sqrt(2)
        mid_y = (y_from + y_to) / 2
        ax_e.plot([3.5, 3.5], [y_from + 0.2, y_to - 0.2],
                 'r--', lw=1.5, alpha=0.6)
        ax_e.text(3.7, mid_y, '✗ blocked', fontsize=10, color='red',
                 va='center', fontweight='bold')

    # κ/g readout
    ax_e.text(0.5, -3.5, rf'$\kappa/g = {kappa_over_g:.2f}$',
             fontsize=11, color='purple',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lavender', alpha=0.8))

    ax_e.set_xlim(-0.8, 6.0); ax_e.set_ylim(-4.5, 15.5)
    ax_e.set_ylabel('Energy  ($E / \\hbar g$)', fontsize=12)
    ax_e.set_title(rf'JC Dressed-State Ladder  ($g/\kappa = {gk:.1f}$)', fontsize=13)
    ax_e.set_xticks([1.0, 3.5])
    ax_e.set_xticklabels(['Bare', 'Dressed'], fontsize=11)
    ax_e.tick_params(axis='y', labelsize=9)

    # ════════════════════════════════════════════════════════════
    # RIGHT: Live transmission spectrum ⟨n⟩(Δ)
    # ════════════════════════════════════════════════════════════
    spec = spectra[i]

    # Fill under the curve
    ax_s.fill_between(delta_grid, 0, spec, color='#ff7f0e', alpha=0.25)
    ax_s.plot(delta_grid, spec, color='#d35400', linewidth=2.5)

    # Mark the ±g positions (expected peak locations)
    if gk > 1.5:
        ax_s.axvline(-g_val, color='red', ls='--', alpha=0.35, lw=1.5)
        ax_s.axvline(g_val, color='red', ls='--', alpha=0.35, lw=1.5)
        # Label splitting
        peak_y = spec_ymax * 0.92
        ax_s.annotate('', xy=(g_val, peak_y), xytext=(-g_val, peak_y),
                     arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
        ax_s.text(0, peak_y + spec_ymax * 0.03, rf'$2g = {2*g_val:.1f}\,\kappa$',
                 fontsize=11, ha='center', color='red', fontweight='bold')

    ax_s.set_xlim(delta_grid[0], delta_grid[-1])
    ax_s.set_ylim(0, spec_ymax)
    ax_s.set_xlabel(r'Detuning  $\Delta / \kappa$', fontsize=13)
    ax_s.set_ylabel(r'$\langle n \rangle$', fontsize=13)
    ax_s.set_title('Cavity Transmission Spectrum', fontsize=13)

    # Numerical readouts in a box
    g2_now = g2_vals[i]
    if g2_now < 0.5:
        g2_color = '#0044cc'
    elif g2_now < 1.2:
        g2_color = '#444444'
    else:
        g2_color = '#cc2200'

    info_text = (rf'$g^{{(2)}}(0) = {g2_now:.3f}$'
                 '\n'
                 rf'$\langle n \rangle_{{res}} = {n_on_res[i]:.4f}$')
    ax_s.text(0.97, 0.75, info_text,
             transform=ax_s.transAxes, fontsize=13, fontweight='bold',
             color=g2_color, ha='right', va='top',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      edgecolor=g2_color, alpha=0.9, linewidth=2))

    # Regime label
    if gk < 1.0:
        regime = 'Weak coupling\n(single peak)'
    elif gk < 3.0:
        regime = 'Onset of\nstrong coupling'
    else:
        regime = 'Strong coupling\n(vacuum Rabi doublet)'
    ax_s.text(0.03, 0.95, regime,
             transform=ax_s.transAxes, fontsize=10, va='top',
             style='italic', color='#555555',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

    fig.suptitle(r'Photon Blockade in the Driven Jaynes–Cummings Cavity'
                 '\n' rf'$\varepsilon/\kappa = {epsilon/kappa:.1f}$, '
                 rf'$\gamma/\kappa = {gamma/kappa:.1f}$',
                 fontsize=14, y=0.98)

    fig.subplots_adjust(top=0.85)
    fig.savefig(f'{frame_dir}/frame_{i:04d}.png', dpi=110,
                facecolor='white', pad_inches=0.3)
    plt.close(fig)

    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{len(g_over_kappa)} frames")

# ── Assemble GIF ─────────────────────────────────────────────
print("Assembling GIF...")
frames = []
for i in range(len(g_over_kappa)):
    frames.append(imageio.v2.imread(f'{frame_dir}/frame_{i:04d}.png'))

durations = [0.22] * len(g_over_kappa)
durations[0] = 1.5   # pause on weak coupling
durations[-1] = 2.5   # pause on strong blockade
# Pause at the splitting onset
split_idx = np.argmin(np.abs(g_over_kappa - 1.0))
durations[split_idx] = 1.2
# Pause when doublet is cleanly resolved
clean_idx = np.argmin(np.abs(g_over_kappa - 3.0))
durations[clean_idx] = 1.0

outpath = f'{fig_dir}/anim_g2_blockade.gif'
imageio.mimsave(outpath, frames, duration=durations, loop=0)
print(f"Saved: {outpath}")
shutil.rmtree(frame_dir)
