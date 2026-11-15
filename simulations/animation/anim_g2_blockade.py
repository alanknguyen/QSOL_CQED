#!/usr/bin/env python3
"""
Animation: photon blockade in the driven Jaynes-Cummings cavity as g/kappa grows.

  - Left: the JC ladder (ground state and the first three doublets, splittings
    2 g sqrt(n+1), level widths set by the polariton linewidth (kappa+gamma)/2).
    The laser is tuned to the lower polariton, omega_L = omega_c - g.  A second
    photon of the same frequency would land at 2 omega_c - 2g, which misses the
    next level 2 omega_c - sqrt(2) g by (2 - sqrt(2)) g: once that exceeds the
    linewidth the second photon is blocked.
  - Right: cavity transmission <n>(Delta), which changes from a single
    Lorentzian into the vacuum Rabi doublet at Delta = +-g, with the values of
    g2(0) and <n> for the drive on the lower polariton (Delta = omega_c - omega_L = g).

Driving on BARE-cavity resonance (Delta = 0) in strong coupling is detuned
from both polaritons, <n> collapses and g2(0) >> 1 (bunching): that is not
the blockade signature.

Generates: animations/anim_g2_blockade.gif
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from qutip import destroy, sigmaz, sigmam, tensor, qeye, steadystate, expect

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import use_style, fig_to_frame, save_gif, COLORS, ANIM_DIR

use_style()

N_cav = 15
kappa, gamma, epsilon = 1.0, 0.1, 0.3
g_over_kappa = np.linspace(0.1, 12, 50)
delta_grid = np.linspace(-18, 18, 361)           # laser detuning in units of kappa

a = tensor(destroy(N_cav), qeye(2))
sm = tensor(qeye(N_cav), sigmam())
sz = tensor(qeye(N_cav), sigmaz())
n_op = a.dag() * a
c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]


def steady(g_val, delta):
    H = delta * n_op + 0.5 * delta * sz + g_val * (a.dag() * sm + a * sm.dag()) + epsilon * (a + a.dag())
    rho = steadystate(H, c_ops)
    n = expect(n_op, rho).real
    return n, expect(a.dag() * a.dag() * a * a, rho).real / n ** 2


print("Computing transmission spectra and g2(0) ...")
spectra, g2_vals, n_pol = [], [], []
for k, gk in enumerate(g_over_kappa):
    g_val = gk * kappa
    n0, g2 = steady(g_val, g_val)                # drive on the lower polariton
    n_pol.append(n0); g2_vals.append(g2)
    spectra.append(np.array([steady(g_val, d)[0] for d in delta_grid]))
    if (k + 1) % 10 == 0:
        print(f"  {k + 1}/{len(g_over_kappa)}  (g/kappa = {gk:.1f}, g2 = {g2:.3f})")
spec_ymax = 1.15 * max(s.max() for s in spectra)

# Ladder geometry: energies in units of g, with one cavity quantum drawn as STEP units
STEP = 5.0
manifold_colors = [COLORS['blue'], COLORS['orange'], COLORS['green']]
y_ground = -STEP
frames = []
print("Rendering frames ...")
for i, gk in enumerate(g_over_kappa):
    g_val = gk * kappa
    half_width = min((kappa + gamma) / (4 * g_val), 2.4)       # polariton HWHM in units of g
    fig, (ax_e, ax_s) = plt.subplots(1, 2, figsize=(8.6, 3.7), gridspec_kw={'width_ratios': [1, 1.25]})

    ax_e.plot([1.0, 4.0], [y_ground, y_ground], color='k', lw=1.5)
    ax_e.text(0.8, y_ground, r'$|g,0\rangle$', ha='right', va='center', fontsize=8)
    for n, color in enumerate(manifold_colors):
        y0 = n * STEP
        split = np.sqrt(n + 1)
        ax_e.plot([1.0, 4.0], [y0, y0], color='0.6', lw=0.5, ls='--')
        for sgn in (+1, -1):
            y = y0 + sgn * split
            ax_e.fill_between([1.0, 4.0], y - half_width, y + half_width, color=color, alpha=0.22, lw=0)
            ax_e.plot([1.0, 4.0], [y, y], color=color, lw=1.6)
        ax_e.annotate('', xy=(4.5, y0 + split), xytext=(4.5, y0 - split),
                      arrowprops=dict(arrowstyle='<->', color=color, lw=0.8, shrinkA=0, shrinkB=0))
        ax_e.text(4.7, y0, rf'$2g\sqrt{{{n + 1}}}$', fontsize=8, color=color, va='center')
        ax_e.text(0.8, y0, rf'$n = {n}$', ha='right', va='center', fontsize=8)

    # first laser photon: |g,0> -> lower polariton; second photon of the same energy
    y1 = -1.0
    y2 = y1 + (STEP - 1.0)                       # = 2 omega_c - 2 g in these units
    ax_e.annotate('', xy=(2.0, y1), xytext=(2.0, y_ground),
                  arrowprops=dict(arrowstyle='->', color=COLORS['red'], lw=1.4, shrinkA=0, shrinkB=0))
    ax_e.text(2.15, 0.5 * (y_ground + y1), r'$\omega_L$', fontsize=8.5, color=COLORS['red'], va='center')
    blocked = (2 - np.sqrt(2)) * g_val > 0.5 * (kappa + gamma)
    ax_e.annotate('', xy=(3.0, y2), xytext=(3.0, y1),
                  arrowprops=dict(arrowstyle='->', color=COLORS['red'], lw=1.0, ls='--',
                                  shrinkA=0, shrinkB=0))
    ax_e.plot([2.6, 3.4], [y2, y2], color=COLORS['red'], lw=0.8, ls=':')
    ax_e.text(3.5, y2, 'blocked' if blocked else 'not resolved', fontsize=7.5,
              color=COLORS['red'], va='center')
    ax_e.set_xlim(-0.6, 6.0); ax_e.set_ylim(y_ground - 1.2, 2 * STEP + 3.0)
    ax_e.set_xticks([]); ax_e.set_yticks([])
    for sp in ax_e.spines.values():
        sp.set_visible(False)
    miss = (2 - np.sqrt(2)) * gk                      # two digits below 1: 0.1 would misstate 0.06
    ax_e.set_title(rf'JC ladder, $g/\kappa = {gk:.1f}$' + '\n'
                   + rf'second photon misses by $(2-\sqrt{{2}})\,g = {miss:.{2 if miss < 1 else 1}f}\,\kappa$',
                   fontsize=9)

    ax_s.fill_between(delta_grid, 0, spectra[i], color=COLORS['orange'], alpha=0.25, lw=0)
    ax_s.plot(delta_grid, spectra[i], color=COLORS['orange'], lw=1.3)
    if gk > 1.5:
        for sgn in (-1, 1):
            ax_s.axvline(sgn * g_val, color=COLORS['red'], ls=':', lw=0.7)
    ax_s.set_xlim(delta_grid[0], delta_grid[-1]); ax_s.set_ylim(0, spec_ymax)
    ax_s.set_xlabel(r'laser detuning $\Delta/\kappa$'); ax_s.set_ylabel(r'$\langle n \rangle$')
    # below g ~ kappa the single peak carries a narrow notch of width gamma (absorption by the atom)
    regime = ('weak coupling: one notched peak' if gk < 1.0 else
              'onset of strong coupling' if gk < 3.0 else r'vacuum Rabi doublet at $\Delta = \pm g$')
    ax_s.set_title(f'cavity transmission: {regime}', fontsize=9)
    ax_s.text(0.97, 0.95, rf'drive at $\Delta = g$:' + '\n' + rf'$g^{{(2)}}(0) = {g2_vals[i]:.3f}$'
              + '\n' + rf'$\langle n \rangle = {n_pol[i]:.3f}$',
              transform=ax_s.transAxes, fontsize=8.5, ha='right', va='top',
              bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='0.5', lw=0.6))
    ax_s.text(0.03, 0.95, rf'$\varepsilon/\kappa = {epsilon / kappa:g}$' + '\n'
              + rf'$\gamma/\kappa = {gamma / kappa:g}$', transform=ax_s.transAxes, fontsize=7.5, va='top')
    frames.append(fig_to_frame(fig, dpi=110))

durations = [220] * len(g_over_kappa)
durations[0] = 1500                                           # weak coupling
durations[int(np.argmin(np.abs(g_over_kappa - 1.0)))] = 1200  # g = kappa
durations[int(np.argmin(np.abs(g_over_kappa - 3.0)))] = 1000  # doublet resolved
durations[-1] = 2500
save_gif(frames, ANIM_DIR / 'anim_g2_blockade.gif', durations)
