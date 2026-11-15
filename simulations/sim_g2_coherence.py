#!/usr/bin/env python3
"""
Second-order coherence g2(0) and photon blockade in the driven, dissipative
Jaynes-Cummings model.

Steady state of
    d rho/dt = -i[H, rho] + kappa D[a] rho + gamma D[sigma-] rho,
    H = Delta a^dag a + (Delta/2) sigma_z + g (a^dag sigma- + a sigma+) + eps (a + a^dag),
with Delta = omega_c - omega_L = omega_a - omega_L (atom and cavity resonant,
laser detuned by Delta), in the frame rotating at the laser frequency.

    g2(0) = <a^dag a^dag a a> / <a^dag a>^2

IMPORTANT: the blockade signature (g2(0) << 1) appears when the laser is
resonant with a polariton, Delta = +-g.  Driving on BARE cavity resonance
(Delta = 0) in strong coupling is detuned from both polaritons by g while
the two-photon transition |0> -> |2,+-> is detuned by only g/sqrt(2); the
result is <n> ~ 1e-8 and g2(0) >> 1 (bunching), not antibunching.
Unless stated otherwise every curve below uses Delta = g (lower polariton).

Generates (in figures/): fig_g2_combined, fig_g2_vs_coupling, fig_g2_vs_drive,
    fig_g2_blockade_spectrum, g2_data.npz

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

from plotstyle import (use_style, save, panel_label, sequential_colors, FULL_W, HALF_W,
                       COLORS, FIG_DIR)

use_style()

kappa = 1.0          # cavity decay sets the unit
gamma = 0.1 * kappa  # atomic decay


def steady_g2(g, eps, delta, N):
    """Return (g2(0), <n>, P(N-1)) in the steady state."""
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    sz = qt.tensor(qt.qeye(N), qt.sigmaz())
    H = (delta * a.dag() * a + 0.5 * delta * sz
         + g * (a.dag() * sm + a * sm.dag()) + eps * (a + a.dag()))
    rho = qt.steadystate(H, [np.sqrt(kappa) * a, np.sqrt(gamma) * sm])
    n = qt.expect(a.dag() * a, rho).real
    g2 = qt.expect(a.dag() * a.dag() * a * a, rho).real / n ** 2 if n > 1e-14 else np.nan
    return g2, n, rho.ptrace(0).diag().real[-1]


def sweep(points, N, tag):
    """Evaluate steady_g2 over a list of (g, eps, delta); warn if the Fock space is too small."""
    out = np.array([steady_g2(g, e, d, N) for g, e, d in points])
    if out[:, 2].max() > 1e-5:
        print(f"  WARNING truncation in {tag}: max P(N-1) = {out[:, 2].max():.1e}")
    return out[:, 0], out[:, 1]


# (a) g2(0) vs g/kappa at Delta = g for several drives; Delta = 0 for contrast
print("(a) g2(0) vs g/kappa ...")
gk_grid = np.linspace(0.1, 15, 60)
eps_list = [0.01, 0.1, 0.5]
g2_a = np.array([sweep([(gk * kappa, e * kappa, gk * kappa) for gk in gk_grid], 20, f'(a) eps={e}')[0]
                 for e in eps_list])
n_a = np.array([sweep([(gk * kappa, e * kappa, gk * kappa) for gk in gk_grid], 20, f'(a) eps={e}')[1]
                for e in eps_list])
g2_a_bare = sweep([(gk * kappa, 0.05 * kappa, 0.0) for gk in gk_grid], 20, '(a) bare')[0]
g2_a_ref = sweep([(gk * kappa, 0.05 * kappa, gk * kappa) for gk in gk_grid], 20, '(a) ref')[0]

# (b) g2(0) vs drive strength at Delta = g
print("(b) g2(0) vs drive ...")
eps_grid = np.logspace(-2, 0.5, 30)
gk_list = [2, 5, 10]
g2_b = np.array([sweep([(gk * kappa, e * kappa, gk * kappa) for e in eps_grid], 60, f'(b) g/k={gk}')[0]
                 for gk in gk_list])

# (c, d) g2(0) and <n> vs laser detuning at g/kappa = 5
print("(c) spectrum ...")
gk_c, eps_c = 5.0, 0.05
g_c = gk_c * kappa
delta_grid = np.linspace(-3, 3, 481) * g_c
g2_c, n_c = sweep([(g_c, eps_c * kappa, d) for d in delta_grid], 15, '(c)')

np.savez(FIG_DIR / 'g2_data.npz', gk_grid=gk_grid, eps_list=np.array(eps_list), g2_vs_coupling=g2_a,
         n_vs_coupling=n_a, g2_bare_resonance=g2_a_bare, eps_grid=eps_grid,
         gk_list=np.array(gk_list), g2_vs_drive=g2_b, delta_grid=delta_grid,
         g2_vs_detuning=g2_c, n_vs_detuning=n_c)

STYLES = ['-', (0, (5, 2)), (0, (1, 1.5))]


def panel_a(ax):
    for e, y, c, ls in zip(eps_list, g2_a, sequential_colors(len(eps_list), 'viridis', 0, 0.75), STYLES):
        ax.semilogy(gk_grid, y, lw=1.1, color=c, ls=ls, label=rf'$\varepsilon/\kappa = {e:g}$')
    ax.semilogy(gk_grid, g2_a_bare, color=COLORS['grey'], ls=(0, (3, 1, 1, 1)), lw=0.9,
                label=r'$\Delta = 0$ (bare cavity)')
    ax.axhline(1, color='0.6', lw=0.4)
    ax.set_xlabel(r'$g/\kappa$'); ax.set_ylabel(r'$g^{(2)}(0)$')
    ax.set_xlim(0, 15); ax.set_ylim(1e-2, 1e3)
    ax.set_xticks([0, 5, 10, 15])
    ax.legend(loc='upper right', fontsize=7.5)


def panel_n(ax):
    for e, y, c, ls in zip(eps_list, n_a, sequential_colors(len(eps_list), 'viridis', 0, 0.75), STYLES):
        ax.semilogy(gk_grid, y, lw=1.1, color=c, ls=ls)
    ax.set_xlabel(r'$g/\kappa$'); ax.set_ylabel(r'$\langle n \rangle$')
    ax.set_xlim(0, 15)
    ax.set_xticks([0, 5, 10, 15])


def panel_b(ax):
    for gk, y, c, ls in zip(gk_list, g2_b, sequential_colors(len(gk_list), 'plasma', 0.1, 0.8), STYLES):
        ax.loglog(eps_grid, y, lw=1.1, color=c, ls=ls, label=rf'$g/\kappa = {gk}$')
    ax.axhline(1, color='0.6', lw=0.4)
    ax.set_xlabel(r'$\varepsilon/\kappa$'); ax.set_ylabel(r'$g^{(2)}(0)$')
    ax.legend(loc='lower right', fontsize=7.5)


def mark_resonances(ax, two_photon):
    for s in (-1, 1):
        ax.axvline(s, color=COLORS['red'], ls=':', lw=0.7)
        if two_photon:
            ax.axvline(s / np.sqrt(2), color=COLORS['orange'], ls='--', lw=0.7)


def panel_c(ax):
    ax.semilogy(delta_grid / g_c, g2_c, color=COLORS['blue'], lw=1.0)
    ax.axhline(1, color='0.6', lw=0.4)
    mark_resonances(ax, True)
    ax.plot([], [], color=COLORS['red'], ls=':', lw=0.7, label=r'polaritons, $\Delta = \pm g$')
    ax.plot([], [], color=COLORS['orange'], ls='--', lw=0.7, label=r'two-photon, $\Delta = \pm g/\sqrt{2}$')
    ax.set_xlabel(r'$\Delta/g$'); ax.set_ylabel(r'$g^{(2)}(0)$')
    ax.set_xlim(-3, 3); ax.set_ylim(3e-2, 1e10)
    ax.set_yticks([1e-1, 1e1, 1e3, 1e5, 1e7])
    ax.legend(loc='upper right', fontsize=7, handlelength=1.8, frameon=True, framealpha=0.92,
              edgecolor='none', borderpad=0.3)


def panel_d(ax):
    ax.plot(delta_grid / g_c, 1e3 * n_c, color=COLORS['blue'], lw=1.0)
    mark_resonances(ax, False)
    ax.set_xlabel(r'$\Delta/g$'); ax.set_ylabel(r'$10^3\,\langle n \rangle$')
    ax.set_xlim(-3, 3); ax.set_ylim(0, None)


fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 4.5))
panel_a(axes[0, 0]); panel_label(axes[0, 0], 'a', r'drive on the lower polariton, $\Delta = g$')
panel_b(axes[0, 1]); panel_label(axes[0, 1], 'b', r'$\Delta = g$')
panel_c(axes[1, 0]); panel_label(axes[1, 0], 'c', rf'$g/\kappa = {gk_c:g}$, $\varepsilon/\kappa = {eps_c:g}$')
panel_d(axes[1, 1]); panel_label(axes[1, 1], 'd', 'cavity photon number')
save(fig, FIG_DIR / 'fig_g2_combined')

fig, (axa, axn) = plt.subplots(2, 1, figsize=(HALF_W, 4.0), sharex=True, height_ratios=[1.5, 1])
panel_a(axa); axa.set_xlabel(''); panel_label(axa, 'a'); panel_n(axn); panel_label(axn, 'b')
save(fig, FIG_DIR / 'fig_g2_vs_coupling')

fig, ax = plt.subplots(figsize=(HALF_W, 2.5)); panel_b(ax)
save(fig, FIG_DIR / 'fig_g2_vs_drive')

fig, (axc, axd) = plt.subplots(2, 1, figsize=(HALF_W, 4.2), sharex=True)
panel_c(axc); axc.set_xlabel(''); panel_label(axc, 'a'); panel_d(axd); panel_label(axd, 'b')
save(fig, FIG_DIR / 'fig_g2_blockade_spectrum')

i5 = int(np.argmin(np.abs(gk_grid - 5))); i10 = int(np.argmin(np.abs(gk_grid - 10)))
print(f"\ng2(0) at Delta = g, eps/kappa = 0.05: g/kappa = {gk_grid[i5]:.2f} -> {g2_a_ref[i5]:.3f}, "
      f"g/kappa = {gk_grid[i10]:.2f} -> {g2_a_ref[i10]:.3f}")
for gk in (5.0, 10.0):
    print(f"  exactly g/kappa = {gk:g}: g2(0) = {steady_g2(gk * kappa, 0.05 * kappa, gk * kappa, 20)[0]:.3f}")
print(f"g2(0) at Delta = 0, eps/kappa = 0.05, g/kappa = {gk_grid[i5]:.2f}: {g2_a_bare[i5]:.3g}")
k = int(np.argmin(np.abs(delta_grid / g_c - 1)))
print(f"spectrum at Delta = g: g2 = {g2_c[k]:.3f}, <n> = {n_c[k]:.5f}; max g2 = {g2_c.max():.3g}")
