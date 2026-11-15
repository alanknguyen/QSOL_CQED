#!/usr/bin/env python3
"""
Cat-state survival in the dissipative Jaynes-Cummings model: Wigner negativity
delta(n_bar, kappa/g) and cat-state fidelity F(n_bar, kappa/g) of the reduced
cavity field at the half-revival time t = t_r/2 = pi sqrt(n_bar)/g.

    delta = int |W| dx dp - 1                  (Kenfack & Zyczkowski 2004)
    F     = max_{beta, theta} <cat(beta, theta)| rho_field |cat(beta, theta)>,
            |cat(beta, theta)> ~ |beta> + e^{i theta} |-beta>,

with beta scanned near i sqrt(n_bar): the two components of the half-revival
state sit on the p axis.  Their relative phase depends on n_bar (the
photon-number parity at t_r/2 is about -0.64 for even integer n_bar and +0.64
for odd), so theta is scanned over [0, 2 pi) rather than fixed.

kappa is sampled on a LOGARITHMIC grid.  delta falls roughly as
exp(-pi kappa n_bar^{3/2} / g) (the fringe-decoherence rate kappa |Delta alpha|^2 / 2
integrated over the time the two components take to separate), so the
delta = 0.05 boundary moves to kappa/g < 0.01 for n_bar > 20 and a linear grid
with step 0.01 cannot resolve it.

Generates (in figures/): fig_phase_combined, fig_phase_diagram, fig_cat_fidelity,
    fig_phase_slices, phase_diagram_data.npz

The sweep (16 x 25 Lindblad runs) takes about six minutes.  To redraw the
figures from the saved data without recomputing:  python sim_phase_diagram.py --replot

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import qutip as qt

from plotstyle import use_style, save, panel_label, sequential_colors, FULL_W, HALF_W, FIG_DIR
from utils import wigner_negativity, cat_fidelity

use_style()

g = 1.0
THRESHOLD = 0.05                                   # delta below this: cat no longer visible
nbar_grid = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25], dtype=float)
kappa_grid = np.concatenate([[0.0], np.geomspace(1e-3, 0.2, 24)])


def half_revival_field(nbar, kappa):
    """Reduced field state at t_r/2 for |e> x |alpha = sqrt(nbar)> with cavity loss kappa."""
    N = int(3 * nbar + 20)
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    H = g * (a.dag() * sm + a * sm.dag())
    psi0 = qt.tensor(qt.coherent(N, np.sqrt(nbar)), qt.basis(2, 0))
    c_ops = [np.sqrt(kappa) * a] if kappa > 0 else []
    res = qt.mesolve(H, psi0, [0, np.pi * np.sqrt(nbar) / g], c_ops=c_ops)
    return res.states[-1].ptrace(0)


DATA = FIG_DIR / 'phase_diagram_data.npz'
if '--replot' in sys.argv and DATA.exists():
    d = np.load(DATA)
    nbar_grid, kappa_grid = d['nbar'], d['kappa']
    delta_map, fid_map, purity_map = d['delta'], d['fidelity'], d['purity']
    print(f"Loaded {DATA.name}; redrawing only.")
else:
    delta_map = np.zeros((len(kappa_grid), len(nbar_grid)))
    fid_map = np.zeros_like(delta_map)
    purity_map = np.zeros_like(delta_map)
    t0 = time.time()
    for j, nbar in enumerate(nbar_grid):
        L = np.sqrt(2 * nbar) + 4.5                    # phase-space window
        xvec = np.linspace(-L, L, int(2 * L / 0.055) | 1)     # ~13 points per fringe at n_bar = 10, ~8 at n_bar = 25
        for i, kap in enumerate(kappa_grid):
            rho = half_revival_field(nbar, kap)
            delta_map[i, j] = max(wigner_negativity(qt.wigner(rho, xvec, xvec), xvec), 0.0)
            fid_map[i, j] = cat_fidelity(rho, nbar)
            purity_map[i, j] = (rho * rho).tr().real
        print(f"  n_bar = {nbar:4.0f} done ({time.time() - t0:4.0f} s): delta(0) = {delta_map[0, j]:.3f}, "
              f"F(0) = {fid_map[0, j]:.3f}, purity(0) = {purity_map[0, j]:.3f}", flush=True)


def boundary(Z, level):
    """kappa at which Z(kappa) first drops below `level`, per n_bar (log-linear interpolation)."""
    out = np.full(len(nbar_grid), np.nan)
    for j in range(len(nbar_grid)):
        col = Z[:, j]
        below = np.where(col < level)[0]
        if len(below) == 0 or below[0] == 0:
            continue
        i = below[0]
        if kappa_grid[i - 1] == 0:                  # crossing inside the first (linear) interval
            out[j] = kappa_grid[i] * (col[i - 1] - level) / (col[i - 1] - col[i])
        else:
            w = (np.log(col[i - 1]) - np.log(level)) / (np.log(col[i - 1]) - np.log(col[i]))
            out[j] = np.exp(np.log(kappa_grid[i - 1]) + w * (np.log(kappa_grid[i]) - np.log(kappa_grid[i - 1])))
    return out


kc_delta = boundary(delta_map, THRESHOLD)
kc_fid = boundary(fid_map, 0.5)
# analytic estimate: delta ~ delta_0 exp(-pi kappa n^{3/2} / g)
kc_model = np.log(delta_map[0] / THRESHOLD) / (np.pi * nbar_grid ** 1.5)

np.savez(DATA, nbar=nbar_grid, kappa=kappa_grid,
         delta=delta_map, fidelity=fid_map, purity=purity_map,
         kappa_c_delta=kc_delta, kappa_c_fidelity=kc_fid)

print("\nBoundary delta = 0.05 (kappa_c/g), with the analytic estimate ln(delta_0/0.05)/(pi n^1.5):")
for nb, kc, km, kf in zip(nbar_grid, kc_delta, kc_model, kc_fid):
    print(f"  n_bar = {nb:4.0f}: kappa_c = {kc:.4f}   estimate = {km:.4f}   (F = 0.5 at {kf:.4f})")
for lo in (3, 10):
    m = (nbar_grid >= lo) & np.isfinite(kc_delta)
    slope, icpt = np.polyfit(np.log(nbar_grid[m]), np.log(kc_delta[m]), 1)
    print(f"  power-law fit for n_bar >= {lo}: kappa_c/g = {np.exp(icpt):.2f} n_bar^{slope:.2f}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
pos = kappa_grid > 0


def draw_map(ax, Z, kc, cmap, vmax, label, level_text, model=None):
    im = ax.pcolormesh(nbar_grid, kappa_grid[pos], Z[pos], cmap=cmap, shading='nearest',
                       vmin=0, vmax=vmax, rasterized=True)
    ax.plot(nbar_grid, kc, color='w', lw=1.6, label=level_text)
    if model is not None:
        ax.plot(nbar_grid[1:], model[1:], color='w', lw=1.0, ls=(0, (4, 2)),
                label=r'$\ln(\delta_0/0.05)\,/\,\pi\bar{n}^{3/2}$')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(1.5, 28); ax.set_ylim(kappa_grid[1] / 1.12, kappa_grid[-1] * 1.12)     # whole cells at the edges
    ax.set_xticks([2, 5, 10, 20]); ax.set_xticklabels(['2', '5', '10', '20'])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel(r'$\bar{n}$'); ax.set_ylabel(r'$\kappa/g$')
    leg = ax.legend(loc='upper right', fontsize=7.5, labelcolor='w', handlelength=1.8)      # the dark corner of both maps
    cb = plt.colorbar(im, ax=ax, pad=0.02, aspect=22)
    cb.set_label(label)


def delta_panel(ax):
    draw_map(ax, delta_map, kc_delta, 'inferno', None, r'$\delta$', r'$\delta = 0.05$', model=kc_model)


def fid_panel(ax):
    draw_map(ax, fid_map, kc_fid, 'viridis', 1.0, r'$F_{\rm cat}$', r'$F_{\rm cat} = 0.5$')


fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 2.6))
delta_panel(axes[0]); panel_label(axes[0], 'a', r'Wigner negativity at $t_r/2$')
fid_panel(axes[1]);   panel_label(axes[1], 'b', r'cat fidelity at $t_r/2$')
save(fig, FIG_DIR / 'fig_phase_combined')

fig, ax = plt.subplots(figsize=(HALF_W, 2.7)); delta_panel(ax); save(fig, FIG_DIR / 'fig_phase_diagram')
fig, ax = plt.subplots(figsize=(HALF_W, 2.7)); fid_panel(ax); save(fig, FIG_DIR / 'fig_cat_fidelity')

# -- slices, scaling collapse and the boundary --
fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(FULL_W, 2.25))
ks = [0, 5, 10, 14, 18]                             # kappa indices: 0, ~0.002, ~0.006, ~0.016, ~0.04
for i, c in zip(ks, ['k'] + sequential_colors(len(ks) - 1, 'plasma', 0.1, 0.85)):
    axa.plot(nbar_grid, delta_map[i], marker='o', ms=2.5, lw=0.9, color=c,
             label='0' if i == 0 else f'{kappa_grid[i]:.3f}')
axa.axhline(THRESHOLD, color='0.5', lw=0.5, ls=':')
axa.set_xlabel(r'$\bar{n}$'); axa.set_ylabel(r'$\delta$')
axa.set_ylim(0, 1.32)
axa.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axa.legend(title=r'$\kappa/g$', fontsize=6.5, title_fontsize=7, loc='upper right', ncol=3,
           labelspacing=0.2, columnspacing=0.7, handlelength=1.2, handletextpad=0.3, borderaxespad=0.2)
panel_label(axa, 'a')

x = np.linspace(0, 6, 200)
axb.semilogy(x, np.exp(-x), color='k', lw=0.8, ls='--', label=r'$e^{-\pi\kappa\bar{n}^{3/2}/g}$')
ns = [4, 6, 10, 16, 25]
for nb, c, mk in zip(ns, sequential_colors(len(ns), 'viridis', 0.0, 0.85), 'osD^v'):
    j = int(np.argmin(np.abs(nbar_grid - nb)))
    xx = np.pi * kappa_grid[pos] * nb ** 1.5
    axb.semilogy(xx, delta_map[pos, j] / delta_map[0, j], marker=mk, ms=2.8, lw=0.6, color=c,
                 label=rf'$\bar{{n}} = {nb}$')
axb.set_xlim(0, 6); axb.set_ylim(3e-3, 1.3)
axb.set_xlabel(r'$\pi\kappa\bar{n}^{3/2}/g$'); axb.set_ylabel(r'$\delta/\delta_0$')
axb.legend(fontsize=6.5, loc='upper right', labelspacing=0.2, handlelength=1.6, borderaxespad=0.2)
panel_label(axb, 'b')

axc.loglog(nbar_grid, kc_delta, 'o', ms=3, color='k', label=r'$\delta = 0.05$')
axc.loglog(nbar_grid[2:], kc_model[2:], color='k', lw=0.8, ls='--', label='estimate')
axc.loglog(nbar_grid, kc_fid, 's', ms=3, mfc='none', color='0.4', label=r'$F_{\rm cat} = 0.5$')
axc.set_xlim(1.5, 30)
axc.set_xticks([2, 5, 10, 20]); axc.set_xticklabels(['2', '5', '10', '20'])
axc.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axc.set_xlabel(r'$\bar{n}$'); axc.set_ylabel(r'$\kappa_c/g$')
axc.legend(fontsize=7, loc='lower left', labelspacing=0.25)
panel_label(axc, 'c')
save(fig, FIG_DIR / 'fig_phase_slices')

j10 = int(np.argmin(np.abs(nbar_grid - 10)))
print(f"\nn_bar = 10: delta(0) = {delta_map[0, j10]:.3f}, F(0) = {fid_map[0, j10]:.3f}, "
      f"kappa_c = {kc_delta[j10]:.4f}")
