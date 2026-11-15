#!/usr/bin/env python3
"""
Animation: a tour of the phase diagram of the half-revival cat state.

Left: the Wigner negativity delta(n_bar, kappa/g) of the reduced cavity field at
t = t_r/2 (the stored sweep of sim_phase_diagram.py), with the simulated boundary
delta = 0.05 (solid white) and the estimate kappa_c = ln(delta_0/0.05)/(pi n_bar^{3/2})
(dashed).  A marker walks two legs:
  1. n_bar = 10, kappa/g from 0.001 to 0.1: at a fixed size, loss erases the fringes
     (boundary at kappa/g = 0.033);
  2. kappa/g = 0.01, n_bar from 2 to 25: at a FIXED loss rate a larger cat is more
     fragile, because its two components are further apart; it crosses the boundary
     near n_bar = 23.
Right: the Wigner function at the marker, recomputed for every frame on one fixed
window, so that the cat is seen to grow along leg 2.

Needs figures/phase_diagram_data.npz (run sim_phase_diagram.py first).
Generates: animations/anim_phase_diagram.gif
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import use_style, wigner_image, fig_to_frame, save_gif, FIG_DIR, ANIM_DIR
from utils import wigner_negativity

use_style()
g = 1.0
THRESHOLD = 0.05

d = np.load(FIG_DIR / 'phase_diagram_data.npz')
nbar_grid, kappa_grid, delta_map, kc_sim = d['nbar'], d['kappa'], d['delta'], d['kappa_c_delta']
pos = kappa_grid > 0


def half_revival_field(nbar, kappa):
    """Reduced field at t_r/2 for |e> x |alpha = sqrt(nbar)> with cavity loss kappa (as in sim_phase_diagram.py)."""
    N = int(3 * nbar + 20)
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    H = g * (a.dag() * sm + a * sm.dag())
    psi0 = qt.tensor(qt.coherent(N, np.sqrt(nbar)), qt.basis(2, 0))
    res = qt.mesolve(H, psi0, [0, np.pi * np.sqrt(nbar) / g], c_ops=[np.sqrt(kappa) * a])
    return res.states[-1].ptrace(0)


leg1 = [(10.0, k) for k in np.geomspace(1e-3, 0.1, 21)]
leg2 = [(float(n), 0.01) for n in range(2, 26)]
LEGS = [(leg1, '1: more loss, same size'), (leg2, '2: larger cat, same loss')]

L = np.sqrt(2 * 25) + 3.8                               # one window for every frame
xvec = np.linspace(-L, L, int(2 * L / 0.055) | 1)       # 8 points per fringe at n_bar = 25, more below

print("Solving the master equation along the two legs ...")
points = []
for leg, label in LEGS:
    for nbar, kap in leg:
        rho = half_revival_field(nbar, kap)
        W = qt.wigner(rho, xvec, xvec)
        points.append(dict(nbar=nbar, kappa=kap, W=W, label=label, delta=max(wigner_negativity(W, xvec), 0.0),
                           purity=(rho * rho).tr().real))
    print(f"  leg '{label}' done: delta from {points[-len(leg)]['delta']:.3f} to {points[-1]['delta']:.3f}", flush=True)

# the frames must agree with the stored map where the two overlap (n_bar = 10, kappa on the sweep grid)
j10 = int(np.argmin(np.abs(nbar_grid - 10)))
for pt in points[:len(leg1)]:
    i = int(np.argmin(np.abs(np.log(kappa_grid[pos]) - np.log(pt['kappa']))))
    if abs(np.log(kappa_grid[pos][i] / pt['kappa'])) < 1e-9:
        assert abs(delta_map[pos][i, j10] - pt['delta']) < 5e-3, (pt['kappa'], delta_map[pos][i, j10], pt['delta'])

nb_line = np.geomspace(1.5, 28, 100)
frames, durations = [], []
print("Rendering ...")
for k, pt in enumerate(points):
    fig, (ax_m, ax_w) = plt.subplots(1, 2, figsize=(8.2, 3.5), gridspec_kw={'width_ratios': [1.12, 1]})
    fig.get_layout_engine().set(wspace=0.08)             # room between the map's colour bar and the p axis
    im = ax_m.pcolormesh(nbar_grid, kappa_grid[pos], delta_map[pos], cmap='inferno', shading='nearest',
                         vmin=0, vmax=0.9, rasterized=True)
    ax_m.plot(nbar_grid, kc_sim, color='w', lw=1.3, label=rf'$\delta = {THRESHOLD}$')
    ax_m.plot(nb_line, np.log(0.85 / THRESHOLD) / (np.pi * nb_line ** 1.5), color='w', lw=0.9, ls='--',
              label=r'$\ln(\delta_0/0.05)/\pi\bar{n}^{3/2}$')
    for leg, _ in LEGS:                                   # the two paths
        ax_m.plot([p[0] for p in leg], [p[1] for p in leg], color='#7fdbff', lw=0.8, alpha=0.9)
    ax_m.plot(pt['nbar'], pt['kappa'], 'o', ms=7, mfc='#7fdbff', mec='k', mew=0.8, zorder=6)
    ax_m.set_xscale('log'); ax_m.set_yscale('log')
    ax_m.set_xlim(1.5, 28); ax_m.set_ylim(kappa_grid[1] / 1.12, kappa_grid[-1] * 1.12)
    ax_m.set_xticks([2, 5, 10, 20]); ax_m.set_xticklabels(['2', '5', '10', '20']); ax_m.minorticks_off()
    ax_m.set_xlabel(r'$\bar{n}$'); ax_m.set_ylabel(r'$\kappa/g$')
    ax_m.legend(loc='upper right', fontsize=7, labelcolor='w', handlelength=1.6, borderaxespad=0.3)
    fig.colorbar(im, ax=ax_m, pad=0.02, aspect=24, ticks=[0, 0.3, 0.6, 0.9], label=r'$\delta$')
    ax_m.set_title(rf"Wigner negativity at $t_r/2$; leg {pt['label']}", fontsize=9, loc='left')

    imw = wigner_image(ax_w, pt['W'], xvec)
    ax_w.set_xlabel(r'$x$'); ax_w.set_ylabel(r'$p$')
    ax_w.set_xticks([-10, -5, 0, 5, 10]); ax_w.set_yticks([-10, -5, 0, 5, 10])
    fig.colorbar(imw, ax=ax_w, shrink=0.9, aspect=22, pad=0.03, ticks=[-0.3, -0.15, 0, 0.15, 0.3], label=r'$W(x,p)$')
    visible = 'cat visible' if pt['delta'] > THRESHOLD else 'fringes lost'
    ax_w.set_title(rf"$\bar{{n}} = {pt['nbar']:.0f}$, $\kappa/g = {pt['kappa']:.4f}$: "
                   rf"$\delta = {pt['delta']:.3f}$ ({visible})", fontsize=9)
    frames.append(fig_to_frame(fig, dpi=115))
    first_of_leg = k in (0, len(leg1))
    last_of_leg = k in (len(leg1) - 1, len(points) - 1)
    durations.append(1600 if first_of_leg else (2000 if last_of_leg else 260))
save_gif(frames, ANIM_DIR / 'anim_phase_diagram.gif', durations)
