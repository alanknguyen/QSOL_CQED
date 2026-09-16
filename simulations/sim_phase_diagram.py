#!/usr/bin/env python3
"""
Cat-state survival phase diagram: Wigner negativity delta(n_bar, kappa/g) and
cat-state fidelity F(n_bar, kappa/g) of the reduced cavity field at the
half-revival time t = t_r/2 = pi sqrt(n_bar)/g.

    delta = int |W| dx dp - 1                  (Kenfack & Zyczkowski 2004)

Cat fidelity.  The JC half-revival cat for an initially excited atom has its
two coherent branches at +- i sqrt(n_bar) (rotated by 90 degrees from the
initial alpha = sqrt(n_bar)) with a relative phase that is NOT the even-cat
phase; the reduced state at t_r/2 is close to an ODD cat for n_bar = 10.
Comparing with the even cat |alpha> + |-alpha> on the real axis gives F = 0
for every n_bar >~ 3 and is meaningless.  We therefore report

    F = max_{beta, theta} <cat(beta, theta)| rho_field |cat(beta, theta)>,
    |cat(beta, theta)> = N (|beta> + e^{i theta} |-beta>),

with beta scanned in a small neighbourhood of i sqrt(n_bar) and theta
scanned over [0, 2 pi).

Generates (in figures/): fig_phase_diagram, fig_cat_fidelity,
    fig_phase_combined, fig_phase_slices, phase_diagram_data.npz

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import qutip as qt

FIG_DIR = Path(__file__).resolve().parents[1] / 'figures'
FIG_DIR.mkdir(exist_ok=True)

g = 1.0
nbar_grid = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25], dtype=float)
kappa_grid = np.linspace(0.0, 0.15, 16)          # step 0.01


def wigner_negativity(W, xvec):
    dx = xvec[1] - xvec[0]
    return float(np.sum(np.abs(W)) * dx * dx - 1.0)


def _scalar(q):
    return complex(q.full()[0, 0]) if isinstance(q, qt.Qobj) else complex(q)


def cat_fidelity(rho, nbar, N):
    """max over beta near i*sqrt(nbar) and relative phase theta of <cat|rho|cat>."""
    best = 0.0
    thetas = np.linspace(0, 2 * np.pi, 361)
    for scale in (0.9, 0.95, 1.0, 1.05, 1.1):
        for dphi in (-0.15, -0.075, 0.0, 0.075, 0.15):
            beta = scale * np.sqrt(nbar) * np.exp(1j * (np.pi / 2 + dphi))
            kp = qt.coherent(N, beta); km = qt.coherent(N, -beta)
            rpp = qt.expect(rho, kp).real
            rmm = qt.expect(rho, km).real
            c = _scalar(kp.dag() * rho * km)     # <beta| rho |-beta>
            ov = _scalar(kp.dag() * km)          # <beta|-beta> = exp(-2|beta|^2)
            num = rpp + rmm + 2 * np.real(np.exp(1j * thetas) * c)
            den = 2 + 2 * np.real(np.exp(1j * thetas) * ov)
            best = max(best, float(np.max(num / den)))
    return best


def half_revival_state(nbar, kappa):
    N = int(3 * nbar + 20)
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    H = g * (a.dag() * sm + a * sm.dag())
    psi0 = qt.tensor(qt.coherent(N, np.sqrt(nbar)), qt.basis(2, 0))
    t_cat = np.pi * np.sqrt(nbar) / g
    c_ops = [np.sqrt(kappa) * a] if kappa > 0 else []
    res = qt.mesolve(H, psi0, [0, t_cat], c_ops, [])
    st = res.states[-1]
    rho = st.ptrace(0)
    return rho, N


delta_map = np.zeros((len(kappa_grid), len(nbar_grid)))
fid_map = np.zeros_like(delta_map)
purity_map = np.zeros_like(delta_map)

t0 = time.time()
for j, nbar in enumerate(nbar_grid):
    L = np.sqrt(2 * nbar) + 4.5
    npts = int(min(400, max(200, 2 * L / 0.07)))
    xvec = np.linspace(-L, L, npts)
    for i, kap in enumerate(kappa_grid):
        rho, N = half_revival_state(nbar, kap)
        W = qt.wigner(rho, xvec, xvec)
        delta_map[i, j] = max(wigner_negativity(W, xvec), 0.0)
        fid_map[i, j] = cat_fidelity(rho, nbar, N)
        purity_map[i, j] = (rho ** 2).tr().real
    print(f"  n_bar = {nbar:4.0f} done  ({time.time() - t0:5.0f} s)  "
          f"delta(k=0) = {delta_map[0, j]:.3f}, F(k=0) = {fid_map[0, j]:.3f}, "
          f"delta(k=0.03) = {delta_map[3, j]:.3f}")

np.savez(FIG_DIR / 'phase_diagram_data.npz', nbar=nbar_grid, kappa=kappa_grid,
         delta=delta_map, fidelity=fid_map, purity=purity_map)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def draw_map(ax, Z, label, cmap, contours=(0.01, 0.05, 0.1, 0.2, 0.3), highlight=0.05,
             vmax=None):
    im = ax.pcolormesh(nbar_grid, kappa_grid, Z, cmap=cmap, shading='nearest',
                       vmin=0, vmax=vmax)
    cs = ax.contour(nbar_grid, kappa_grid, Z, levels=[c for c in contours if c != highlight],
                    colors='white', linewidths=0.8, linestyles='--')
    ax.clabel(cs, fontsize=7, fmt='%.2f')
    if highlight is not None:
        ax.contour(nbar_grid, kappa_grid, Z, levels=[highlight], colors='cyan', linewidths=2.2)
    ax.set_xlabel(r'$\bar{n}$', fontsize=12)
    ax.set_ylabel(r'$\kappa/g$', fontsize=12)
    plt.colorbar(im, ax=ax, label=label)


fig, ax = plt.subplots(figsize=(8, 6))
draw_map(ax, delta_map, r'$\delta$', 'inferno')
ax.set_title(r'(a) Wigner negativity $\delta$ at $t = t_r/2$ (cyan: $\delta = 0.05$)', fontsize=12)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_phase_diagram.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_phase_diagram.pdf', bbox_inches='tight'); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 6))
draw_map(ax, fid_map, r'$F_{\rm cat}$', 'viridis', contours=(0.3, 0.5, 0.7), highlight=0.5, vmax=1)
ax.set_title(r'(b) Cat-state fidelity $F = \max_{\beta,\theta}\langle{\rm cat}|\rho_{\rm field}|{\rm cat}\rangle$ (cyan: 0.5)', fontsize=11)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_cat_fidelity.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_cat_fidelity.pdf', bbox_inches='tight'); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
draw_map(axes[0], delta_map, r'$\delta$', 'inferno')
axes[0].set_title(r'(a) Wigner negativity $\delta$ (cyan: $\delta = 0.05$)', fontsize=12)
draw_map(axes[1], fid_map, r'$F_{\rm cat}$', 'viridis', contours=(0.3, 0.5, 0.7), highlight=0.5, vmax=1)
axes[1].set_title(r'(b) Cat-state fidelity $F_{\rm cat}$ (cyan: $F = 0.5$)', fontsize=12)
fig.suptitle(r'Cat-state survival in the dissipative Jaynes--Cummings model at $t = t_r/2$', fontsize=14)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_phase_combined.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_phase_combined.pdf', bbox_inches='tight'); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
for kap in (0.0, 0.03, 0.06, 0.10, 0.15):
    i = np.argmin(np.abs(kappa_grid - kap))
    axes[0].plot(nbar_grid, delta_map[i], 'o-', ms=4, label=rf'$\kappa/g = {kappa_grid[i]:.2f}$')
axes[0].axhline(0.05, color='gray', ls=':')
axes[0].set_xlabel(r'$\bar{n}$', fontsize=12); axes[0].set_ylabel(r'$\delta$', fontsize=12)
axes[0].set_title(r'Wigner negativity vs $\bar{n}$ (fixed $\kappa/g$)', fontsize=12)
axes[0].legend(fontsize=9)
for nb in (3, 5, 10, 14, 20):
    j = np.argmin(np.abs(nbar_grid - nb))
    axes[1].plot(kappa_grid, delta_map[:, j], 'o-', ms=4, label=rf'$\bar{{n}} = {nbar_grid[j]:.0f}$')
axes[1].set_xlabel(r'$\kappa/g$', fontsize=12); axes[1].set_ylabel(r'$\delta$', fontsize=12)
axes[1].set_title(r'Wigner negativity vs $\kappa/g$ (fixed $\bar{n}$)', fontsize=12)
axes[1].legend(fontsize=9)
fig.suptitle('Cat-state survival: parameter slices', fontsize=14)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_phase_slices.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_phase_slices.pdf', bbox_inches='tight'); plt.close(fig)

j10 = np.argmin(np.abs(nbar_grid - 10))
print(f"\nn_bar = 10: delta(k=0) = {delta_map[0, j10]:.3f}, F(k=0) = {fid_map[0, j10]:.3f}, "
      f"delta(k=0.02) = {delta_map[2, j10]:.3f}")
print("Saved: fig_phase_diagram, fig_cat_fidelity, fig_phase_combined, fig_phase_slices")
