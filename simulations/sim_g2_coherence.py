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
Unless stated otherwise every curve below uses Delta = g.

Generates (in figures/):
    fig_g2_vs_coupling, fig_g2_vs_drive, fig_g2_blockade_spectrum, fig_g2_combined

Author: Nguyen Khoi Nguyen (Alan), Boston University
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import qutip as qt

FIG_DIR = Path(__file__).resolve().parents[1] / 'figures'
FIG_DIR.mkdir(exist_ok=True)

kappa = 1.0          # cavity decay sets the unit
gamma = 0.1 * kappa  # atomic decay


def steady_g2(g, eps, delta, N, kappa=kappa, gamma=gamma):
    """Return (g2(0), <n>, P(N-1)) in the steady state."""
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    sz = qt.tensor(qt.qeye(N), qt.sigmaz())
    H = (delta * a.dag() * a + 0.5 * delta * sz
         + g * (a.dag() * sm + a * sm.dag()) + eps * (a + a.dag()))
    c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma) * sm]
    rho = qt.steadystate(H, c_ops)
    n = qt.expect(a.dag() * a, rho).real
    g2 = qt.expect(a.dag() * a.dag() * a * a, rho).real / n**2 if n > 1e-14 else np.nan
    p_top = rho.ptrace(0).diag().real[-1]
    return g2, n, p_top


def check_trunc(p_top, tag):
    if p_top > 1e-5:
        print(f"  WARNING truncation: {tag} P(N-1) = {p_top:.1e}")


# ---------------------------------------------------------------------------
# (a) g2(0) vs g/kappa at Delta = g for several drives; Delta = 0 for contrast
# ---------------------------------------------------------------------------
print("(a) g2(0) vs g/kappa ...")
N_a = 20
gk_grid = np.linspace(0.1, 15, 50)
eps_list = [0.01, 0.05, 0.1, 0.5]
g2_a = {e: [] for e in eps_list}
n_a = {e: [] for e in eps_list}
g2_a_res = []          # Delta = 0 reference for eps/kappa = 0.05
for gk in gk_grid:
    g = gk * kappa
    for e in eps_list:
        g2v, nv, pt = steady_g2(g, e * kappa, delta=g, N=N_a)
        check_trunc(pt, f"(a) g/k={gk:.1f} eps={e}")
        g2_a[e].append(g2v); n_a[e].append(nv)
    g2_a_res.append(steady_g2(g, 0.05 * kappa, delta=0.0, N=N_a)[0])

# ---------------------------------------------------------------------------
# (b) g2(0) vs drive strength at Delta = g
# ---------------------------------------------------------------------------
print("(b) g2(0) vs drive ...")
N_b = 60
eps_grid = np.logspace(-2, 0.5, 30)      # eps/kappa from 0.01 to ~3
gk_list = [2, 5, 10]
g2_b = {gk: [] for gk in gk_list}
n_b = {gk: [] for gk in gk_list}
for gk in gk_list:
    g = gk * kappa
    for e in eps_grid:
        g2v, nv, pt = steady_g2(g, e * kappa, delta=g, N=N_b)
        check_trunc(pt, f"(b) g/k={gk} eps={e:.2f}")
        g2_b[gk].append(g2v); n_b[gk].append(nv)

# ---------------------------------------------------------------------------
# (c) blockade spectrum: g2(0) and <n> vs Delta at g/kappa = 5
# ---------------------------------------------------------------------------
print("(c) spectrum ...")
N_c = 15
gk_c, eps_c = 5.0, 0.05
g_c = gk_c * kappa
delta_grid = np.linspace(-3, 3, 361) * g_c
g2_c, n_c = [], []
for d in delta_grid:
    g2v, nv, pt = steady_g2(g_c, eps_c * kappa, delta=d, N=N_c)
    g2_c.append(g2v); n_c.append(nv)
g2_c = np.array(g2_c); n_c = np.array(n_c)

np.savez(FIG_DIR / 'g2_data.npz', gk_grid=gk_grid, eps_list=eps_list,
         g2_a={e: np.array(v) for e, v in g2_a.items()}, g2_a_res=np.array(g2_a_res),
         eps_grid=eps_grid, gk_list=gk_list, g2_b={k: np.array(v) for k, v in g2_b.items()},
         delta_grid=delta_grid, g2_c=g2_c, n_c=n_c, allow_pickle=True)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def panel_a(ax, ax_n=None):
    for e in eps_list:
        ax.semilogy(gk_grid, g2_a[e], lw=2, label=rf'$\varepsilon/\kappa = {e}$')
    ax.semilogy(gk_grid, g2_a_res, color='gray', ls='--', lw=1.2,
                label=r'$\Delta = 0$ (bare resonance), $\varepsilon/\kappa = 0.05$')
    ax.axhline(1, color='k', ls=':', lw=0.8)
    ax.axvline(1, color='red', ls='--', alpha=0.4, lw=1.2)
    ax.set_xlabel(r'$g/\kappa$', fontsize=12)
    ax.set_ylabel(r'$g^{(2)}(0)$', fontsize=12)
    ax.set_title(r'(a) Blockade transition, drive on polariton $\Delta = g$', fontsize=11)
    ax.set_ylim(1e-2, 1e4)
    ax.legend(fontsize=8, loc='upper right')
    if ax_n is not None:
        for e in eps_list:
            ax_n.semilogy(gk_grid, n_a[e], lw=2, label=rf'$\varepsilon/\kappa = {e}$')
        ax_n.set_xlabel(r'$g/\kappa$', fontsize=12)
        ax_n.set_ylabel(r'$\langle n \rangle$', fontsize=12)
        ax_n.legend(fontsize=8)


def panel_b(ax):
    for gk in gk_list:
        ax.loglog(eps_grid, g2_b[gk], lw=2, label=rf'$g/\kappa = {gk}$')
    ax.axhline(1, color='k', ls=':', lw=0.8)
    ax.set_xlabel(r'$\varepsilon/\kappa$', fontsize=12)
    ax.set_ylabel(r'$g^{(2)}(0)$', fontsize=12)
    ax.set_title(r'(b) $g^{(2)}(0)$ vs drive strength, $\Delta = g$', fontsize=11)
    ax.legend(fontsize=9)


def panel_c(ax):
    ax.semilogy(delta_grid / g_c, g2_c, 'b-', lw=1.8)
    ax.axhline(1, color='k', ls=':', lw=0.8)
    for s in (-1, 1):
        ax.axvline(s, color='red', ls='--', alpha=0.4, lw=1.2)
        ax.axvline(s / np.sqrt(2), color='orange', ls='--', alpha=0.5, lw=1.2)
    ax.set_xlabel(r'$\Delta/g$', fontsize=12)
    ax.set_ylabel(r'$g^{(2)}(0)$', fontsize=12)
    ax.set_title(rf'(c) Blockade spectrum, $g/\kappa = {gk_c:.0f}$, $\varepsilon/\kappa = {eps_c}$', fontsize=11)
    ax.text(0.02, 0.95, 'red: polaritons $\\Delta = \\pm g$ (antibunching)\n'
            'orange: two-photon resonance $\\Delta = \\pm g/\\sqrt{2}$ (bunching)',
            transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))


def panel_d(ax):
    ax.plot(delta_grid / g_c, n_c, 'r-', lw=1.8)
    for s in (-1, 1):
        ax.axvline(s, color='red', ls='--', alpha=0.4, lw=1.2)
    ax.set_xlabel(r'$\Delta/g$', fontsize=12)
    ax.set_ylabel(r'$\langle n \rangle$', fontsize=12)
    ax.set_title(r'(d) Vacuum Rabi doublet in transmission', fontsize=11)


# individual figures
fig, (axa, axn) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                               gridspec_kw={'height_ratios': [2, 1]})
panel_a(axa, axn)
fig.suptitle(r'Photon blockade: $g^{(2)}(0)$ vs coupling strength ($\gamma/\kappa = 0.1$)', fontsize=13)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_g2_vs_coupling.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_g2_vs_coupling.pdf', bbox_inches='tight'); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 5)); panel_b(ax); fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_g2_vs_drive.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_g2_vs_drive.pdf', bbox_inches='tight'); plt.close(fig)

fig, (axc, axd) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
panel_c(axc); panel_d(axd); fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_g2_blockade_spectrum.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_g2_blockade_spectrum.pdf', bbox_inches='tight'); plt.close(fig)

fig, axes = plt.subplots(2, 2, figsize=(15, 11))
panel_a(axes[0, 0]); panel_b(axes[0, 1]); panel_c(axes[1, 0]); panel_d(axes[1, 1])
fig.suptitle('Photon blockade and second-order coherence in the driven JC model', fontsize=14)
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig_g2_combined.png', dpi=200, bbox_inches='tight')
fig.savefig(FIG_DIR / 'fig_g2_combined.pdf', bbox_inches='tight'); plt.close(fig)

i5 = np.argmin(np.abs(gk_grid - 5)); i10 = np.argmin(np.abs(gk_grid - 10))
print(f"\ng2(0) at Delta=g, eps/k=0.05: g/k=5 -> {g2_a[0.05][i5]:.3f}, g/k=10 -> {g2_a[0.05][i10]:.3f}")
print(f"g2(0) at Delta=0, eps/k=0.05: g/k=5 -> {g2_a_res[i5]:.3g}")
print("Saved: fig_g2_vs_coupling, fig_g2_vs_drive, fig_g2_blockade_spectrum, fig_g2_combined")
