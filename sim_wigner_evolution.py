"""
Wigner Function Time Evolution in the Jaynes-Cummings Model
============================================================
Simulates a two-level atom coupled to a single cavity mode.
Initial state: |e> ⊗ |α> (atom excited, field in coherent state).

Key physics: During the collapse of Rabi oscillations, the field state
splits into a superposition of two coherent components — a mesoscopic
Schrödinger cat state. The Wigner function develops negative regions,
a definitive signature of non-classicality.

Author: Nguyen Khoi Nguyen (Alan), Boston University
Dependencies: qutip >= 5.0, numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from qutip import (
    basis, tensor, destroy, sigmax, sigmaz, sigmam, sigmap,
    mesolve, wigner, ptrace, expect, coherent, Qobj
)

# ============================================================
# Physical parameters
# ============================================================
N_cav = 35              # Fock space truncation
g = 1.0                 # Vacuum Rabi coupling (sets time unit)
omega_c = 0.0           # Rotating frame → cavity freq = 0
omega_a = 0.0           # Resonant: Δ = ω_a - ω_c = 0
alpha = np.sqrt(10.0)   # Coherent state amplitude → <n> = 10
kappa = 0.0             # Cavity decay rate (set nonzero for dissipation study)
gamma = 0.0             # Spontaneous emission rate

# ============================================================
# Operators in joint atom ⊗ cavity Hilbert space
# ============================================================
# Cavity operators
a = tensor(destroy(N_cav), Qobj(np.eye(2)))       # annihilation
# Atom operators (2-level: |e>=|0>, |g>=|1> in qutip convention)
sm = tensor(Qobj(np.eye(N_cav)), sigmam())         # |g><e| lowering
sp = tensor(Qobj(np.eye(N_cav)), sigmap())         # |e><g| raising
sz = tensor(Qobj(np.eye(N_cav)), sigmaz())         # σ_z

# Jaynes-Cummings Hamiltonian (rotating wave approximation)
# H_JC = ω_c a†a + (ω_a/2)σ_z + g(a†σ⁻ + a σ⁺)
H_JC = omega_c * a.dag() * a + 0.5 * omega_a * sz + g * (a.dag() * sm + a * sp)

# Collapse operators for Lindblad master equation
c_ops = []
if kappa > 0:
    c_ops.append(np.sqrt(kappa) * a)       # cavity photon loss
if gamma > 0:
    c_ops.append(np.sqrt(gamma) * sm)      # spontaneous emission

# ============================================================
# Initial state: |excited> ⊗ |α>
# ============================================================
psi_atom_e = basis(2, 0)      # excited state |e>
psi_field_coh = coherent(N_cav, alpha)
psi0 = tensor(psi_field_coh, psi_atom_e)

# ============================================================
# Time evolution
# ============================================================
# Characteristic times (on resonance, coherent state with <n> = |α|²):
#   t_collapse ≈ 1/g              (dephasing time)
#   t_revival  ≈ 2π√<n>/g         (first revival)
n_bar = np.abs(alpha)**2
t_collapse = 1.0 / g
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# Dense time grid for inversion plot
tlist_dense = np.linspace(0, 1.5 * t_revival, 2000)

# Snapshot times for Wigner function
t_snapshots = [
    0.0,                          # t=0: initial coherent state
    0.5 * t_collapse,             # early oscillation
    2.0 * t_collapse,             # entering collapse
    0.5 * t_revival,              # mid-collapse (cat state!)
    0.75 * t_revival,             # approaching revival
    t_revival,                    # first revival
]
snapshot_labels = [
    r'$t = 0$' + '\n(coherent state)',
    r'$t = 0.5\,t_c$' + '\n(early Rabi)',
    r'$t = 2\,t_c$' + '\n(collapse onset)',
    r'$t = 0.5\,t_r$' + '\n(cat state)',
    r'$t = 0.75\,t_r$' + '\n(pre-revival)',
    r'$t = t_r$' + '\n(first revival)',
]

# Combine all times, solve once
all_times = sorted(set(list(tlist_dense) + t_snapshots))
tlist_all = np.array(all_times)

print("Solving Jaynes-Cummings dynamics...")
print(f"  <n> = {n_bar:.1f}, g = {g:.2f}")
print(f"  t_collapse ≈ {t_collapse:.3f}, t_revival ≈ {t_revival:.3f}")
print(f"  Fock space truncation N = {N_cav}")

result = mesolve(H_JC, psi0, tlist_all, c_ops, [sz],
                 options={'store_states': True})

# ============================================================
# Extract atomic inversion for dense time grid
# ============================================================
# Map all_times indices back to dense grid
inversion = np.array(result.expect[0])

# ============================================================
# Figure 1: Atomic inversion with snapshot markers
# ============================================================
fig1, ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(tlist_all * g, inversion, 'b-', linewidth=0.6, alpha=0.9)
ax1.set_xlabel(r'$gt$', fontsize=14)
ax1.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=14)
ax1.set_title(
    rf'Atomic Inversion: Jaynes–Cummings Model ($\bar{{n}} = {n_bar:.0f}$, $\Delta = 0$)',
    fontsize=13
)
ax1.set_ylim(-1.1, 1.1)
ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')

# Mark snapshot times
for i, ts in enumerate(t_snapshots):
    idx = np.argmin(np.abs(tlist_all - ts))
    color = plt.cm.plasma(i / len(t_snapshots))
    ax1.axvline(ts * g, color=color, linewidth=1.2, linestyle=':', alpha=0.7)
    ax1.plot(ts * g, inversion[idx], 'o', color=color, markersize=8, zorder=5)

# Annotate collapse and revival
ax1.annotate('collapse', xy=(2 * t_collapse * g, 0), fontsize=10,
             ha='center', va='bottom', color='gray',
             arrowprops=dict(arrowstyle='->', color='gray'),
             xytext=(2 * t_collapse * g, 0.5))
ax1.annotate('revival', xy=(t_revival * g, 0.5), fontsize=10,
             ha='center', va='bottom', color='gray',
             arrowprops=dict(arrowstyle='->', color='gray'),
             xytext=(t_revival * g, 0.9))

ax1.set_xlim(0, 1.5 * t_revival * g)
fig1.tight_layout()
fig1.savefig('/home/claude/fig_inversion_snapshots.png', dpi=200, bbox_inches='tight')
fig1.savefig('/home/claude/fig_inversion_snapshots.pdf', bbox_inches='tight')
print("Saved: fig_inversion_snapshots.png/pdf")

# ============================================================
# Figure 2: Wigner function snapshots of the cavity field
# ============================================================
xvec = np.linspace(-7, 7, 200)

fig2, axes = plt.subplots(2, 3, figsize=(14, 9))
axes_flat = axes.flatten()

print("Computing Wigner functions at snapshot times...")
for i, ts in enumerate(t_snapshots):
    idx = np.argmin(np.abs(tlist_all - ts))
    # Get the full state at this time
    psi_t = result.states[idx]
    # Partial trace over atom → reduced density matrix of cavity field
    rho_field = ptrace(psi_t, 0)  # trace out atom (index 1), keep field (index 0)
    # Compute Wigner function
    W = wigner(rho_field, xvec, xvec)

    ax = axes_flat[i]
    wlim = np.max(np.abs(W))
    im = ax.contourf(xvec, xvec, W, levels=100, cmap='RdBu_r',
                      vmin=-wlim, vmax=wlim)
    ax.set_xlabel(r'$x$', fontsize=11)
    ax.set_ylabel(r'$p$', fontsize=11)
    ax.set_title(snapshot_labels[i], fontsize=10)
    ax.set_aspect('equal')
    fig2.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Report negativity
    W_neg = np.sum(W[W < 0]) * (xvec[1] - xvec[0])**2
    purity = (rho_field * rho_field).tr()
    print(f"  t={ts:.3f}: Wigner negativity volume = {abs(W_neg):.4f}, "
          f"field purity Tr(ρ²) = {purity:.4f}")

fig2.suptitle(
    r'Wigner Function Evolution of Cavity Field — Jaynes–Cummings Model'
    + f'\n' + rf'($\bar{{n}}={n_bar:.0f}$, $g={g}$, resonant)',
    fontsize=13, y=1.02
)
fig2.tight_layout()
fig2.savefig('/home/claude/fig_wigner_evolution.png', dpi=200, bbox_inches='tight')
fig2.savefig('/home/claude/fig_wigner_evolution.pdf', bbox_inches='tight')
print("Saved: fig_wigner_evolution.png/pdf")

# ============================================================
# Figure 3: Wigner cross-sections at cat-state time
# ============================================================
cat_idx = np.argmin(np.abs(tlist_all - 0.5 * t_revival))
psi_cat = result.states[cat_idx]
rho_cat = ptrace(psi_cat, 0)
W_cat = wigner(rho_cat, xvec, xvec)

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5))

# 2D contour
wlim = np.max(np.abs(W_cat))
im3 = ax3a.contourf(xvec, xvec, W_cat, levels=100, cmap='RdBu_r',
                      vmin=-wlim, vmax=wlim)
ax3a.set_xlabel(r'$x$', fontsize=13)
ax3a.set_ylabel(r'$p$', fontsize=13)
ax3a.set_title(r'Wigner function at $t = t_r/2$ (cat state)', fontsize=12)
ax3a.set_aspect('equal')
fig3.colorbar(im3, ax=ax3a)

# Cross-section through p=0
mid_idx = len(xvec) // 2
ax3b.plot(xvec, W_cat[mid_idx, :], 'b-', linewidth=1.5, label=r'$W(x, p=0)$')
ax3b.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax3b.fill_between(xvec, W_cat[mid_idx, :], 0,
                   where=(W_cat[mid_idx, :] < 0),
                   color='red', alpha=0.3, label='Negative region')
ax3b.set_xlabel(r'$x$', fontsize=13)
ax3b.set_ylabel(r'$W(x, 0)$', fontsize=13)
ax3b.set_title(r'Cross-section at $p = 0$: negative regions confirm non-classicality',
               fontsize=11)
ax3b.legend(fontsize=11)

fig3.tight_layout()
fig3.savefig('/home/claude/fig_cat_state_detail.png', dpi=200, bbox_inches='tight')
fig3.savefig('/home/claude/fig_cat_state_detail.pdf', bbox_inches='tight')
print("Saved: fig_cat_state_detail.png/pdf")

# ============================================================
# Figure 4: Dissipative JC — Wigner at cat time for various κ
# ============================================================
kappa_values = [0.0, 0.02 * g, 0.05 * g, 0.1 * g]

fig4, axes4 = plt.subplots(1, 4, figsize=(16, 4))

print("\nComputing dissipative Wigner functions...")
for j, kap in enumerate(kappa_values):
    c_ops_diss = [np.sqrt(kap) * a] if kap > 0 else []
    tlist_short = np.linspace(0, 0.5 * t_revival, 500)
    res_diss = mesolve(H_JC, psi0, tlist_short, c_ops_diss, [])
    rho_f = ptrace(res_diss.states[-1], 0)
    W_d = wigner(rho_f, xvec, xvec)
    wlim_d = max(np.max(np.abs(W_d)), 1e-10)

    im4 = axes4[j].contourf(xvec, xvec, W_d, levels=80, cmap='RdBu_r',
                              vmin=-wlim_d, vmax=wlim_d)
    axes4[j].set_title(rf'$\kappa/g = {kap/g:.2f}$', fontsize=11)
    axes4[j].set_xlabel(r'$x$', fontsize=10)
    if j == 0:
        axes4[j].set_ylabel(r'$p$', fontsize=10)
    axes4[j].set_aspect('equal')
    fig4.colorbar(im4, ax=axes4[j], fraction=0.046, pad=0.04)

    W_neg_d = np.sum(W_d[W_d < 0]) * (xvec[1] - xvec[0])**2
    purity_d = (rho_f * rho_f).tr()
    print(f"  κ/g = {kap/g:.2f}: negativity = {abs(W_neg_d):.4f}, purity = {purity_d:.4f}")

fig4.suptitle(
    rf'Decoherence of Schrödinger Cat State at $t = t_r/2$ ($\bar{{n}}={n_bar:.0f}$)',
    fontsize=13, y=1.04
)
fig4.tight_layout()
fig4.savefig('/home/claude/fig_wigner_decoherence.png', dpi=200, bbox_inches='tight')
fig4.savefig('/home/claude/fig_wigner_decoherence.pdf', bbox_inches='tight')
print("Saved: fig_wigner_decoherence.png/pdf")

print("\n=== Wigner evolution simulations complete ===")
