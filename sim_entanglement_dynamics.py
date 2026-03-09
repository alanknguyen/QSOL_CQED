"""
Entanglement Dynamics in the Jaynes-Cummings Model
====================================================
Computes the von Neumann entropy of the reduced atomic state S(ρ_atom)
as a measure of atom-field entanglement during JC evolution.

Key result: Entanglement entropy peaks during the collapse period
(when the atom and field are maximally entangled / field is in a cat state)
and drops at revival times (when the system approximately refactorizes).

We compare coherent, thermal, and squeezed initial field states to show
how the field's quantum statistics shape entanglement dynamics.

Author: Nguyen Khoi Nguyen (Alan), Boston University
Dependencies: qutip >= 5.0, numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap,
    mesolve, ptrace, entropy_vn, coherent, thermal_dm,
    squeeze, Qobj, ket2dm
)

# ============================================================
# System parameters
# ============================================================
N_cav = 50           # Fock space truncation (needs to be larger for thermal)
g = 1.0              # Vacuum Rabi coupling
n_bar = 10.0         # Mean photon number for all initial states
alpha = np.sqrt(n_bar)

# Operators
a = tensor(destroy(N_cav), Qobj(np.eye(2)))
sm = tensor(Qobj(np.eye(N_cav)), sigmam())
sp = tensor(Qobj(np.eye(N_cav)), sigmap())
sz = tensor(Qobj(np.eye(N_cav)), sigmaz())

# JC Hamiltonian (resonant, rotating frame)
H_JC = g * (a.dag() * sm + a * sp)

# ============================================================
# Initial states: atom excited, field in various states
# ============================================================
psi_e = basis(2, 0)  # |excited>
rho_e = ket2dm(psi_e)

# 1) Coherent state |α>
rho_coh = ket2dm(coherent(N_cav, alpha))

# 2) Thermal state with same <n>
rho_th = thermal_dm(N_cav, n_bar)

# 3) Squeezed vacuum: <n> = sinh²(r) → r = arcsinh(√n_bar)
r_sq = np.arcsinh(np.sqrt(n_bar))
psi_sq = squeeze(N_cav, r_sq) * basis(N_cav, 0)
rho_sq = ket2dm(psi_sq)

# 4) Fock state |n=10>  (exact photon number, for comparison)
rho_fock = ket2dm(basis(N_cav, int(n_bar)))

initial_states = {
    'Coherent': rho_coh,
    'Thermal': rho_th,
    'Squeezed vacuum': rho_sq,
    'Fock $|10\\rangle$': rho_fock,
}

# Time grid — go to ~1.5 revival times
t_revival = 2 * np.pi * np.sqrt(n_bar) / g
tlist = np.linspace(0, 1.8 * t_revival, 2000)

# ============================================================
# Solve and compute entanglement entropy
# ============================================================
results = {}

for label, rho_field in initial_states.items():
    print(f"Solving JC dynamics for {label} state...")

    # Joint initial state
    rho0 = tensor(rho_field, rho_e)

    # Solve master equation (no dissipation for now)
    result = mesolve(H_JC, rho0, tlist, [], [sz], options={'store_states': True})

    # Compute von Neumann entropy of reduced atomic state at each time
    S_atom = np.zeros(len(tlist))
    purity_field = np.zeros(len(tlist))
    inversion = np.array(result.expect[0])

    for i, t in enumerate(tlist):
        rho_t = result.states[i]
        # Reduced state of atom (trace out field, index 0)
        rho_atom = ptrace(rho_t, 1)
        S_atom[i] = entropy_vn(rho_atom, 2)  # base-2 entropy (bits)
        # Also track field purity
        rho_f = ptrace(rho_t, 0)
        purity_field[i] = (rho_f * rho_f).tr().real

    results[label] = {
        'inversion': inversion,
        'entropy': S_atom,
        'purity': purity_field,
    }
    print(f"  Max entropy: {np.max(S_atom):.4f} bits, "
          f"Min field purity: {np.min(purity_field):.4f}")

# ============================================================
# Figure 1: Entanglement entropy comparison
# ============================================================
colors = {'Coherent': '#2166ac', 'Thermal': '#b2182b',
          'Squeezed vacuum': '#1b7837', 'Fock $|10\\rangle$': '#e08214'}

fig1, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

# Top: inversion
ax1 = axes[0]
for label, data in results.items():
    ax1.plot(tlist * g, data['inversion'], color=colors[label],
             linewidth=0.7, alpha=0.85, label=label)
ax1.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax1.set_title(rf'Jaynes–Cummings Dynamics: $\bar{{n}} = {n_bar:.0f}$, $\Delta = 0$',
              fontsize=13)
ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax1.set_ylim(-1.1, 1.1)
ax1.legend(fontsize=10, loc='upper right', ncol=2)

# Bottom: entropy
ax2 = axes[1]
for label, data in results.items():
    ax2.plot(tlist * g, data['entropy'], color=colors[label],
             linewidth=1.0, alpha=0.9, label=label)
ax2.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax2.set_xlabel(r'$gt$', fontsize=13)
ax2.set_ylim(-0.05, 1.15)
ax2.axhline(1.0, color='gray', linewidth=0.5, linestyle=':', label=r'$S_{\max}=1$ bit')
ax2.legend(fontsize=10, loc='upper right', ncol=2)

# Mark revival time
for ax in axes:
    ax.axvline(t_revival * g, color='purple', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.axvline(2 * t_revival * g * 0.5, color='orange', linewidth=0.8,
               linestyle='--', alpha=0.4)

fig1.tight_layout()
fig1.savefig('/home/claude/fig_entanglement_comparison.png', dpi=200, bbox_inches='tight')
fig1.savefig('/home/claude/fig_entanglement_comparison.pdf', bbox_inches='tight')
print("\nSaved: fig_entanglement_comparison.png/pdf")

# ============================================================
# Figure 2: Coherent state — entropy + inversion + purity
# ============================================================
fig2, (ax2a, ax2b, ax2c) = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
data_coh = results['Coherent']

# Inversion
ax2a.plot(tlist * g, data_coh['inversion'], 'b-', linewidth=0.7)
ax2a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax2a.set_title(rf'Coherent State Field ($\alpha = {alpha:.1f}$, $\bar{{n}} = {n_bar:.0f}$)',
               fontsize=13)
ax2a.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax2a.set_ylim(-1.1, 1.1)

# Entropy
ax2b.plot(tlist * g, data_coh['entropy'], 'r-', linewidth=1.0)
ax2b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax2b.set_ylim(-0.05, 1.15)
ax2b.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')

# Annotate correlation
ax2b.annotate('max entanglement\n(cat state forms)',
              xy=(0.5 * t_revival * g, 1.0), fontsize=9,
              ha='center', va='top',
              arrowprops=dict(arrowstyle='->', color='red'),
              xytext=(0.5 * t_revival * g + 3, 0.6),
              color='red')

# Field purity
ax2c.plot(tlist * g, data_coh['purity'], color='#1b7837', linewidth=1.0)
ax2c.set_ylabel(r'$\mathrm{Tr}(\rho_{\rm field}^2)$', fontsize=13)
ax2c.set_xlabel(r'$gt$', fontsize=13)
ax2c.set_ylim(-0.05, 1.05)

# Mark revival
for ax in [ax2a, ax2b, ax2c]:
    ax.axvline(t_revival * g, color='purple', linewidth=0.8, linestyle='--', alpha=0.5,
               label=r'$t_{\rm rev}$')

ax2a.legend(fontsize=9)
fig2.tight_layout()
fig2.savefig('/home/claude/fig_coherent_entropy_purity.png', dpi=200, bbox_inches='tight')
fig2.savefig('/home/claude/fig_coherent_entropy_purity.pdf', bbox_inches='tight')
print("Saved: fig_coherent_entropy_purity.png/pdf")

# ============================================================
# Figure 3: Entropy scaling — how max entropy depends on <n>
# ============================================================
n_bar_values = [1, 4, 9, 16, 25, 36]
max_entropies = {}

fig3, axes3 = plt.subplots(2, 3, figsize=(14, 7), sharex=False)
axes3_flat = axes3.flatten()

for idx, nb in enumerate(n_bar_values):
    print(f"Computing entropy for <n> = {nb}...")
    N_cav_local = max(int(3 * nb + 10), 25)
    al = np.sqrt(nb)
    t_rev_local = 2 * np.pi * np.sqrt(nb) / g if nb > 0 else 10.0
    tl = np.linspace(0, 1.5 * t_rev_local, 800)

    a_loc = tensor(destroy(N_cav_local), Qobj(np.eye(2)))
    sm_loc = tensor(Qobj(np.eye(N_cav_local)), sigmam())
    sp_loc = tensor(Qobj(np.eye(N_cav_local)), sigmap())
    H_loc = g * (a_loc.dag() * sm_loc + a_loc * sp_loc)
    rho0_loc = tensor(ket2dm(coherent(N_cav_local, al)), rho_e)

    res_loc = mesolve(H_loc, rho0_loc, tl, [], [], options={'store_states': True})
    S_loc = np.zeros(len(tl))
    for i in range(len(tl)):
        rho_atom_loc = ptrace(res_loc.states[i], 1)
        S_loc[i] = entropy_vn(rho_atom_loc, 2)

    max_entropies[nb] = np.max(S_loc)

    ax = axes3_flat[idx]
    ax.plot(tl * g, S_loc, 'r-', linewidth=0.8)
    ax.set_title(rf'$\bar{{n}} = {nb}$', fontsize=11)
    ax.set_ylim(-0.05, 1.15)
    ax.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
    if idx >= 3:
        ax.set_xlabel(r'$gt$', fontsize=11)
    if idx % 3 == 0:
        ax.set_ylabel(r'$S$ [bits]', fontsize=11)

fig3.suptitle('Atom–Field Entanglement Entropy vs. Mean Photon Number', fontsize=13, y=1.01)
fig3.tight_layout()
fig3.savefig('/home/claude/fig_entropy_nbar_scaling.png', dpi=200, bbox_inches='tight')
fig3.savefig('/home/claude/fig_entropy_nbar_scaling.pdf', bbox_inches='tight')
print("Saved: fig_entropy_nbar_scaling.png/pdf")

# ============================================================
# Figure 4: Dissipative entanglement — entropy vs κ
# ============================================================
kappa_values = [0.0, 0.01, 0.05, 0.1, 0.2]

fig4, (ax4a, ax4b) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
tlist_diss = np.linspace(0, 1.5 * t_revival, 1500)

for kap in kappa_values:
    print(f"Solving dissipative JC with κ/g = {kap:.2f}...")
    c_ops_d = [np.sqrt(kap) * a] if kap > 0 else []
    rho0_d = tensor(ket2dm(coherent(N_cav, alpha)), rho_e)
    res_d = mesolve(H_JC, rho0_d, tlist_diss, c_ops_d, [sz], options={'store_states': True})

    inv_d = np.array(res_d.expect[0])
    S_d = np.zeros(len(tlist_diss))
    for i in range(len(tlist_diss)):
        rho_at = ptrace(res_d.states[i], 1)
        S_d[i] = entropy_vn(rho_at, 2)

    lbl = rf'$\kappa/g = {kap:.2f}$'
    ax4a.plot(tlist_diss * g, inv_d, linewidth=0.7, alpha=0.85, label=lbl)
    ax4b.plot(tlist_diss * g, S_d, linewidth=1.0, alpha=0.9, label=lbl)

ax4a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax4a.set_title(rf'Effect of Cavity Decay on JC Dynamics ($\bar{{n}} = {n_bar:.0f}$)',
               fontsize=13)
ax4a.set_ylim(-1.1, 1.1)
ax4a.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax4a.legend(fontsize=9, ncol=3)

ax4b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax4b.set_xlabel(r'$gt$', fontsize=13)
ax4b.set_ylim(-0.05, 1.15)
ax4b.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
ax4b.legend(fontsize=9, ncol=3)

fig4.tight_layout()
fig4.savefig('/home/claude/fig_dissipative_entanglement.png', dpi=200, bbox_inches='tight')
fig4.savefig('/home/claude/fig_dissipative_entanglement.pdf', bbox_inches='tight')
print("Saved: fig_dissipative_entanglement.png/pdf")

print("\n=== Entanglement dynamics simulations complete ===")
