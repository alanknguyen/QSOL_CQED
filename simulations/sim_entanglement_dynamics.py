# sim_entanglement_dynamics.py
# Von Neumann entropy of the reduced atomic state S(rho_atom)
# as a measure of atom-field entanglement in the Jaynes-Cummings model.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Physics, Boston University
# Advised by Prof. Luca Dal Negro
#
# Generates: fig_entanglement_comparison, fig_coherent_entropy_purity,
#            fig_entropy_nbar_scaling, fig_dissipative_entanglement
#
# S = 0: product state (no entanglement)
# S = 1 bit: maximally entangled atom-field state
# Entropy peaks at collapse, dips at revival (coherent field only).

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap,
    mesolve, ptrace, entropy_vn, coherent, thermal_dm,
    squeeze, Qobj, ket2dm
)

# -- Parameters --
N_cav = 50
g = 1.0
n_bar = 10.0
alpha = np.sqrt(n_bar)
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# -- Joint operators --
a  = tensor(destroy(N_cav), Qobj(np.eye(2)))
sm = tensor(Qobj(np.eye(N_cav)), sigmam())
sp = tensor(Qobj(np.eye(N_cav)), sigmap())
sz = tensor(Qobj(np.eye(N_cav)), sigmaz())

# Resonant JC Hamiltonian (rotating frame, Delta = 0)
H_JC = g * (a.dag() * sm + a * sp)

# -- Atom initial state --
psi_e = basis(2, 0)   # |excited>
rho_e = ket2dm(psi_e)

# -- Field initial states, all with <n> = 10 --
rho_coh  = ket2dm(coherent(N_cav, alpha))          # coherent |alpha>
rho_th   = thermal_dm(N_cav, n_bar)                 # thermal (Bose-Einstein)
r_sq     = np.arcsinh(np.sqrt(n_bar))               # <n> = sinh^2(r)
rho_sq   = ket2dm(squeeze(N_cav, r_sq) * basis(N_cav, 0))  # squeezed vacuum
rho_fock = ket2dm(basis(N_cav, int(n_bar)))          # Fock |10>

initial_states = {
    'Coherent': rho_coh,
    'Thermal': rho_th,
    'Squeezed vacuum': rho_sq,
    'Fock $|10\\rangle$': rho_fock,
}

tlist = np.linspace(0, 1.8 * t_revival, 2000)


# -- Solve for each initial field state --
results = {}

for label, rho_field in initial_states.items():
    print(f"Solving: {label}")
    rho0 = tensor(rho_field, rho_e)
    result = mesolve(H_JC, rho0, tlist, [], [sz], options={'store_states': True})

    S_atom = np.zeros(len(tlist))
    purity_field = np.zeros(len(tlist))
    inversion = np.array(result.expect[0])

    for i in range(len(tlist)):
        rho_atom = ptrace(result.states[i], 1)     # trace out field
        S_atom[i] = entropy_vn(rho_atom, 2)        # base-2 (bits)
        rho_f = ptrace(result.states[i], 0)
        purity_field[i] = (rho_f * rho_f).tr().real

    results[label] = {
        'inversion': inversion,
        'entropy': S_atom,
        'purity': purity_field,
    }
    print(f"  max S = {np.max(S_atom):.4f} bits, min purity = {np.min(purity_field):.4f}")


# Fig 5: Inversion and entropy for all four field states
colors = {'Coherent': '#2166ac', 'Thermal': '#b2182b',
          'Squeezed vacuum': '#1b7837', 'Fock $|10\\rangle$': '#e08214'}

fig1, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

ax1 = axes[0]
for label, data in results.items():
    ax1.plot(tlist * g, data['inversion'], color=colors[label],
             linewidth=0.7, alpha=0.85, label=label)
ax1.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax1.set_title(rf'JC Dynamics: $\bar{{n}} = {n_bar:.0f}$, $\Delta = 0$', fontsize=13)
ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax1.set_ylim(-1.1, 1.1)
ax1.legend(fontsize=10, loc='upper right', ncol=2)

ax2 = axes[1]
for label, data in results.items():
    ax2.plot(tlist * g, data['entropy'], color=colors[label],
             linewidth=1.0, alpha=0.9, label=label)
ax2.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax2.set_xlabel(r'$gt$', fontsize=13)
ax2.set_ylim(-0.05, 1.15)
ax2.axhline(1.0, color='gray', linewidth=0.5, linestyle=':', label=r'$S_{\max}=1$ bit')
ax2.legend(fontsize=10, loc='upper right', ncol=2)

for ax in axes:
    ax.axvline(t_revival * g, color='purple', linewidth=0.8, linestyle='--', alpha=0.5)

fig1.tight_layout()
fig1.savefig('fig_entanglement_comparison.png', dpi=200, bbox_inches='tight')
fig1.savefig('fig_entanglement_comparison.pdf', bbox_inches='tight')
print("\nSaved: fig_entanglement_comparison")


# Fig 6: Coherent state detail: inversion, entropy, field purity
fig2, (ax2a, ax2b, ax2c) = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
data_coh = results['Coherent']

ax2a.plot(tlist * g, data_coh['inversion'], 'b-', linewidth=0.7)
ax2a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax2a.set_title(rf'Coherent State ($\alpha = {alpha:.1f}$, $\bar{{n}} = {n_bar:.0f}$)', fontsize=13)
ax2a.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax2a.set_ylim(-1.1, 1.1)

ax2b.plot(tlist * g, data_coh['entropy'], 'r-', linewidth=1.0)
ax2b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax2b.set_ylim(-0.05, 1.15)
ax2b.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
ax2b.annotate('max entanglement\n(cat state forms)',
              xy=(0.5 * t_revival * g, 1.0), fontsize=9,
              ha='center', va='top',
              arrowprops=dict(arrowstyle='->', color='red'),
              xytext=(0.5 * t_revival * g + 3, 0.6), color='red')

# Purity Tr(rho_field^2): drops to ~0.5 when field is a two-component cat
ax2c.plot(tlist * g, data_coh['purity'], color='#1b7837', linewidth=1.0)
ax2c.set_ylabel(r'$\mathrm{Tr}(\rho_{\rm field}^2)$', fontsize=13)
ax2c.set_xlabel(r'$gt$', fontsize=13)
ax2c.set_ylim(-0.05, 1.05)

for ax in [ax2a, ax2b, ax2c]:
    ax.axvline(t_revival * g, color='purple', linewidth=0.8, linestyle='--', alpha=0.5,
               label=r'$t_{\rm rev}$')

ax2a.legend(fontsize=9)
fig2.tight_layout()
fig2.savefig('fig_coherent_entropy_purity.png', dpi=200, bbox_inches='tight')
fig2.savefig('fig_coherent_entropy_purity.pdf', bbox_inches='tight')
print("Saved: fig_coherent_entropy_purity")


# Fig 7: Entropy vs time for varying <n>
# t_c ~ 1/g (independent of n_bar), t_r = 2pi sqrt(n_bar)/g
n_bar_values = [1, 4, 9, 16, 25, 36]

fig3, axes3 = plt.subplots(2, 3, figsize=(14, 7))
axes3_flat = axes3.flatten()

for idx, nb in enumerate(n_bar_values):
    print(f"Entropy for <n> = {nb}...")
    N_loc = max(int(3 * nb + 10), 25)
    t_rev_loc = 2 * np.pi * np.sqrt(nb) / g if nb > 0 else 10.0
    tl = np.linspace(0, 1.5 * t_rev_loc, 800)

    a_loc  = tensor(destroy(N_loc), Qobj(np.eye(2)))
    sm_loc = tensor(Qobj(np.eye(N_loc)), sigmam())
    sp_loc = tensor(Qobj(np.eye(N_loc)), sigmap())
    H_loc  = g * (a_loc.dag() * sm_loc + a_loc * sp_loc)
    rho0_loc = tensor(ket2dm(coherent(N_loc, np.sqrt(nb))), rho_e)

    res = mesolve(H_loc, rho0_loc, tl, [], [], options={'store_states': True})
    S_loc = np.zeros(len(tl))
    for i in range(len(tl)):
        S_loc[i] = entropy_vn(ptrace(res.states[i], 1), 2)

    ax = axes3_flat[idx]
    ax.plot(tl * g, S_loc, 'r-', linewidth=0.8)
    ax.set_title(rf'$\bar{{n}} = {nb}$', fontsize=11)
    ax.set_ylim(-0.05, 1.15)
    ax.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
    if idx >= 3:
        ax.set_xlabel(r'$gt$', fontsize=11)
    if idx % 3 == 0:
        ax.set_ylabel(r'$S$ [bits]', fontsize=11)

fig3.suptitle('Entanglement Entropy vs. Mean Photon Number', fontsize=13, y=1.01)
fig3.tight_layout()
fig3.savefig('fig_entropy_nbar_scaling.png', dpi=200, bbox_inches='tight')
fig3.savefig('fig_entropy_nbar_scaling.pdf', bbox_inches='tight')
print("Saved: fig_entropy_nbar_scaling")


# Fig 8: Dissipative dynamics, Lindblad with cavity loss kappa
# Shows damping of revivals and entropy saturation
kappa_values = [0.0, 0.01, 0.05, 0.1, 0.2]

fig4, (ax4a, ax4b) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
tlist_diss = np.linspace(0, 1.5 * t_revival, 1500)

for kap in kappa_values:
    print(f"Dissipative JC: kappa/g = {kap:.2f}")
    c_ops = [np.sqrt(kap) * a] if kap > 0 else []
    rho0 = tensor(ket2dm(coherent(N_cav, alpha)), rho_e)
    res = mesolve(H_JC, rho0, tlist_diss, c_ops, [sz], options={'store_states': True})

    inv = np.array(res.expect[0])
    S = np.zeros(len(tlist_diss))
    for i in range(len(tlist_diss)):
        S[i] = entropy_vn(ptrace(res.states[i], 1), 2)

    lbl = rf'$\kappa/g = {kap:.2f}$'
    ax4a.plot(tlist_diss * g, inv, linewidth=0.7, alpha=0.85, label=lbl)
    ax4b.plot(tlist_diss * g, S, linewidth=1.0, alpha=0.9, label=lbl)

ax4a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax4a.set_title(rf'Cavity Decay Effect ($\bar{{n}} = {n_bar:.0f}$)', fontsize=13)
ax4a.set_ylim(-1.1, 1.1)
ax4a.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax4a.legend(fontsize=9, ncol=3)

ax4b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax4b.set_xlabel(r'$gt$', fontsize=13)
ax4b.set_ylim(-0.05, 1.15)
ax4b.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
ax4b.legend(fontsize=9, ncol=3)

fig4.tight_layout()
fig4.savefig('fig_dissipative_entanglement.png', dpi=200, bbox_inches='tight')
fig4.savefig('fig_dissipative_entanglement.pdf', bbox_inches='tight')
print("Saved: fig_dissipative_entanglement")