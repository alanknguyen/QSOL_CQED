# sim_entanglement_dynamics.py
# Von Neumann entropy of the reduced atomic state S(rho_atom) and field
# purity Tr[rho_field^2] in the Jaynes-Cummings model.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Physics, Boston University
# Advised by Prof. Luca Dal Negro
#
# Generates: fig_entanglement_comparison, fig_coherent_entropy_purity,
#            fig_entropy_nbar_scaling, fig_dissipative_entanglement
#
# Interpretation (coherent field, atom initially excited):
#   * early collapse (gt ~ 1-3): S -> 1 bit, atom and field strongly
#     entangled; the reduced field is a two-branch MIXTURE (P ~ 0.5).
#   * half revival t_r/2: S reaches a MINIMUM (~0.15 bit) and the field
#     purity a maximum (~0.96): the atom and field approximately
#     disentangle and the field is a nearly pure Schrodinger cat
#     (Gea-Banacloche, PRL 65, 3385 (1990); Phoenix & Knight 1991).
#   * full revival t_r: the branches rephase, S is again close to 1 bit.
#
# S(rho_atom) is an entanglement measure only for a globally PURE
# atom-field state (coherent, Fock, squeezed-vacuum initial fields with
# no dissipation).  For the thermal field and for kappa > 0 the global
# state is mixed and S(rho_atom) measures total (quantum + classical)
# correlations.  We therefore also compute the logarithmic negativity
#   E_N = log2 || rho^{T_atom} ||_1,
# which is a proper entanglement monotone for mixed states (E_N = 1 for a
# maximally entangled atom-field pair, 0 for separable states).

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap, mesolve, ptrace,
    entropy_vn, coherent, thermal_dm, squeeze, qeye, ket2dm, expect, num,
    negativity,
)

PAPER_FIG = Path(__file__).resolve().parents[1] / 'paper' / 'figures'
PAPER_FIG.mkdir(parents=True, exist_ok=True)

# -- Parameters --
g = 1.0
n_bar = 10.0
alpha = np.sqrt(n_bar)
t_revival = 2 * np.pi * np.sqrt(n_bar) / g
t_collapse = 1.0 / g

# Fock-space truncation per initial state.  Coherent and Fock states with
# n_bar = 10 have negligible weight above n = 50.  Thermal and squeezed-vacuum
# photon distributions have geometric / heavy tails: at N = 50 the truncated
# states have <n> = 9.57 (thermal) and 9.15 (squeezed) instead of 10.
# N = 150 restores <n> = 10.000 and keeps the truncated weight below 1e-6.
N_CAV = {
    'Coherent': 50,
    'Thermal': 150,
    'Squeezed vacuum': 150,
    'Fock $|10\\rangle$': 50,
}
COLORS = {'Coherent': '#2166ac', 'Thermal': '#b2182b',
          'Squeezed vacuum': '#1b7837', 'Fock $|10\\rangle$': '#e08214'}


def jc_ops(N):
    a = tensor(destroy(N), qeye(2))
    sm = tensor(qeye(N), sigmam())
    sp = tensor(qeye(N), sigmap())
    sz = tensor(qeye(N), sigmaz())
    H = g * (a.dag() * sm + a * sp)        # resonant JC, rotating frame
    return a, sz, H


def field_state(label, N):
    if label == 'Coherent':
        return coherent(N, alpha)
    if label == 'Thermal':
        return thermal_dm(N, n_bar)
    if label == 'Squeezed vacuum':
        r = np.arcsinh(np.sqrt(n_bar))       # <n> = sinh^2(r)
        return squeeze(N, r) * basis(N, 0)
    if label.startswith('Fock'):
        return basis(N, int(n_bar))
    raise ValueError(label)


def run_jc(label, N, tlist, kappa=0.0):
    """Evolve |e> x rho_field; return inversion, S(rho_atom), Tr[rho_field^2]."""
    a, sz, H = jc_ops(N)
    phi = field_state(label, N)
    atom = basis(2, 0)
    if phi.isket:
        rho0 = tensor(phi, atom)              # ket: sesolve path, much faster
    else:
        rho0 = tensor(phi, ket2dm(atom))
    n0 = expect(num(N), phi)
    c_ops = [np.sqrt(kappa) * a] if kappa > 0 else []
    res = mesolve(H, rho0, tlist, c_ops, [sz], options={'store_states': True})
    S = np.array([entropy_vn(ptrace(st, 1), 2) for st in res.states])
    P = np.array([(ptrace(st, 0) ** 2).tr().real for st in res.states])
    EN = np.array([negativity(st if st.isoper else ket2dm(st), 1, logarithmic=True)
                   for st in res.states])
    return dict(inversion=np.array(res.expect[0]), entropy=S, purity=P,
                logneg=EN, n0=n0)


# ---------------------------------------------------------------------------
# Solve for the four initial field states
# ---------------------------------------------------------------------------
tlist = np.linspace(0, 1.8 * t_revival, 2000)
results = {}
for label, N in N_CAV.items():
    print(f"Solving: {label}  (N_cav = {N})")
    results[label] = run_jc(label, N, tlist)
    d = results[label]
    print(f"  <n>(0) = {d['n0']:.4f}   max S = {d['entropy'].max():.3f} bits   "
          f"max E_N = {d['logneg'].max():.3f}   min P = {d['purity'].min():.3f}")

i_half = np.argmin(np.abs(tlist - 0.5 * t_revival))
i_rev = np.argmin(np.abs(tlist - t_revival))
coh = results['Coherent']
print(f"Coherent: S(t_r/2) = {coh['entropy'][i_half]:.3f}, S(t_r) = {coh['entropy'][i_rev]:.3f}, "
      f"P(t_r/2) = {coh['purity'][i_half]:.3f}, P(t_r) = {coh['purity'][i_rev]:.3f}, "
      f"E_N(t_r/2) = {coh['logneg'][i_half]:.3f}, E_N(t_r) = {coh['logneg'][i_rev]:.3f}")
th = results['Thermal']
print(f"Thermal: S(gt>2) mean = {th['entropy'][tlist > 2].mean():.3f}, "
      f"E_N mean (gt>2) = {th['logneg'][tlist > 2].mean():.3f}, E_N max = {th['logneg'].max():.3f}")


def mark_times(ax, labels=False):
    """Vertical guides at t_r/2 (dotted) and t_r (dashed); labels go in the legend."""
    ax.axvline(0.5 * t_revival * g, color='crimson', lw=0.8, ls=':', alpha=0.7,
               label=r'$t_r/2$' if labels else None)
    ax.axvline(t_revival * g, color='purple', lw=0.8, ls='--', alpha=0.5,
               label=r'$t_r$' if labels else None)


# ---------------------------------------------------------------------------
# Fig 5: inversion and entropy for the four field states
# ---------------------------------------------------------------------------
fock = 'Fock $|10\\rangle$'
fig1 = plt.figure(figsize=(11, 9.5))
gs = fig1.add_gridspec(3, 1, height_ratios=[3, 1.3, 3], hspace=0.32)
ax1 = fig1.add_subplot(gs[0])
axf = fig1.add_subplot(gs[1])
ax2 = fig1.add_subplot(gs[2], sharex=ax1)

# (a) inversion: the three states with a spread of Rabi frequencies
for label in ('Coherent', 'Thermal', 'Squeezed vacuum'):
    ax1.plot(tlist * g, results[label]['inversion'], color=COLORS[label],
             lw=0.7, alpha=0.9, label=label)
ax1.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax1.set_title(rf'(a) Atomic inversion, $\bar{{n}} = {n_bar:.0f}$, $\Delta = 0$', fontsize=12)
ax1.axhline(0, color='gray', lw=0.5, ls='--')
ax1.set_ylim(-1.1, 1.1)
mark_times(ax1, labels=True)
ax1.legend(fontsize=8, loc='upper right', ncol=5, framealpha=0.95)
ax1.tick_params(labelbottom=False)

# (b) Fock |10>: a single Rabi frequency 2g sqrt(11), no collapse; short window
mwin = tlist * g < 3.0
axf.plot(tlist[mwin] * g, results[fock]['inversion'][mwin], color=COLORS[fock], lw=1.0,
         label=r'$\langle\sigma_z\rangle$, period $\pi/(g\sqrt{11})$')
axf.plot(tlist[mwin] * g, 2 * results[fock]['entropy'][mwin] - 1, color='#7f3b08', lw=0.9, ls='--',
         label=r'$2S - 1$, period $\pi/(2g\sqrt{11})$')
axf.axhline(0, color='gray', lw=0.4, ls='--')
axf.set_ylim(-1.15, 1.15)
axf.set_xlim(0, 3.0)
axf.set_xlabel(r'$gt$', fontsize=10)
axf.set_title(r'(b) Fock $|10\rangle$: periodic Rabi oscillation and entropy (no collapse)', fontsize=10)
axf.legend(fontsize=7, loc='upper right', ncol=2, framealpha=0.95)
axf.tick_params(labelsize=8)

# (c) entropy (solid) and logarithmic negativity (dashed)
for label in ('Coherent', 'Thermal', 'Squeezed vacuum'):
    ax2.plot(tlist * g, results[label]['entropy'], color=COLORS[label],
             lw=1.0, alpha=0.9, label=label)
    ax2.plot(tlist * g, results[label]['logneg'], color=COLORS[label],
             lw=0.9, ls='--', alpha=0.8)
ax2.plot([], [], color='k', lw=1.0, label=r'$S(\rho_{\rm atom})$')
ax2.plot([], [], color='k', lw=0.9, ls='--', label=r'$E_{\mathcal{N}}$')
ax2.set_ylabel(r'$S(\rho_{\rm atom})$ [bits],  $E_{\mathcal{N}}$', fontsize=12)
ax2.set_xlabel(r'$gt$', fontsize=13)
ax2.set_title(r'(c) Reduced-atom entropy $S$ (solid) and logarithmic negativity $E_{\mathcal{N}}$ (dashed)',
              fontsize=11)
ax2.set_ylim(-0.05, 1.15)
ax2.axhline(1.0, color='gray', lw=0.5, ls=':')
mark_times(ax2, labels=True)
ax2.legend(fontsize=7, loc='lower right', ncol=4, framealpha=0.95)
ax2.text(0.01, 0.03, r'thermal: global state mixed, $S \approx 1$ but $E_{\mathcal{N}} \approx 0.2$',
         transform=ax2.transAxes, fontsize=7, color=COLORS['Thermal'])

fig1.savefig(PAPER_FIG / 'fig_entanglement_comparison.png', dpi=200, bbox_inches='tight')
fig1.savefig(PAPER_FIG / 'fig_entanglement_comparison.pdf', bbox_inches='tight')
print("\nSaved: fig_entanglement_comparison")


# ---------------------------------------------------------------------------
# Fig 6: coherent state detail: inversion, entropy, field purity
# ---------------------------------------------------------------------------
fig2, (ax2a, ax2b, ax2c) = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
S = coh['entropy']; P = coh['purity']

ax2a.plot(tlist * g, coh['inversion'], 'b-', lw=0.7)
ax2a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax2a.set_title(rf'Coherent-state field ($\alpha = {alpha:.2f}$, $\bar{{n}} = {n_bar:.0f}$)', fontsize=13)
ax2a.axhline(0, color='gray', lw=0.5, ls='--')
ax2a.set_ylim(-1.1, 1.1)
mark_times(ax2a, labels=True)
ax2a.legend(fontsize=9, loc='upper right', ncol=2, framealpha=0.95)

ax2b.plot(tlist * g, S, 'r-', lw=1.0)
ax2b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax2b.set_ylim(-0.05, 1.15)
ax2b.axhline(1.0, color='gray', lw=0.5, ls=':')
i_peak = np.argmax(S[tlist * g < 5])
ax2b.annotate('near-maximal entanglement\nduring collapse',
              xy=(tlist[i_peak] * g, S[i_peak]), xytext=(6.5, 0.95),
              fontsize=9, color='red', ha='left', va='center',
              arrowprops=dict(arrowstyle='->', color='red'))
ax2b.annotate(rf'entropy minimum $S \approx {S[i_half]:.2f}$ at $t_r/2$:' '\n'
              'atom and field disentangle,\nfield $\\approx$ pure cat state',
              xy=(tlist[i_half] * g, S[i_half]), xytext=(13.5, 0.30),
              fontsize=9, color='red', ha='left', va='center',
              arrowprops=dict(arrowstyle='->', color='red'))
mark_times(ax2b)

ax2c.plot(tlist * g, P, color='#1b7837', lw=1.0)
ax2c.set_ylabel(r'$\mathcal{P} = \mathrm{Tr}(\rho_{\rm field}^2)$', fontsize=13)
ax2c.set_xlabel(r'$gt$', fontsize=13)
ax2c.set_ylim(-0.05, 1.05)
ax2c.axhline(0.5, color='gray', lw=0.5, ls=':')
ax2c.annotate(rf'$\mathcal{{P}} \approx {P[i_half]:.2f}$: nearly pure cat',
              xy=(tlist[i_half] * g, P[i_half]), xytext=(13.5, 0.92),
              fontsize=9, color='#1b7837', ha='left',
              arrowprops=dict(arrowstyle='->', color='#1b7837'))
i_pmin = np.argmin(P[tlist * g < 5])
ax2c.annotate(r'$\mathcal{P} \to 0.5$: two-branch mixture' '\n(field entangled with atom)',
              xy=(tlist[i_pmin] * g, P[i_pmin]), xytext=(6.5, 0.25),
              fontsize=9, color='#1b7837', ha='left', va='center',
              arrowprops=dict(arrowstyle='->', color='#1b7837'))
mark_times(ax2c)

fig2.tight_layout()
fig2.savefig(PAPER_FIG / 'fig_coherent_entropy_purity.png', dpi=200, bbox_inches='tight')
fig2.savefig(PAPER_FIG / 'fig_coherent_entropy_purity.pdf', bbox_inches='tight')
print("Saved: fig_coherent_entropy_purity")


# ---------------------------------------------------------------------------
# Fig 7: entropy vs time for varying <n>  (t_c ~ 1/g, t_r = 2 pi sqrt(n)/g)
# ---------------------------------------------------------------------------
n_bar_values = [1, 4, 9, 16, 25, 36]
fig3, axes3 = plt.subplots(2, 3, figsize=(14, 7))

for idx, nb in enumerate(n_bar_values):
    print(f"Entropy for <n> = {nb}...")
    N_loc = int(3 * nb + 20)
    t_rev_loc = 2 * np.pi * np.sqrt(nb) / g
    tl = np.linspace(0, 1.5 * t_rev_loc, 800)
    a_loc, sz_loc, H_loc = jc_ops(N_loc)
    psi0_loc = tensor(coherent(N_loc, np.sqrt(nb)), basis(2, 0))
    res = mesolve(H_loc, psi0_loc, tl, [], [], options={'store_states': True})
    S_loc = np.array([entropy_vn(ptrace(st, 1), 2) for st in res.states])

    ax = axes3.flat[idx]
    ax.plot(tl * g, S_loc, 'r-', lw=0.8)
    ax.axvline(0.5 * t_rev_loc * g, color='crimson', lw=0.7, ls=':', alpha=0.7)
    ax.axvline(t_rev_loc * g, color='purple', lw=0.7, ls='--', alpha=0.5)
    ax.set_title(rf'$\bar{{n}} = {nb}$', fontsize=11)
    ax.set_ylim(-0.05, 1.15)
    ax.axhline(1.0, color='gray', lw=0.5, ls=':')
    if idx >= 3:
        ax.set_xlabel(r'$gt$', fontsize=11)
    if idx % 3 == 0:
        ax.set_ylabel(r'$S$ [bits]', fontsize=11)

fig3.suptitle('Entanglement entropy vs. mean photon number '
              r'(dotted: $t_r/2$, dashed: $t_r$)', fontsize=13, y=1.01)
fig3.tight_layout()
fig3.savefig(PAPER_FIG / 'fig_entropy_nbar_scaling.png', dpi=200, bbox_inches='tight')
fig3.savefig(PAPER_FIG / 'fig_entropy_nbar_scaling.pdf', bbox_inches='tight')
print("Saved: fig_entropy_nbar_scaling")


# ---------------------------------------------------------------------------
# Fig 8: dissipative dynamics, Lindblad with cavity loss kappa
# ---------------------------------------------------------------------------
kappa_values = [0.0, 0.01, 0.05, 0.1, 0.2]
tlist_diss = np.linspace(0, 1.5 * t_revival, 1500)
fig4, (ax4a, ax4b, ax4c) = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
ih = np.argmin(np.abs(tlist_diss - 0.5 * t_revival))

for kap in kappa_values:
    print(f"Dissipative JC: kappa/g = {kap:.2f}")
    d = run_jc('Coherent', N_CAV['Coherent'], tlist_diss, kappa=kap)
    lbl = rf'$\kappa/g = {kap:.2f}$'
    ax4a.plot(tlist_diss * g, d['inversion'], lw=0.7, alpha=0.85, label=lbl)
    ax4b.plot(tlist_diss * g, d['entropy'], lw=1.0, alpha=0.9, label=lbl)
    ax4c.plot(tlist_diss * g, d['logneg'], lw=1.0, alpha=0.9, label=lbl)
    print(f"   S(t_r/2) = {d['entropy'][ih]:.3f}, E_N(t_r/2) = {d['logneg'][ih]:.3f}, "
          f"E_N(t_end) = {d['logneg'][-1]:.3f}, S(t_end) = {d['entropy'][-1]:.3f}")

ax4a.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=13)
ax4a.set_title(rf'Effect of cavity decay ($\bar{{n}} = {n_bar:.0f}$)', fontsize=13)
ax4a.set_ylim(-1.1, 1.1)
ax4a.axhline(0, color='gray', lw=0.5, ls='--')
mark_times(ax4a, labels=True)
ax4a.legend(fontsize=9, ncol=4, loc='upper right', framealpha=0.95)

ax4b.set_ylabel(r'$S(\rho_{\rm atom})$ [bits]', fontsize=13)
ax4b.set_ylim(-0.05, 1.15)
ax4b.axhline(1.0, color='gray', lw=0.5, ls=':')
ax4b.legend(fontsize=9, ncol=3, loc='lower right', framealpha=0.95)
ax4b.text(0.01, 0.03, r'for $\kappa > 0$ the global state is mixed: $S(\rho_{\rm atom})$ = total correlations',
          transform=ax4b.transAxes, fontsize=7, color='gray')
mark_times(ax4b)

ax4c.set_ylabel(r'$E_{\mathcal{N}}$ (log-negativity)', fontsize=12)
ax4c.set_xlabel(r'$gt$', fontsize=13)
ax4c.set_ylim(-0.05, 1.15)
ax4c.axhline(1.0, color='gray', lw=0.5, ls=':')
ax4c.legend(fontsize=9, ncol=3, loc='upper right', framealpha=0.95)
ax4c.text(0.01, 0.03, r'$E_{\mathcal{N}} = \log_2\|\rho^{T_{\rm atom}}\|_1$: entanglement monotone valid for mixed states',
          transform=ax4c.transAxes, fontsize=7, color='gray')
mark_times(ax4c)

fig4.tight_layout()
fig4.savefig(PAPER_FIG / 'fig_dissipative_entanglement.png', dpi=200, bbox_inches='tight')
fig4.savefig(PAPER_FIG / 'fig_dissipative_entanglement.pdf', bbox_inches='tight')
print("Saved: fig_dissipative_entanglement")
