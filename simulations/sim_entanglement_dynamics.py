# sim_entanglement_dynamics.py
# Atom-field entanglement in the Jaynes-Cummings model: von Neumann entropy of
# the reduced atom S(rho_atom), logarithmic negativity E_N, and field purity
# Tr[rho_field^2].
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Physics, Boston University
# Began as a course project taught by Prof. Luca Dal Negro
#
# Generates (paper/figures/): fig_entanglement_comparison,
#     fig_coherent_entropy_purity, fig_entropy_nbar_scaling,
#     fig_dissipative_entanglement
#
# Interpretation (coherent field, atom initially excited):
#   * collapse (gt ~ 1-3): S -> 1 bit, atom and field strongly entangled; the
#     reduced field is a two-branch MIXTURE (purity ~ 0.5).
#   * half revival t_r/2: S reaches a MINIMUM (0.14 bit) and the field purity
#     a maximum (0.96): the atom and field approximately disentangle and the
#     field is a nearly pure Schrodinger cat
#     (Gea-Banacloche, PRL 65, 3385 (1990); Phoenix & Knight 1991).
#   * revival t_r: the branches meet again and S is back near 0.8 bit.
#
# S(rho_atom) is an entanglement measure only for a globally PURE atom-field
# state (coherent, Fock, squeezed-vacuum initial fields without dissipation).
# For the thermal field and for kappa > 0 the global state is mixed and
# S(rho_atom) also counts classical correlations.  We therefore also compute
#   E_N = log2 || rho^{T_atom} ||_1,
# an entanglement monotone for mixed states (E_N = 1 for a maximally entangled
# atom-field pair, 0 for separable states).  For a pure state with Schmidt
# probabilities (l1, l2) it reduces to E_N = log2(1 + 2 sqrt(l1 l2)).

import numpy as np
import matplotlib.pyplot as plt
from qutip import (basis, tensor, destroy, sigmaz, sigmam, mesolve, entropy_vn, coherent,
                   thermal_dm, qeye, ket2dm, expect, num, negativity)

from plotstyle import (use_style, save, panel_label, time_guides, sequential_colors,
                       FULL_W, COLORS, GUIDE_HALF, GUIDE_FULL, PAPER_FIG)
from utils import squeezed_vacuum

use_style()

# -- Parameters --
g = 1.0
n_bar = 10.0
alpha = np.sqrt(n_bar)
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# Fock-space truncation per initial state.  Coherent and Fock states with
# n_bar = 10 have negligible weight above n = 50.  Thermal and squeezed-vacuum
# photon distributions have heavy tails: at N = 50 the truncated states have
# <n> = 9.57 (thermal) and 9.15 (squeezed) instead of 10.  The weight beyond
# N is (10/11)^N for the thermal state (6e-7 at N = 150) and falls more slowly
# for the squeezed vacuum (1.6e-4 at N = 150, 9e-8 at N = 300).
STATES = {
    #  key        label              N_cav  colour
    'coherent': ('coherent',          50, COLORS['blue']),
    'thermal':  ('thermal',          150, COLORS['red']),
    'squeezed': ('squeezed vacuum',  300, COLORS['green']),
    'fock':     (r'Fock $|10\rangle$', 50, COLORS['orange']),
}


def jc_ops(N):
    a = tensor(destroy(N), qeye(2))
    sm = tensor(qeye(N), sigmam())
    sz = tensor(qeye(N), sigmaz())
    H = g * (a.dag() * sm + a * sm.dag())       # resonant JC, interaction picture
    return a, sz, H


def field_state(key, N):
    if key == 'coherent':
        return coherent(N, alpha)
    if key == 'thermal':
        return thermal_dm(N, n_bar)
    if key == 'squeezed':
        return squeezed_vacuum(N, n_bar)        # <n> = sinh^2(r) = n_bar
    if key == 'fock':
        return basis(N, int(n_bar))
    raise ValueError(key)


def log_negativity(state):
    """E_N = log2 ||rho^{T_atom}||_1; closed form for kets, trace norm otherwise."""
    if state.isket:
        lam = np.clip(np.linalg.eigvalsh(state.ptrace(1).full()), 0.0, 1.0)
        return float(np.log2(1.0 + 2.0 * np.sqrt(lam[0] * lam[1])))
    return float(negativity(state, 1, logarithmic=True))


def run_jc(key, N, tlist, kappa=0.0):
    """Evolve |e> x rho_field and return inversion, S(rho_atom), E_N, field purity."""
    a, sz, H = jc_ops(N)
    phi = field_state(key, N)
    atom = basis(2, 0)
    # kets are propagated as kets (fast); the thermal field needs the density matrix
    rho0 = tensor(phi, atom) if phi.isket else tensor(phi, ket2dm(atom))
    c_ops = [np.sqrt(kappa) * a] if kappa > 0 else []
    res = mesolve(H, rho0, tlist, c_ops=c_ops, e_ops=[sz], options={'store_states': True})
    return dict(
        inversion=np.array(res.expect[0]),
        entropy=np.array([entropy_vn(st.ptrace(1), 2) for st in res.states]),
        purity=np.array([(st.ptrace(0) ** 2).tr().real for st in res.states]),
        logneg=np.array([log_negativity(st) for st in res.states]),
        n0=expect(num(N), phi),
    )


# ---------------------------------------------------------------------------
# Solve for the four initial field states
# ---------------------------------------------------------------------------
tlist = np.linspace(0, 1.8 * t_revival, 2400)
t_fock = np.linspace(0, 3.0 / g, 600)
results = {}
for key, (label, N, _) in STATES.items():
    print(f"Solving: {key}  (N_cav = {N})")
    results[key] = run_jc(key, N, t_fock if key == 'fock' else tlist)
    d = results[key]
    print(f"  <n>(0) = {d['n0']:.5f}   max S = {d['entropy'].max():.3f} bit   "
          f"max E_N = {d['logneg'].max():.3f}   min purity = {d['purity'].min():.3f}")

# closed-form E_N for kets agrees with the general trace-norm routine
_a, _sz, _H = jc_ops(50)
_st = mesolve(_H, tensor(coherent(50, alpha), basis(2, 0)), [0, 3.3]).states[-1]
assert abs(log_negativity(_st) - negativity(ket2dm(_st), 1, logarithmic=True)) < 1e-7

i_half = int(np.argmin(np.abs(tlist - 0.5 * t_revival)))
i_rev = int(np.argmin(np.abs(tlist - t_revival)))
coh, th, sq = results['coherent'], results['thermal'], results['squeezed']
print(f"Coherent: S(t_r/2) = {coh['entropy'][i_half]:.3f}, S(t_r) = {coh['entropy'][i_rev]:.3f}, "
      f"P(t_r/2) = {coh['purity'][i_half]:.3f}, P(t_r) = {coh['purity'][i_rev]:.3f}, "
      f"E_N(t_r/2) = {coh['logneg'][i_half]:.3f}, E_N(t_r) = {coh['logneg'][i_rev]:.3f}")
late = tlist * g > 2
print(f"Thermal: S(gt>2) mean = {th['entropy'][late].mean():.3f}, min = {th['entropy'][late].min():.3f}; "
      f"E_N(gt>2) mean = {th['logneg'][late].mean():.3f}, max = {th['logneg'].max():.3f}")
print(f"Squeezed: S(gt>2) min = {sq['entropy'][late].min():.3f}, mean = {sq['entropy'][late].mean():.3f}")


# ---------------------------------------------------------------------------
# Fig 5: inversion and entanglement for the four field states
#        rows = initial field, left = <sigma_z>, right = S (coloured) and E_N (thin black)
# ---------------------------------------------------------------------------
fig5 = plt.figure(figsize=(FULL_W, 6.1))
gs5 = fig5.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1.2])
ax_t = np.array([[fig5.add_subplot(gs5[r, c]) for c in range(2)] for r in range(3)])
ax_b = [fig5.add_subplot(gs5[3, c]) for c in range(2)]
for r in range(2):                       # rows (a)-(f) share the time axis
    for c in range(2):
        ax_t[r, c].sharex(ax_t[2, c])
        ax_t[r, c].tick_params(labelbottom=False)
letters = iter('abcdefgh')


def draw_row(ax_inv, ax_ent, key, t):
    label, _, color = STATES[key]
    d = results[key]
    ax_inv.plot(t * g, d['inversion'], color=color, lw=0.6)
    ax_inv.axhline(0, color='0.6', lw=0.4)
    ax_inv.set_ylim(-1.1, 1.1)
    ax_inv.set_yticks([-1, 0, 1])
    ax_inv.set_ylabel(r'$\langle \sigma_z \rangle$')
    panel_label(ax_inv, next(letters), f'{label} field')
    ax_ent.plot(t * g, d['entropy'], color=color, lw=1.0, label=r'$S(\rho_{\rm atom})$')
    # thin solid black: on the rapidly fluctuating thermal trace a dashed line reads as noise
    ax_ent.plot(t * g, d['logneg'], color='k', lw=0.55, label=r'$E_{\mathcal{N}}$')
    ax_ent.axhline(1.0, color='0.6', lw=0.4)
    ax_ent.set_ylim(-0.05, 1.1)
    ax_ent.set_yticks([0, 0.5, 1])
    ax_ent.set_ylabel(r'$S$ (bit), $E_{\mathcal{N}}$')
    panel_label(ax_ent, next(letters))


for row, key in enumerate(('coherent', 'thermal', 'squeezed')):
    draw_row(ax_t[row, 0], ax_t[row, 1], key, tlist)
    for ax in ax_t[row]:
        time_guides(ax, t_revival * g)
        ax.set_xlim(0, tlist[-1] * g)
ax_t[0, 0].text(0.5 * t_revival * g + 0.4, -1.02, r'$t_r/2$', fontsize=8, va='bottom',
                color=GUIDE_HALF['color'])
ax_t[0, 0].text(t_revival * g + 0.4, -1.02, r'$t_r$', fontsize=8, va='bottom',
                color=GUIDE_FULL['color'])
ax_t[0, 1].legend(loc='lower right', bbox_to_anchor=(1.0, 0.98), ncol=2, columnspacing=1.0,
                  borderaxespad=0.0)
for ax in ax_t[-1]:
    ax.set_xlabel(r'$gt$')

draw_row(ax_b[0], ax_b[1], 'fock', t_fock)
for ax in ax_b:
    ax.set_xlim(0, t_fock[-1] * g)
    ax.set_xlabel(r'$gt$')
save(fig5, PAPER_FIG / 'fig_entanglement_comparison')


# ---------------------------------------------------------------------------
# Fig 6: coherent field in detail: inversion, entropy, field purity
# ---------------------------------------------------------------------------
S, P = coh['entropy'], coh['purity']
fig6, (ax6a, ax6b, ax6c) = plt.subplots(3, 1, figsize=(FULL_W, 4.9), sharex=True)
arrow = dict(arrowstyle='->', lw=0.6, shrinkA=1, shrinkB=2)
# white backing, so the t_r/2 and t_r guide lines do not strike through the annotations
note_box = dict(boxstyle='square,pad=0.15', fc='white', ec='none')

ax6a.plot(tlist * g, coh['inversion'], color=COLORS['blue'], lw=0.7)
ax6a.axhline(0, color='0.6', lw=0.4)
ax6a.set_ylabel(r'$\langle \sigma_z \rangle$')
ax6a.set_ylim(-1.1, 1.1)
ax6a.set_yticks([-1, 0, 1])
time_guides(ax6a, t_revival * g, labels=True)
panel_label(ax6a, 'a')

ax6b.plot(tlist * g, S, color=COLORS['red'], lw=1.0)
ax6b.axhline(1.0, color='0.6', lw=0.4)
ax6b.set_ylabel(r'$S(\rho_{\rm atom})$ (bit)')
ax6b.set_ylim(-0.05, 1.18)
ax6b.set_yticks([0, 0.5, 1])
time_guides(ax6b, t_revival * g)
i_peak = int(np.argmax(S[tlist * g < 5]))
ax6b.annotate('collapse: atom and field entangled',
              xy=(tlist[i_peak] * g, S[i_peak]), xytext=(3.6, 1.09),
              fontsize=8, ha='left', va='center', arrowprops=arrow, bbox=note_box)
ax6b.annotate(rf'$S = {S[i_half]:.2f}$ bit: disentangled',
              xy=(tlist[i_half] * g, S[i_half]), xytext=(12.6, 0.07),
              fontsize=8, ha='left', va='center', arrowprops=arrow, bbox=note_box)
panel_label(ax6b, 'b')

ax6c.plot(tlist * g, P, color=COLORS['green'], lw=1.0)
ax6c.axhline(0.5, color='0.6', lw=0.4)
ax6c.set_ylabel(r'$\mathcal{P} = \mathrm{Tr}\,\rho_{\rm field}^2$')
ax6c.set_xlabel(r'$gt$')
ax6c.set_ylim(0.38, 1.06)
ax6c.set_yticks([0.5, 0.75, 1.0])
time_guides(ax6c, t_revival * g)
i_pmin = int(np.argmin(P[tlist * g < 5]))
ax6c.annotate(rf'$\mathcal{{P}} = {P[i_half]:.2f}$: nearly pure cat state',
              xy=(tlist[i_half] * g, P[i_half]), xytext=(12.2, 1.005),
              fontsize=8, ha='left', va='center', arrowprops=arrow, bbox=note_box)
ax6c.annotate(r'$\mathcal{P} \approx 1/2$: two-branch mixture',
              xy=(tlist[i_pmin] * g, P[i_pmin]), xytext=(4.6, 0.44),
              fontsize=8, ha='left', va='center', arrowprops=arrow, bbox=note_box)
panel_label(ax6c, 'c')
ax6c.set_xlim(0, tlist[-1] * g)
save(fig6, PAPER_FIG / 'fig_coherent_entropy_purity')


# ---------------------------------------------------------------------------
# Fig 7: entropy for increasing <n>  (collapse ~ 1/g, t_r = 2 pi sqrt(n)/g)
# ---------------------------------------------------------------------------
n_bar_values = [1, 4, 9, 16, 25, 36]
fig7, axes7 = plt.subplots(2, 3, figsize=(FULL_W, 3.5), sharey=True)
print("Entropy minimum near t_r/2 for a coherent field:")
for ax, letter, nb in zip(axes7.flat, 'abcdef', n_bar_values):
    N_loc = int(3 * nb + 20)
    t_rev_loc = 2 * np.pi * np.sqrt(nb) / g
    # sample the entropy oscillation (period pi / (2 g sqrt(n+1))) at >= 25 points
    tl = np.linspace(0, 1.5 * t_rev_loc, int(1.5 * t_rev_loc / 0.01) + 1)
    _, _, H_loc = jc_ops(N_loc)
    psi0_loc = tensor(coherent(N_loc, np.sqrt(nb)), basis(2, 0))
    res = mesolve(H_loc, psi0_loc, tl, options={'store_states': True})
    S_loc = np.array([entropy_vn(st.ptrace(1), 2) for st in res.states])
    win = (tl > 0.3 * t_rev_loc) & (tl < 0.7 * t_rev_loc)
    k = int(np.argmin(np.where(win, S_loc, np.inf)))
    print(f"  n_bar = {nb:2d}: S_min = {S_loc[k]:.3f} at t = {tl[k] / t_rev_loc:.3f} t_r, "
          f"S(t_r/2) = {S_loc[int(np.argmin(np.abs(tl - 0.5 * t_rev_loc)))]:.3f}")

    ax.plot(tl * g, S_loc, color=COLORS['red'], lw=0.6)
    ax.axhline(1.0, color='0.6', lw=0.4)
    time_guides(ax, t_rev_loc * g)
    ax.set_xlim(0, tl[-1] * g)
    ax.set_ylim(-0.05, 1.1)
    ax.set_yticks([0, 0.5, 1])
    ax.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
    panel_label(ax, letter, rf'$\bar{{n}} = {nb}$')
for ax in axes7.flat:               # every panel has its own time range
    ax.set_xlabel(r'$gt$')
for ax in axes7[:, 0]:
    ax.set_ylabel(r'$S(\rho_{\rm atom})$ (bit)')
save(fig7, PAPER_FIG / 'fig_entropy_nbar_scaling')


# ---------------------------------------------------------------------------
# Fig 8: cavity loss kappa (Lindblad): inversion, reduced entropy, negativity
# ---------------------------------------------------------------------------
kappa_values = [0.0, 0.01, 0.05, 0.1, 0.2]
styles = ['-', (0, (6, 1.5)), (0, (3, 1.5)), '-.', ':']      # also distinguishable in greyscale
colors8 = ['k'] + sequential_colors(len(kappa_values) - 1, 'plasma', 0.12, 0.85)
tlist_diss = np.linspace(0, 30.0, 1801)        # to gt = 30, the time quoted in the text
ih = int(np.argmin(np.abs(tlist_diss - 0.5 * t_revival)))
fig8, (ax8a, ax8b, ax8c) = plt.subplots(3, 1, figsize=(FULL_W, 5.6), sharex=True)

for kap, color, ls in zip(kappa_values, colors8, styles):
    d = run_jc('coherent', STATES['coherent'][1], tlist_diss, kappa=kap)
    lbl = '$0$' if kap == 0 else rf'${kap:.2f}$'       # two decimals, as in the text and tables
    ax8a.plot(tlist_diss * g, d['inversion'], lw=0.7, color=color, ls=ls, label=lbl)
    ax8b.plot(tlist_diss * g, d['entropy'], lw=1.0, color=color, ls=ls)
    ax8c.plot(tlist_diss * g, d['logneg'], lw=1.0, color=color, ls=ls)
    print(f"kappa/g = {kap:.2f}: S(t_r/2) = {d['entropy'][ih]:.3f}, E_N(t_r/2) = {d['logneg'][ih]:.3f}, "
          f"S(t_end) = {d['entropy'][-1]:.3f}, E_N(t_end) = {d['logneg'][-1]:.3f}, "
          f"<sz>(t_end) = {d['inversion'][-1]:+.3f}")

ax8a.set_ylabel(r'$\langle \sigma_z \rangle$')
ax8a.set_ylim(-1.1, 1.1)
ax8a.set_yticks([-1, 0, 1])
ax8a.axhline(0, color='0.6', lw=0.4)
handles, labels8 = ax8a.get_legend_handles_labels()
ax8b.legend(handles, labels8, title=r'$\kappa/g$', ncol=5, loc='lower right', columnspacing=1.2,
            handlelength=1.6, borderaxespad=0.4)
ax8b.set_ylabel(r'$S(\rho_{\rm atom})$ (bit)')
ax8c.set_ylabel(r'$E_{\mathcal{N}}$')
ax8c.set_xlabel(r'$gt$')
for ax, letter in zip((ax8a, ax8b, ax8c), 'abc'):
    time_guides(ax, t_revival * g, labels=(letter == 'a'))
    panel_label(ax, letter)
    ax.set_xlim(0, tlist_diss[-1] * g)
for ax in (ax8b, ax8c):
    ax.set_ylim(-0.05, 1.1)
    ax.set_yticks([0, 0.5, 1])
    ax.axhline(1.0, color='0.6', lw=0.4)
save(fig8, PAPER_FIG / 'fig_dissipative_entanglement')
