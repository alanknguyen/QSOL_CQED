# sim_wigner_evolution.py
# Wigner function time evolution in the Jaynes-Cummings model.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Physics, Boston University
# Advised by Prof. Luca Dal Negro,
#
# Generates: fig_inversion_snapshots, fig_wigner_evolution,
#            fig_cat_state_detail, fig_wigner_decoherence
#
# Initial state: |e> x |alpha>, resonant (Delta = 0).
# The cavity field evolves from a coherent state into a Schrodinger
# cat state at t = t_r/2. Wigner negativity quantifies non-classicality.

import numpy as np
import matplotlib.pyplot as plt
from qutip import (
    basis, tensor, destroy, sigmaz, sigmam, sigmap,
    mesolve, wigner, ptrace, coherent, Qobj
)

# -- Parameters --
N_cav = 35
g = 1.0
omega_c = 0.0
omega_a = 0.0       # resonant
alpha = np.sqrt(10.0)  # <n> = 10
n_bar = np.abs(alpha)**2
t_collapse = 1.0 / g
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# -- Joint Hilbert space operators --
a  = tensor(destroy(N_cav), Qobj(np.eye(2)))
sm = tensor(Qobj(np.eye(N_cav)), sigmam())
sp = tensor(Qobj(np.eye(N_cav)), sigmap())
sz = tensor(Qobj(np.eye(N_cav)), sigmaz())

# H_JC = omega_c a^dag a + (omega_a/2) sigma_z + g(a^dag sigma^- + a sigma^+)
H_JC = omega_c * a.dag() * a + 0.5 * omega_a * sz + g * (a.dag() * sm + a * sp)

# -- Initial state: atom excited, field coherent --
psi0 = tensor(coherent(N_cav, alpha), basis(2, 0))

# -- Time grid --
# Wigner snapshots at six characteristic times
t_snapshots = [
    0.0,
    0.5 * t_collapse,
    2.0 * t_collapse,
    0.5 * t_revival,      # cat state forms here
    0.75 * t_revival,
    t_revival,             # first revival
]
snapshot_labels = [
    r'$t = 0$' + '\n(coherent state)',
    r'$t = 0.5\,t_c$' + '\n(early Rabi)',
    r'$t = 2\,t_c$' + '\n(collapse onset)',
    r'$t = 0.5\,t_r$' + '\n(cat state)',
    r'$t = 0.75\,t_r$' + '\n(pre-revival)',
    r'$t = t_r$' + '\n(first revival)',
]

tlist_dense = np.linspace(0, 1.5 * t_revival, 2000)
tlist_all = np.array(sorted(set(list(tlist_dense) + t_snapshots)))

# -- Solve --
print(f"Solving JC: <n>={n_bar:.0f}, g={g}, N_cav={N_cav}")
print(f"  t_collapse={t_collapse:.3f}, t_revival={t_revival:.3f}")

result = mesolve(H_JC, psi0, tlist_all, [], [sz], options={'store_states': True})
inversion = np.array(result.expect[0])


# Fig 1: Atomic inversion <sigma_z>(t) with Wigner snapshot markers
fig1, ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(tlist_all * g, inversion, 'b-', linewidth=0.6, alpha=0.9)
ax1.set_xlabel(r'$gt$', fontsize=14)
ax1.set_ylabel(r'$\langle \sigma_z \rangle$', fontsize=14)
ax1.set_title(rf'Atomic Inversion ($\bar{{n}} = {n_bar:.0f}$, $\Delta = 0$)', fontsize=13)
ax1.set_ylim(-1.1, 1.1)
ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')

for i, ts in enumerate(t_snapshots):
    idx = np.argmin(np.abs(tlist_all - ts))
    color = plt.cm.plasma(i / len(t_snapshots))
    ax1.axvline(ts * g, color=color, linewidth=1.2, linestyle=':', alpha=0.7)
    ax1.plot(ts * g, inversion[idx], 'o', color=color, markersize=8, zorder=5)

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
fig1.savefig('fig_inversion_snapshots.png', dpi=200, bbox_inches='tight')
fig1.savefig('fig_inversion_snapshots.pdf', bbox_inches='tight')
print("Saved: fig_inversion_snapshots")


# Fig 2: Wigner function W(x,p) of the reduced cavity field
# at the six snapshot times
xvec = np.linspace(-7, 7, 200)

fig2, axes = plt.subplots(2, 3, figsize=(14, 9))
axes_flat = axes.flatten()

print("Computing Wigner functions...")
for i, ts in enumerate(t_snapshots):
    idx = np.argmin(np.abs(tlist_all - ts))
    psi_t = result.states[idx]
    rho_field = ptrace(psi_t, 0)   # trace out atom, keep field
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

    # Wigner negativity: integral of |W| - 1 (= 0 for classical states)
    W_neg = np.sum(W[W < 0]) * (xvec[1] - xvec[0])**2
    purity = (rho_field * rho_field).tr()
    print(f"  t={ts:.3f}: negativity={abs(W_neg):.4f}, purity={purity:.4f}")

fig2.suptitle(
    r'Wigner Function Evolution of Cavity Field'
    + '\n' + rf'($\bar{{n}}={n_bar:.0f}$, $g={g}$, resonant)',
    fontsize=13, y=1.02
)
fig2.tight_layout()
fig2.savefig('fig_wigner_evolution.png', dpi=200, bbox_inches='tight')
fig2.savefig('fig_wigner_evolution.pdf', bbox_inches='tight')
print("Saved: fig_wigner_evolution")


# Fig 3: Cat state at t = t_r/2, with W(x, p=0) cross-section
cat_idx = np.argmin(np.abs(tlist_all - 0.5 * t_revival))
rho_cat = ptrace(result.states[cat_idx], 0)
W_cat = wigner(rho_cat, xvec, xvec)

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5))

wlim = np.max(np.abs(W_cat))
im3 = ax3a.contourf(xvec, xvec, W_cat, levels=100, cmap='RdBu_r',
                      vmin=-wlim, vmax=wlim)
ax3a.set_xlabel(r'$x$', fontsize=13)
ax3a.set_ylabel(r'$p$', fontsize=13)
ax3a.set_title(r'Wigner function at $t = t_r/2$ (cat state)', fontsize=12)
ax3a.set_aspect('equal')
fig3.colorbar(im3, ax=ax3a)

# Cross-section at p = 0: negative values prove non-classicality
mid_idx = len(xvec) // 2
ax3b.plot(xvec, W_cat[mid_idx, :], 'b-', linewidth=1.5, label=r'$W(x, p=0)$')
ax3b.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax3b.fill_between(xvec, W_cat[mid_idx, :], 0,
                   where=(W_cat[mid_idx, :] < 0),
                   color='red', alpha=0.3, label='Negative region')
ax3b.set_xlabel(r'$x$', fontsize=13)
ax3b.set_ylabel(r'$W(x, 0)$', fontsize=13)
ax3b.set_title(r'Cross-section at $p = 0$', fontsize=11)
ax3b.legend(fontsize=11)

fig3.tight_layout()
fig3.savefig('fig_cat_state_detail.png', dpi=200, bbox_inches='tight')
fig3.savefig('fig_cat_state_detail.pdf', bbox_inches='tight')
print("Saved: fig_cat_state_detail")


# Fig 4: Decoherence of the cat state under cavity loss kappa.
#         Lindblad: d rho/dt = -i[H,rho] + kappa D[a] rho
kappa_values = [0.0, 0.02 * g, 0.05 * g, 0.1 * g]

fig4, axes4 = plt.subplots(1, 4, figsize=(16, 4))

print("\nDissipative Wigner functions at t = t_r/2:")
for j, kap in enumerate(kappa_values):
    c_ops = [np.sqrt(kap) * a] if kap > 0 else []
    tlist_short = np.linspace(0, 0.5 * t_revival, 500)
    res = mesolve(H_JC, psi0, tlist_short, c_ops, [],
                  options={'store_states': True})
    rho_f = ptrace(res.states[-1], 0)
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
    print(f"  kappa/g={kap/g:.2f}: negativity={abs(W_neg_d):.4f}, purity={purity_d:.4f}")

fig4.suptitle(
    rf'Decoherence of Cat State at $t = t_r/2$ ($\bar{{n}}={n_bar:.0f}$)',
    fontsize=13, y=1.04
)
fig4.tight_layout()
fig4.savefig('fig_wigner_decoherence.png', dpi=200, bbox_inches='tight')
fig4.savefig('fig_wigner_decoherence.pdf', bbox_inches='tight')
print("Saved: fig_wigner_decoherence")