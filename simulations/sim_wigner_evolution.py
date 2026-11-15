# sim_wigner_evolution.py
# Wigner function time evolution in the Jaynes-Cummings model.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Physics, Boston University
# Began as a course project taught by Prof. Luca Dal Negro
#
# Generates (paper/figures/): fig_inversion_snapshots, fig_wigner_evolution,
#                             fig_cat_state_detail, fig_wigner_decoherence
#
# Initial state: |e> x |alpha>, resonant (Delta = 0).
# The cavity field splits into two counter-rotating components.  At
# t = t_r/2 they sit on opposite sides of phase space (on the p axis for real
# alpha), the atom and field approximately disentangle, and the field is left
# in a Schrodinger cat state (Gea-Banacloche 1990).  Wigner negativity
#   delta = int |W| dx dp - 1        (Kenfack & Zyczkowski 2004)
# quantifies non-classicality; delta = 2 * (integrated negative volume).

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from qutip import basis, tensor, destroy, sigmaz, sigmam, qeye, mesolve, wigner, coherent

from plotstyle import (use_style, save, panel_label, wigner_image, FULL_W, W_MAX,
                       COLORS, PAPER_FIG)
from utils import wigner_negativity, cat_fidelity, photon_parity

use_style()

# -- Parameters --
N_cav = 35
g = 1.0
n_bar = 10.0
alpha = np.sqrt(n_bar)
t_collapse = 1.0 / g                       # operational collapse time, t_c = 1/g
t_revival = 2 * np.pi * np.sqrt(n_bar) / g

# -- Operators and Hamiltonian (resonant, interaction picture) --
a = tensor(destroy(N_cav), qeye(2))
sm = tensor(qeye(N_cav), sigmam())
sz = tensor(qeye(N_cav), sigmaz())
H_JC = g * (a.dag() * sm + a * sm.dag())

# -- Initial state: atom excited, field coherent --
psi0 = tensor(coherent(N_cav, alpha), basis(2, 0))

# -- Six snapshot times; panel letters are shared by Figs. 1 and 2 --
snapshots = [
    ('a', 0.0,               r'$t = 0$'),
    ('b', 0.5 * t_collapse,  r'$t = 0.5\,t_c$'),
    ('c', 2.0 * t_collapse,  r'$t = 2\,t_c$'),
    ('d', 0.5 * t_revival,   r'$t = t_r/2$'),
    ('e', 0.75 * t_revival,  r'$t = 0.75\,t_r$'),
    ('f', t_revival,         r'$t = t_r$'),
]
t_snap = [s[1] for s in snapshots]
tlist = np.array(sorted(set(np.linspace(0, 1.5 * t_revival, 3000)) | set(t_snap)))

print(f"Solving JC: <n> = {n_bar:.0f}, g = {g}, N_cav = {N_cav}")
print(f"  t_collapse = {t_collapse:.3f}, t_revival = {t_revival:.3f}")
result = mesolve(H_JC, psi0, tlist, e_ops=[sz], options={'store_states': True})
inversion = np.array(result.expect[0])
state_at = lambda t: result.states[int(np.argmin(np.abs(tlist - t)))]

# Phase-space grid: 281 x 281 points on [-7, 7]^2, spacing 0.05, i.e. 14 points
# per interference fringe (fringe period pi / sqrt(2 n_bar) = 0.70).
xvec = np.linspace(-7, 7, 281)


# ---------------------------------------------------------------------------
# Fig 1: atomic inversion with the snapshot times marked
# ---------------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(FULL_W, 2.35))
ax1.plot(tlist * g, inversion, color=COLORS['blue'], lw=0.8)
ax1.axhline(0, color='0.6', lw=0.4)

# Collapse envelope from linearising sqrt(n+1) about n_bar and summing over the
# Poisson weights: exp[-g^2 t^2 n_bar / 2(n_bar + 1)], which tends to the Gaussian
# exp(-g^2 t^2 / 2) for large n_bar.  With the factor n_bar/(n_bar+1) it bounds
# every extremum of the oscillation; without it the curve overshoots it visibly.
t_env = np.linspace(0, 4.5, 200)
env = np.exp(-(g * t_env) ** 2 * n_bar / (2 * (n_bar + 1)))
ax1.plot(t_env, env, color=COLORS['grey'], lw=0.7, ls='--', label='collapse envelope')
ax1.plot(t_env, -env, color=COLORS['grey'], lw=0.7, ls='--')

trans = ax1.get_xaxis_transform()
for k, (letter, ts, _) in enumerate(snapshots):
    idx = int(np.argmin(np.abs(tlist - ts)))
    ax1.axvline(ts * g, color='0.45', lw=0.5, ls=':')
    ax1.plot(ts * g, inversion[idx], 'o', ms=4.5, mfc='white', mec='k', mew=0.8,
             zorder=5, clip_on=False)
    # (a), (b), (c) are close together: stagger their labels
    ax1.text(ts * g, 1.10 if letter == 'b' else 1.02, f'({letter})', transform=trans,
             ha='center', va='bottom', fontsize=8)
ax1.text(0.5 * t_revival * g + 0.25, -0.98, r'$t_r/2$', fontsize=8, va='bottom')
ax1.text(t_revival * g + 0.25, -0.98, r'$t_r$', fontsize=8, va='bottom')

ax1.set_xlabel(r'$gt$')
ax1.set_ylabel(r'$\langle \sigma_z \rangle$')
ax1.set_xlim(0, 1.5 * t_revival * g)
ax1.set_ylim(-1.08, 1.08)
ax1.set_yticks([-1, -0.5, 0, 0.5, 1])
ax1.legend(loc='upper left', bbox_to_anchor=(0.16, 1.0))
save(fig1, PAPER_FIG / 'fig_inversion_snapshots')


# ---------------------------------------------------------------------------
# Fig 2: Wigner function of the reduced cavity field at the six snapshot times
# ---------------------------------------------------------------------------
fig2, axes = plt.subplots(2, 3, figsize=(FULL_W, 3.95), sharex=True, sharey=True)
phi = np.linspace(0, 2 * np.pi, 400)
radius = np.sqrt(2 * n_bar)        # |alpha| = sqrt(n_bar)  <->  x^2 + p^2 = 2 n_bar

print("Wigner functions of the reduced field:")
for ax, (letter, ts, label) in zip(axes.flat, snapshots):
    rho_field = state_at(ts).ptrace(0)
    W = wigner(rho_field, xvec, xvec)
    im = wigner_image(ax, W, xvec)
    ax.plot(radius * np.cos(phi), radius * np.sin(phi), color='0.35', lw=0.4, ls=(0, (4, 4)))
    delta = wigner_negativity(W, xvec)
    purity = (rho_field * rho_field).tr().real
    panel_label(ax, letter, label)
    ax.text(0.04, 0.04, rf'$\delta = {max(delta, 0):.2f}$', transform=ax.transAxes,
            fontsize=8, va='bottom')
    ax.set_xticks([-5, 0, 5]); ax.set_yticks([-5, 0, 5])
    print(f"  ({letter}) t = {ts:6.3f}: delta = {delta:.4f}, purity = {purity:.4f}, "
          f"W in [{W.min():.3f}, {W.max():.3f}]")
for ax in axes[-1]:
    ax.set_xlabel(r'$x$')
for ax in axes[:, 0]:
    ax.set_ylabel(r'$p$')
cb = fig2.colorbar(im, ax=axes, shrink=0.85, aspect=30, pad=0.02,
                   ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb.set_label(r'$W(x,p)$')
save(fig2, PAPER_FIG / 'fig_wigner_evolution')


# ---------------------------------------------------------------------------
# Fig 3: the cat state at t = t_r/2 and its cut along p = 0
# ---------------------------------------------------------------------------
rho_cat = state_at(0.5 * t_revival).ptrace(0)
W_cat = wigner(rho_cat, xvec, xvec)
delta_cat = wigner_negativity(W_cat, xvec)
x_line = np.linspace(-7, 7, 2801)
W_line = wigner(rho_cat, x_line, np.array([0.0]))[0]       # exactly on p = 0

peaks, _ = find_peaks(W_line, height=0.05)
troughs, _ = find_peaks(-W_line, height=0.05)
spacing = np.mean(np.diff(x_line[peaks]))
F_cat, beta_fit, theta_fit = cat_fidelity(rho_cat, n_bar, return_params=True)
parity = photon_parity(rho_cat)
print(f"Cat state at t_r/2: delta = {delta_cat:.4f}, purity = {(rho_cat * rho_cat).tr().real:.4f}")
print(f"  W(x, 0): min = {W_line.min():.4f} at x = {x_line[np.argmin(W_line)]:.3f}, "
      f"max = {W_line.max():.4f}; mean fringe spacing = {spacing:.3f} "
      f"(pi/sqrt(2 n_bar) = {np.pi / np.sqrt(2 * n_bar):.3f})")
print(f"  best-fit cat: F = {F_cat:.3f}, beta = {beta_fit:.3f}, theta = {theta_fit / np.pi:.2f} pi; "
      f"photon-number parity = {parity:+.3f}")

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(FULL_W, 2.75),
                                  gridspec_kw={'width_ratios': [1, 1.3]})
im3 = wigner_image(ax3a, W_cat, xvec)
ax3a.axhline(0, color='k', lw=0.5, ls=(0, (4, 3)))
ax3a.set_xlim(-6, 6); ax3a.set_ylim(-6, 6)
ax3a.set_xticks([-4, 0, 4]); ax3a.set_yticks([-4, 0, 4])
ax3a.set_xlabel(r'$x$'); ax3a.set_ylabel(r'$p$')
panel_label(ax3a, 'a', rf'$W(x,p)$ at $t_r/2$, $\delta = {delta_cat:.2f}$')
cb3 = fig3.colorbar(im3, ax=ax3a, shrink=0.9, aspect=25, pad=0.03,
                    ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb3.set_label(r'$W(x,p)$')

ax3b.axhline(0, color='0.6', lw=0.4)
# negative regions in the blue of the map in (a)
ax3b.fill_between(x_line, W_line, 0, where=W_line < 0, color=COLORS['blue'], alpha=0.30,
                  lw=0, label=r'$W < 0$')
ax3b.plot(x_line, W_line, color=COLORS['blue'], lw=1.0)
# fringe-spacing marker between the two tallest neighbouring maxima
k = int(np.argmax(W_line[peaks]))
k2 = k + 1 if k + 1 < len(peaks) and (k == 0 or W_line[peaks[k + 1]] > W_line[peaks[k - 1]]) else k - 1
xa, xb = sorted((x_line[peaks[k]], x_line[peaks[k2]]))
y_arrow = W_line.max() + 0.035
ax3b.annotate('', xy=(xa, y_arrow), xytext=(xb, y_arrow),
              arrowprops=dict(arrowstyle='<->', lw=0.6, shrinkA=0, shrinkB=0))
ax3b.text(xb + 0.15, y_arrow, r'$\pi/\sqrt{2\bar{n}}$', fontsize=8, va='center', ha='left')
ax3b.set_xlim(-4.5, 5.5)
ax3b.set_ylim(-0.31, 0.34)
ax3b.set_xlabel(r'$x$'); ax3b.set_ylabel(r'$W(x,\,p=0)$')
ax3b.legend(loc='lower right')
panel_label(ax3b, 'b', r'cut along $p = 0$')
# give (b) exactly the height of the square map in (a), so frames and titles line up
fig3.canvas.draw()
fig3.set_layout_engine('none')
_pa, _pb = ax3a.get_position(), ax3b.get_position()
ax3b.set_position([_pb.x0, _pa.y0, _pb.width, _pa.height])
save(fig3, PAPER_FIG / 'fig_cat_state_detail')


# ---------------------------------------------------------------------------
# Fig 4: decoherence of the cat under cavity loss kappa
#        d rho/dt = -i[H, rho] + kappa D[a] rho
# ---------------------------------------------------------------------------
kappa_values = [0.0, 0.02, 0.05, 0.10]
fig4, axes4 = plt.subplots(1, 4, figsize=(FULL_W, 2.05), sharex=True, sharey=True)

print("Reduced field at t = t_r/2 with cavity loss:")
for ax, letter, kap in zip(axes4, 'abcd', kappa_values):
    c_ops = [np.sqrt(kap * g) * a] if kap > 0 else []
    res = mesolve(H_JC, psi0, [0, 0.5 * t_revival], c_ops=c_ops)
    rho_f = res.states[-1].ptrace(0)
    W_d = wigner(rho_f, xvec, xvec)
    im4 = wigner_image(ax, W_d, xvec)
    delta_d = wigner_negativity(W_d, xvec)
    purity_d = (rho_f * rho_f).tr().real
    panel_label(ax, letter, rf"$\kappa/g = {'0' if kap == 0 else f'{kap:.2f}'}$", y=1.17)
    ax.text(0.0, 1.03, rf'$\delta = {delta_d:.3f}$', transform=ax.transAxes, fontsize=8.5,
            va='bottom')
    ax.set_xlim(-6.5, 6.5); ax.set_ylim(-6.5, 6.5)
    ax.set_xticks([-5, 0, 5]); ax.set_yticks([-5, 0, 5])
    ax.set_xlabel(r'$x$')
    print(f"  kappa/g = {kap:.2f}: delta = {delta_d:.4f}, purity = {purity_d:.4f}, "
          f"W in [{W_d.min():.3f}, {W_d.max():.3f}]")
axes4[0].set_ylabel(r'$p$')
cb4 = fig4.colorbar(im4, ax=axes4, shrink=0.95, aspect=18, pad=0.015, ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb4.set_label(r'$W(x,p)$')
save(fig4, PAPER_FIG / 'fig_wigner_decoherence')
