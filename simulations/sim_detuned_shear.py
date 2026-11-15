# sim_detuned_shear.py
# Rotation and shear off resonance: detuning, heralded cat states and cavity loss.
#
# With an atom-cavity detuning Delta the Rabi frequencies are
#   Omega_n = sqrt(Delta^2 + 4 g^2 (n+1)).
# Write x = Delta / (2 g sqrt(n_bar)), sin(theta) = 1/sqrt(1 + x^2) and
# n_crit = Delta^2 / 4 g^2 (the "critical photon number" of the dispersive
# regime), so that x^2 = n_crit / n_bar.  Expanding Omega_n about n_bar, exactly
# as sim_rotation_shear.py does on resonance, predicts
#   rotation   phi = g^2 t / Omega : the two field components sit at angles -/+ phi;
#   shear      s = phi sin^2(theta) = phi n_bar / (n_bar + n_crit);
#   purity     (1 + cos^2 theta)/2 + sin^4(theta) sin^2(phi) / 2   between collapse and revival;
#   heralding  finding the atom in |g> leaves a BALANCED superposition of the two
#              components with probability sin^2(theta)/2, whose fidelity with an
#              ideal superposition of two coherent states is (1 + s^2/4)^(-1/2)
#              for n_bar -> infinity (the overlap of each component with the coherent
#              state at its centre).  At finite n_bar the best two coherent states sit
#              slightly off the centres and the fidelity is a little higher: by 2% at
#              n_bar = 10, 0.2% at n_bar = 100 for phi = pi/2;
#   revivals   amplitude sin^2(theta) (1 + pi^2 k^2 sin^4 theta)^(-1/4) at t_k = k pi Omega / g^2;
#   loss       the coherence between the components decays as exp(-Gamma),
#              Gamma = (kappa n_bar Omega / g^2) (phi - sin(2 phi)/2),
#              so the heralded-cat fidelity is (1 + s^2/4)^(-1/2) (1 + exp(-Gamma))/2.
# For a cat of fixed size D = |beta_1 - beta_2| = 2 sqrt(n_bar) sin(phi):
#   * a LARGER resonant field needs a smaller rotation, so the shear vanishes
#     (F -> (1 + D^2/16 n_bar)^(-1/2)) and Gamma falls from pi kappa D^3 / 8g
#     (phi = pi/2, the half-revival cat) to kappa D^3 / 6g;
#   * at fixed n_bar, detuning lowers the shear but lengthens the interaction;
#     it raises the fidelity only if  Gamma_0 / (exp(Gamma_0) + 1) < 2 phi^2 / (4 + phi^2),
#     Gamma_0 being the resonant exponent (kappa D n_bar / g < 3/2 for small phi),
#     and it always lowers the heralding probability.
# Every statement is checked here against the exact amplitudes (no loss) or
# against the master equation (loss).
#
# Generates: paper/figures/fig_detuned_shear  and prints the comparison tables.
#
# Nguyen Khoi Nguyen (Alan), Boston University

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.linalg import eigh
from scipy.optimize import differential_evolution, minimize
from scipy.special import gammaln
from scipy.stats import poisson

from plotstyle import (use_style, save, panel_label, wigner_image, sequential_colors,
                       FULL_W, COLORS, PAPER_FIG)

use_style()
g = 1.0
PI = np.pi


# ---------------------------------------------------------------------------
# Exact solution with detuning (no loss)
# ---------------------------------------------------------------------------
def fock_grid(nb):
    """Fock numbers and the coherent-state amplitudes c_n; the neglected Poisson weight is < 1e-30."""
    n = np.arange(int(nb + 14 * np.sqrt(nb) + 40))
    return n, np.sqrt(poisson.pmf(n, nb))


def angles(nb, x):
    """Detuning, mean Rabi frequency and sin(theta) for x = Delta / (2 g sqrt(nb))."""
    Delta = 2 * g * np.sqrt(nb) * x
    Om = np.sqrt(Delta ** 2 + 4 * g * g * nb)
    return Delta, Om, 2 * g * np.sqrt(nb) / Om


def amplitudes(nb, x, phi):
    """Exact amplitudes of |e,n> and |g,n> for |e> x |alpha> at the time when the rotation angle is phi.

    H = (Delta/2) sigma_z + g (a^dag sigma_- + a sigma_+) couples |e,n> to |g,n+1> only, so
      c_e,n   = c_n [cos(Omega_n t/2) - i (Delta/Omega_n) sin(Omega_n t/2)],
      c_g,n+1 = -i c_n (2 g sqrt(n+1)/Omega_n) sin(Omega_n t/2).
    """
    Delta, Om, _ = angles(nb, x)
    t = phi * Om / g ** 2
    n, c = fock_grid(nb)
    O = np.sqrt(Delta ** 2 + 4 * g * g * (n + 1))
    ce = c * (np.cos(O * t / 2) - 1j * (Delta / O) * np.sin(O * t / 2))
    cg = np.zeros(len(n), complex)
    cg[1:] = (-1j * c * (2 * g * np.sqrt(n + 1) / O) * np.sin(O * t / 2))[:-1]
    return n, ce, cg


def shear_exact(nb, x, phi):
    """Shear of one branch |Phi_+> = sum_n c_n exp(i Omega_n t/2)|n>, read off its covariance matrix.

    A sheared coherent state has det V = 1/4 and Tr V = 1 + s^2/2, i.e. s^2 = 4 (<n> - |<a>|^2).
    """
    Delta, Om, _ = angles(nb, x)
    t = phi * Om / g ** 2
    n, c = fock_grid(nb)
    psi = c * np.exp(0.5j * np.sqrt(Delta ** 2 + 4 * g * g * (n + 1)) * t)
    a_mean = np.sum(np.conj(psi[:-1]) * psi[1:] * np.sqrt(n[1:]))
    return 2 * np.sqrt(np.sum(n * abs(psi) ** 2) - abs(a_mean) ** 2)


def coherent_amplitudes(n, beta):
    return np.exp(-abs(beta) ** 2 / 2 + n * np.log(complex(beta)) - 0.5 * gammaln(n + 1))


def cat_overlap(state, n):
    """Function z = (Re beta1, Im beta1, Re beta2, Im beta2) -> largest fidelity of `state` (vector or
    density matrix on the Fock numbers n) with a superposition a1|beta1> + a2|beta2>; the best
    weights a1, a2 solve a 2x2 generalised eigenproblem."""
    pure = state.ndim == 1

    def F(z):
        K = np.stack([coherent_amplitudes(n, complex(z[0], z[1])),
                      coherent_amplitudes(n, complex(z[2], z[3]))], axis=1)
        if pure:
            v = K.conj().T @ state
            M = np.outer(v, v.conj())
        else:
            M = K.conj().T @ state @ K
        return eigh(M, K.conj().T @ K, eigvals_only=True)[-1]
    return F


def cat_fidelity(state, n, b1, b2, step=0.3):
    """Largest fidelity of a normalised state with any superposition a1|beta1> + a2|beta2>, the
    amplitudes beta1, beta2 being optimised by Nelder-Mead starting from b1, b2.

    The starting simplex is given explicitly, with side `step` in every direction.  SciPy's default
    moves each coordinate by 5% of its value, and a coordinate that is exactly zero by only 2.5e-4:
    at phi = pi/2, where Re b1 = Re b2 = 0, that simplex is flat in two directions and the search
    stops short of the maximum (0.780 instead of 0.802 for n_bar = 10 on resonance).  A second pass
    from the optimum with a smaller simplex guards against premature termination.
    """
    F = cat_overlap(state, n)
    z = np.array([b1.real, b1.imag, b2.real, b2.imag])
    for side in (step, step / 10):
        res = minimize(lambda v: -F(v), z, method='Nelder-Mead',
                       options={'initial_simplex': np.vstack([z, z + side * np.eye(4)]),
                                'xatol': 1e-8, 'fatol': 1e-13, 'maxiter': 20000, 'maxfev': 20000})
        z = res.x
    return -res.fun


def heralded(nb, x, phi):
    """No loss: (probability of finding |g>, purity of the unheralded field, fidelity of the heralded cat)."""
    n, ce, cg = amplitudes(nb, x, phi)
    pe, pg = np.vdot(ce, ce).real, np.vdot(cg, cg).real
    purity = pe ** 2 + pg ** 2 + 2 * abs(np.vdot(ce, cg)) ** 2
    al = np.sqrt(nb)
    return pg, purity, cat_fidelity(cg / np.sqrt(pg), n, al * np.exp(-1j * phi), al * np.exp(1j * phi))


# ---------------------------------------------------------------------------
# Master equation with cavity loss
# ---------------------------------------------------------------------------
def heralded_with_loss(nb, x, kappa, phi, width=7.0):
    """(probability of |g>, fidelity of the heralded field with the best cat, purity of the heralded field).

    The Fock space is a window n_bar -/+ width sqrt(n_bar) (plus a margin): photons
    are only lost, a few per cent of them at most, so nothing reaches the edges.
    """
    Delta, Om, _ = angles(nb, x)
    lo = max(0, int(nb - width * np.sqrt(nb) - 15))
    hi = int(nb + width * np.sqrt(nb) + 15)
    N = hi - lo
    a = qt.tensor(qt.destroy(N, offset=lo), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    sz = qt.tensor(qt.qeye(N), qt.sigmaz())
    H = 0.5 * Delta * sz + g * (a.dag() * sm + a * sm.dag())
    psi0 = qt.tensor(qt.coherent(N, np.sqrt(nb), offset=lo, method='analytic'), qt.basis(2, 0))
    t = phi * Om / g ** 2
    opts = {'atol': 1e-10, 'rtol': 1e-8, 'nsteps': 2000000, 'store_states': True}
    final = qt.mesolve(H, psi0, [0.0, t], c_ops=[np.sqrt(kappa) * a] if kappa > 0 else [],
                       options=opts).states[-1]
    rho = final.full() if final.isoper else np.outer(final.full().ravel(), final.full().ravel().conj())
    rho_g = rho[1::2, 1::2]                      # atom index 1 = |g>
    pg = np.trace(rho_g).real
    rho_g = rho_g / pg
    amp = np.sqrt(nb) * np.exp(-kappa * t / 2)
    F = cat_fidelity(rho_g, np.arange(lo, hi), amp * np.exp(-1j * phi), amp * np.exp(1j * phi))
    return pg, F, np.trace(rho_g @ rho_g).real


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------
def shear(x, phi):
    return phi / (1 + x * x)


def loss_exponent(nb, x, kappa, phi):
    return (2 * kappa * nb ** 1.5 / g) * np.sqrt(1 + x * x) * (phi - 0.5 * np.sin(2 * phi))


def fidelity_model(nb, x, kappa, phi):
    return (1 + shear(x, phi) ** 2 / 4) ** -0.5 * 0.5 * (1 + np.exp(-loss_exponent(nb, x, kappa, phi)))


def model_ellipses(nb, phi, s, level=2.0):
    """1/e contours of the two Gaussian components: centre sqrt(2 nb) at angle -/+ phi, shear +/- s."""
    tau = np.linspace(0, 2 * PI, 200)
    out = []
    for sign in (+1, -1):
        ang = -sign * phi
        V_loc = 0.5 * np.array([[1.0, sign * s], [sign * s, 1 + s * s]])     # (radial, tangential)
        R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        Lc = np.linalg.cholesky(R @ V_loc @ R.T)
        centre = np.sqrt(2 * nb) * np.array([np.cos(ang), np.sin(ang)])
        out.append(np.sqrt(level) * Lc @ np.vstack([np.cos(tau), np.sin(tau)]) + centre[:, None])
    return out


# ---------------------------------------------------------------------------
# Tables: no loss
# ---------------------------------------------------------------------------
print("Detuned Jaynes-Cummings model, n_bar = 100, exact / closed form")
print("  x = Delta/(2 g sqrt(n_bar)), phi = rotation angle")
print("     x   phi/pi   field purity      P(atom in g)      heralded-cat fidelity")
for x in (0.0, 0.5, 1.0, 2.0):
    sth2 = 1 / (1 + x * x)
    for phi in (PI / 8, PI / 4, PI / 2):
        pg, pur, F = heralded(100, x, phi)
        pur_m = 0.5 * (2 - sth2) + 0.5 * sth2 ** 2 * np.sin(phi) ** 2
        print(f"   {x:3.1f}  {phi / PI:5.3f}   {pur:.4f} / {pur_m:.4f}   {pg:.4f} / {0.5 * sth2:.4f}   "
              f"{F:.4f} / {(1 + shear(x, phi) ** 2 / 4) ** -0.5:.4f}")

print("\nRevival amplitude |sum_n p_n sin^2(theta_n) exp(i Omega_n t_k - 2 pi i k n)| at t_k = k pi Omega/g^2,")
print("exact / sin^2(theta) (1 + pi^2 k^2 sin^4 theta)^(-1/4), n_bar = 100")
n100, c100 = fock_grid(100)
for x in (0.0, 0.5, 1.0, 2.0):
    Delta, Om, sth = angles(100, x)
    O = np.sqrt(Delta ** 2 + 4 * g * g * (n100 + 1))
    row = []
    for k in (1, 2, 3):
        Z = np.sum(c100 ** 2 * (4 * g * g * (n100 + 1) / O ** 2) * np.exp(1j * (O * k * PI * Om / g ** 2 - 2 * PI * k * n100)))
        row.append(f"k={k}: {abs(Z):.4f} / {sth ** 2 * (1 + (PI * k * sth ** 2) ** 2) ** -0.25:.4f}")
    print(f"   x = {x:3.1f}   " + "   ".join(row))

# the master-equation route must reproduce the exact amplitudes when kappa = 0
_pg, _, _F = heralded(25, 0.5, PI / 4)
_pg2, _F2, _ = heralded_with_loss(25, 0.5, 0.0, PI / 4)
assert abs(_pg - _pg2) < 1e-6 and abs(_F - _F2) < 1e-5, (_pg, _pg2, _F, _F2)

# the fidelity search must reach the global optimum, also at phi = pi/2 where two of the starting
# coordinates vanish: compare with a differential-evolution search over a box around the start
for _nb, _x in ((10, 0.0), (10, 1.0), (100, 0.0)):
    _n, _, _cg = amplitudes(_nb, _x, PI / 2)
    _psi = _cg / np.linalg.norm(_cg)
    _F = cat_fidelity(_psi, _n, -1j * np.sqrt(_nb), 1j * np.sqrt(_nb))
    _obj = cat_overlap(_psi, _n)
    _ref = differential_evolution(lambda v: -_obj(v), [(-1.6, 1.6), (-np.sqrt(_nb) - 1.6, -np.sqrt(_nb) + 1.6),
                                                       (-1.6, 1.6), (np.sqrt(_nb) - 1.6, np.sqrt(_nb) + 1.6)],
                                  seed=1, tol=1e-10, popsize=30)
    assert abs(_F + _ref.fun) < 1e-6, (_nb, _x, _F, -_ref.fun)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
D_CAT = 2 * np.sqrt(10.0)                 # |beta_1 - beta_2| of the half-revival cat of n_bar = 10
fig = plt.figure(figsize=(FULL_W, 4.9))
gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1.0])

# (a)-(c): the same cat made three ways
cases = [(10, 0.0, r'$\bar{n} = 10$, $\Delta = 0$'),
         (10, 1.0, r'$\bar{n} = 10$, $\Delta = 2g\sqrt{\bar{n}}$'),
         (100, 0.0, r'$\bar{n} = 100$, $\Delta = 0$')]
top_axes = []
print(f"\nWigner panels: cat size D = {D_CAT:.3f}")
for col, (nb, x, lab) in enumerate(cases):
    phi = np.arcsin(min(D_CAT / (2 * np.sqrt(nb)), 1.0))
    n, ce, cg = amplitudes(nb, x, phi)
    pg = np.vdot(cg, cg).real
    psi = cg / np.sqrt(pg)
    F = cat_fidelity(psi, n, np.sqrt(nb) * np.exp(-1j * phi), np.sqrt(nb) * np.exp(1j * phi))
    xc = np.sqrt(2 * nb) * np.cos(phi)
    pv = np.linspace(-7.5, 7.5, 301)
    xv = pv + xc
    W = qt.wigner(qt.Qobj(psi), xv, pv)
    ax = fig.add_subplot(gs[0, col])
    im = wigner_image(ax, W, xv, pvec=pv)
    for pts, ls in zip(model_ellipses(nb, phi, shear(x, phi)), ('-', '--')):
        ax.plot(pts[0], pts[1], color='k', lw=0.7, ls=ls)
    ax.set_xlim(xv[0], xv[-1]); ax.set_ylim(pv[0], pv[-1])
    ax.set_xticks([-5, 0, 5] if xc < 1 else [10, 15, 20]); ax.set_yticks([-5, 0, 5])
    ax.tick_params(labelleft=(col == 0))
    ax.set_xlabel(r'$x$')
    panel_label(ax, 'abc'[col], lab, fontsize=8)
    ax.text(0.04, 0.04, rf'$F = {F:.2f}$', transform=ax.transAxes, ha='left', va='bottom', fontsize=7.5)
    ax.text(0.96, 0.04, rf'$gt = {phi * angles(nb, x)[1] / g ** 2:.1f}$', transform=ax.transAxes, ha='right',
            va='bottom', fontsize=7.5)
    top_axes.append(ax)
    print(f"   n_bar = {nb:3d}, x = {x:.1f}: phi = {phi:.3f}, g t = {phi * angles(nb, x)[1] / g ** 2:5.2f}, "
          f"P_g = {pg:.3f}, F = {F:.4f} (closed form {(1 + shear(x, phi) ** 2 / 4) ** -0.5:.4f})")
top_axes[0].set_ylabel(r'$p$')
cb = fig.colorbar(im, ax=top_axes, shrink=0.9, aspect=24, pad=0.02, ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb.set_label(r'$W(x,p)$')

# (d): shear per unit rotation against n_crit / n_bar
ax = fig.add_subplot(gs[1, 0])
ratio = np.geomspace(1e-2, 1e2, 200)
ax.loglog(ratio, 1 / (1 + ratio), color='k', lw=0.9, label=r'$\bar{n}/(\bar{n}+n_{\rm crit})$')
ax.loglog(ratio[ratio > 1.5], 1 / ratio[ratio > 1.5], color='0.5', lw=0.7, ls='--', label=r'$\bar{n}/n_{\rm crit}$')
pts_ratio = np.geomspace(1e-2, 1e2, 9)
for nb, mk, col in zip((25, 100, 400), ('o', 's', '^'), sequential_colors(3, 'viridis')):
    s_ex = [shear_exact(nb, np.sqrt(r), PI / 2) / (PI / 2) for r in pts_ratio]
    ax.loglog(pts_ratio, s_ex, mk, ms=3.6, mfc='none', mec=col, mew=0.9, label=rf'$\bar{{n}} = {nb}$')
ax.text(0.0095, 0.55, 'resonant', fontsize=7, va='top')
ax.text(70, 1.25, 'dispersive\n(Kerr)', fontsize=7, ha='right', va='top', color='0.35')
ax.set_xlabel(r'$n_{\rm crit}/\bar{n} = \Delta^2/4g^2\bar{n}$')
ax.set_ylabel(r'shear per rotation, $s/\varphi$')
ax.set_xlim(7e-3, 1.4e2); ax.set_ylim(6e-3, 1.6)
ax.legend(loc='lower left', fontsize=7, handlelength=1.4, borderaxespad=0.3)
panel_label(ax, 'd')

# (e): fixed cat size, resonance, larger and larger fields, with loss
ax = fig.add_subplot(gs[1, 1])
kap_line = np.linspace(0, 0.02, 200)
kap_pts = np.array([0.0, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02])
print(f"\nFixed cat size D = {D_CAT:.2f} on resonance, with cavity loss: master equation / closed form")
for nb, mk, col, ls in zip((10, 25, 100), ('o', 's', '^'), sequential_colors(3, 'viridis'), ('-', '--', '-.')):
    phi = np.arcsin(min(D_CAT / (2 * np.sqrt(nb)), 1.0))
    ax.plot(kap_line, fidelity_model(nb, 0.0, kap_line, phi), color=col, ls=ls, lw=0.9)
    F_me = [heralded_with_loss(nb, 0.0, k, phi)[1] for k in kap_pts]
    ax.plot(kap_pts, F_me, mk, ms=3.6, mfc='none', mec=col, mew=0.9, label=rf'${nb}$')
    print(f"   n_bar = {nb:3d} (phi = {phi:.3f}): " + "  ".join(
        f"{k:.4f}: {f:.3f}/{fidelity_model(nb, 0.0, k, phi):.3f}" for k, f in zip(kap_pts, F_me)))
ax.plot(kap_line, 0.5 * (1 + np.exp(-kap_line * D_CAT ** 3 / (6 * g))), color='0.5', lw=0.7, ls=':',
        label=r'$\infty$')
ax.set_xlabel(r'$\kappa/g$'); ax.set_ylabel(r'heralded-cat fidelity $F$')
ax.set_xlim(-0.0008, 0.0208); ax.set_ylim(0.4, 1.02)
ax.set_xticks([0, 0.01, 0.02])
ax.legend(loc='upper right', fontsize=7, handlelength=1.4, borderaxespad=0.2, labelspacing=0.25,
          title=r'$\bar{n}$', title_fontsize=7.5)
panel_label(ax, 'e', r'$\Delta = 0$', fontsize=8)

# (f): fixed cat size and fixed field, detuning against loss
ax = fig.add_subplot(gs[1, 2])
nbF = 25
phiF = np.arcsin(D_CAT / (2 * np.sqrt(nbF)))
x_line = np.linspace(0, 2, 201)
x_pts = np.linspace(0, 2, 9)
kappas = (0.0, 0.002, 0.005, 0.01)
print(f"\nFixed cat size and field (n_bar = {nbF}, phi = {phiF:.3f}), fidelity against detuning: master equation / closed form")
for kap, mk, col, ls in zip(kappas, ('o', 's', '^', 'D'), sequential_colors(4, 'plasma', 0.0, 0.8), ('-', '--', '-.', ':')):
    Fm = fidelity_model(nbF, x_line, kap, phiF)
    ax.plot(x_line, Fm, color=col, ls=ls, lw=0.9)
    F_me = [heralded_with_loss(nbF, xx, kap, phiF, width=8.0)[1] for xx in x_pts]
    ax.plot(x_pts, F_me, mk, ms=3.4, mfc='none', mec=col, mew=0.9)
    ax.text(2.0, min(Fm[-1], F_me[-1]) - 0.012, (r'$\kappa/g = 0$' if kap == 0 else rf'${kap:g}$'), fontsize=7,
            ha='right', va='top', color=col)
    if 0 < kap and np.argmax(Fm) > 0:
        ax.plot(x_line[np.argmax(Fm)], Fm.max(), '|', color=col, ms=9, mew=1.2)
    print(f"   kappa/g = {kap:5.3f}: " + "  ".join(f"x={xx:.2f}: {f:.3f}/{fidelity_model(nbF, xx, kap, phiF):.3f}"
                                                 for xx, f in zip(x_pts, F_me))
          + f"   best x (closed form) = {x_line[np.argmax(Fm)]:.2f}")
G0 = loss_exponent(nbF, 0.0, 1.0, phiF)          # resonant exponent per unit kappa
kk = np.linspace(1e-5, 0.05, 20000)
k_th = kk[np.argmax(G0 * kk / (np.exp(G0 * kk) + 1) >= 2 * phiF ** 2 / (4 + phiF ** 2))]     # first crossing
print(f"   detuning stops helping at kappa/g = {k_th:.4f}   (small-angle estimate 3 g/(2 D n_bar) = {1.5 / (D_CAT * nbF):.4f})")
ax.set_xlabel(r'$\Delta/2g\sqrt{\bar{n}}$'); ax.set_ylabel(r'heralded-cat fidelity $F$')
ax.set_xlim(-0.08, 2.08); ax.set_ylim(0.6, 1.02)
panel_label(ax, 'f', rf'$\bar{{n}} = {nbF}$', fontsize=8)

save(fig, PAPER_FIG / 'fig_detuned_shear')
