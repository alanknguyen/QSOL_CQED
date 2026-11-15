#!/usr/bin/env python3
"""
Animation for "Rotation and shear in the Jaynes-Cummings model" (paper/ejp/jc_beyond_linear).

A larger field does not make a better cat state.  One of the two field components,
|Phi_-> = sum_n c_n exp(-i g t sqrt(n+1)) |n>, is drawn at t = t_r/2 in its own frame
(tangential coordinate horizontal, radial coordinate vertical, vacuum units) while the
mean photon number n_bar grows from 4 to 400.
  - left: its exact Wigner function.  The crescent of a small field straightens, but it
    does not shrink to a coherent state (dashed circle): it converges to a FIXED sheared
    Gaussian, covariance V(s) with s = pi/2 (solid contour), because the shear
    accumulated by t_r/2 is the same in vacuum units for every n_bar.
  - right, top: squeezing of the component, 10 log10(2 lambda_min) with lambda_min the
    smallest eigenvalue of its covariance matrix; limit 6.26 dB, reached slowly.
  - right, bottom: fidelity of the half-revival field with the best cat state
    |beta> + e^{i theta}|-beta> (the quantity of table II of the paper); limit
    (1 + pi^2/16)^(-1/2) = 0.786, reached quickly and never exceeded.

The Wigner function is obtained from the wavefunction, W(x,p) = (1/pi) int psi*(x+y)
psi(x-y) e^{2ipy} dy, with the component turned onto the p axis so that psi(x) stays
within |x| < 12 for every n_bar; this is checked against QuTiP for n_bar = 10.

Generates: animations/anim_cat_limit.gif
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import use_style, wigner_image, fig_to_frame, save_gif, COLORS, ANIM_DIR

use_style()
g = 1.0
PI = np.pi
S_HALF = PI / 2                                           # shear at t_r/2
F_LIMIT = (1 + PI ** 2 / 16) ** -0.5
DB_LIMIT = -10 * np.log10(2 * np.linalg.eigvalsh(0.5 * np.array([[1, S_HALF], [S_HALF, 1 + S_HALF ** 2]]))[0])

H = 0.025                                                 # spacing of the wavefunction grid
X = np.arange(-480, 481) * H                              # |x| <= 12
TAU = X[240:721:3]                                        # tangential window, |tau| <= 6, step 0.075
RHO = np.linspace(-3.6, 3.6, 97)                          # radial window
J = np.arange(-240, 241)                                  # y = J H, |y| <= 6


def fock(nb):
    n = np.arange(int(nb + 8 * np.sqrt(nb) + 20))
    return n, np.sqrt(poisson.pmf(n, nb))


def hermite_functions(nmax, x):
    """Oscillator eigenfunctions phi_n(x), n = 0 .. nmax-1, by the stable three-term recurrence."""
    out = np.empty((nmax, len(x)))
    out[0] = PI ** -0.25 * np.exp(-x * x / 2)
    out[1] = np.sqrt(2.0) * x * out[0]
    for m in range(1, nmax - 1):
        out[m + 1] = np.sqrt(2.0 / (m + 1)) * x * out[m] - np.sqrt(m / (m + 1.0)) * out[m - 1]
    return out


def component(nb):
    """Fock amplitudes of |Phi_-> at t_r/2, turned so that its centre lies on the +p axis; returns
    (n, amplitudes, distance of the centre from the origin, covariance matrix in (radial, tangential))."""
    n, c = fock(nb)
    psi = c * np.exp(-1j * PI * np.sqrt(nb) * np.sqrt(n + 1.0))
    a1 = np.sum(np.conj(psi[:-1]) * psi[1:] * np.sqrt(n[1:]))
    psi = psi * np.exp(-1j * np.angle(a1) * n)            # centre on the +x axis: x radial, p tangential
    a1 = np.sum(np.conj(psi[:-1]) * psi[1:] * np.sqrt(n[1:]))
    a2 = np.sum(np.conj(psi[:-2]) * psi[2:] * np.sqrt(n[2:] * (n[2:] - 1.0)))
    dn, d2 = np.sum(n * abs(psi) ** 2) - abs(a1) ** 2, a2 - a1 * a1
    V = np.array([[0.5 + dn + d2.real, d2.imag], [d2.imag, 0.5 + dn - d2.real]])
    return n, psi * np.exp(0.5j * PI * n), np.sqrt(2) * abs(a1), V


def wigner_local(n, amp, p0):
    """W on the (TAU, RHO) window: tangential tau = -x, radial rho = p - p0."""
    psi = amp @ hermite_functions(len(n), X)
    assert abs(np.sum(abs(psi) ** 2) * H - 1) < 1e-6 and abs(psi[0]) + abs(psi[-1]) < 1e-8
    kernel = np.exp(2j * np.outer(p0 + RHO, J * H))       # (rho, y)
    W = np.empty((len(RHO), len(TAU)))
    for k, tau in enumerate(TAU):
        i = int(round(-tau / H)) + 480                     # x = -tau
        W[:, k] = (kernel @ (np.conj(psi[i + J]) * psi[i - J])).real * H / PI
    return W


def squeezing_db(V):
    return -10 * np.log10(2 * np.linalg.eigvalsh(V)[0])


def cat_fidelity(nb):
    """max over beta, theta of <cat|rho_field|cat> at t_r/2, |cat> ~ |beta> + e^{i theta}|-beta>."""
    n, c = fock(nb)
    th = PI * np.sqrt(nb) * np.sqrt(n + 1.0)
    ce = c * np.cos(th)
    cg = np.concatenate([[0.0], (-1j * c * np.sin(th))[:-1]])
    lg, thetas = 0.5 * gammaln(n + 1), np.linspace(0, 2 * PI, 2001)

    def neg_f(z):
        beta = complex(z[0], z[1])
        kp = np.exp(-abs(beta) ** 2 / 2 + n * np.log(beta) - lg)
        km = kp * (-1.0) ** n
        num = 0.0
        for v in (ce, cg):
            a, b = np.vdot(kp, v), np.vdot(km, v)
            num = num + abs(a) ** 2 + abs(b) ** 2 + 2 * np.real(np.exp(1j * thetas) * np.conj(a) * b)
        return -np.max(num / (2 + 2 * np.real(np.exp(1j * thetas) * np.exp(-2 * abs(beta) ** 2))))

    z = np.array([0.0, np.sqrt(nb)])
    for side in (0.3, 0.03):                               # explicit simplex: a start coordinate is zero
        z = minimize(neg_f, z, method='Nelder-Mead', options={'initial_simplex': np.vstack([z, z + side * np.eye(2)]),
                                                              'xatol': 1e-8, 'fatol': 1e-13}).x
    return -neg_f(z)


# the wavefunction route must reproduce QuTiP's Wigner function (n_bar = 10)
_n, _amp, _p0, _V = component(10)
_ref = qt.wigner(qt.Qobj(_amp), -TAU[::8], _p0 + RHO[::8])     # rows: p, columns: x = -tau
assert np.abs(wigner_local(_n, _amp, _p0)[::8, ::8] - _ref).max() < 1e-6

nb_frames = np.geomspace(4, 400, 41)
nb_line = np.geomspace(4, 400, 61)
print("squeezing and cat fidelity against n_bar ...")
db_line = np.array([squeezing_db(component(x)[3]) for x in nb_line])
F_line = np.array([cat_fidelity(x) for x in nb_line])
for x in (10, 100, 400):
    print(f"   n_bar = {x:3d}: squeezing {squeezing_db(component(x)[3]):.2f} dB (limit {DB_LIMIT:.2f}), "
          f"F_cat = {cat_fidelity(x):.4f} (limit {F_LIMIT:.4f})")

tau_c = np.linspace(0, 2 * PI, 200)
V_lim = 0.5 * np.array([[1 + S_HALF ** 2, S_HALF], [S_HALF, 1.0]])           # (tangential, radial) order, Phi_-
ell = np.sqrt(2) * np.linalg.cholesky(V_lim) @ np.vstack([np.cos(tau_c), np.sin(tau_c)])

frames, durations = [], []
for nb in nb_frames:
    n, amp, p0, V = component(nb)
    W = wigner_local(n, amp, p0)
    fig = plt.figure(figsize=(8.4, 3.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1])
    ax = fig.add_subplot(gs[:, 0])
    im = wigner_image(ax, W, TAU, pvec=RHO)
    ax.plot(ell[0], ell[1], color='k', lw=0.9, label=r'sheared Gaussian, $s = \pi/2$')
    ax.plot(np.cos(tau_c), np.sin(tau_c), color='k', lw=0.8, ls='--', label='coherent state')
    ax.set_xlim(TAU[0], TAU[-1]); ax.set_ylim(RHO[0], RHO[-1])
    ax.set_xticks([-4, -2, 0, 2, 4]); ax.set_yticks([-2, 0, 2])
    ax.set_xlabel('tangential coordinate'); ax.set_ylabel('radial coordinate')
    ax.legend(loc='upper center', fontsize=7.5, ncol=2, handlelength=1.8, columnspacing=1.5)
    fig.colorbar(im, ax=ax, shrink=0.75, aspect=20, pad=0.03, ticks=[-0.3, 0, 0.3], label=r'$W$')
    ax.set_title(rf'one field component at $t_r/2$, $\bar{{n}} = {nb:.0f}$', fontsize=9)

    for row, (y, y_now, lim, lab, col) in enumerate((
            (db_line, squeezing_db(V), DB_LIMIT, 'squeezing (dB)', COLORS['blue']),
            (F_line, cat_fidelity(nb), F_LIMIT, r'cat fidelity $F_{\rm cat}$', COLORS['orange']))):
        a2 = fig.add_subplot(gs[row, 1])
        a2.axhline(lim, color='0.35', lw=0.7, ls='--')
        a2.semilogx(nb_line, y, color=col, lw=1.1)
        a2.plot(nb, y_now, 'o', color=col, ms=5, mec='k', mew=0.5, zorder=5)
        a2.set_xlim(3.5, 460); a2.set_ylabel(lab)
        if row == 0:
            a2.set_ylim(0, 7.2); a2.tick_params(labelbottom=False)
            a2.text(0.97, 0.10, rf'{y_now:.1f} dB (limit {lim:.1f} dB)', transform=a2.transAxes, ha='right', fontsize=8)
            a2.set_title('a larger field does not make a better cat', fontsize=9)
        else:
            a2.set_ylim(0.70, 0.80); a2.set_xlabel(r'mean photon number $\bar{n}$')
            a2.text(0.97, 0.10, rf'{y_now:.3f} (limit {lim:.3f})', transform=a2.transAxes, ha='right', fontsize=8)
    frames.append(fig_to_frame(fig, dpi=110))
    durations.append(230)
durations[0] = 1500
durations[-1] = 2500
save_gif(frames, ANIM_DIR / 'anim_cat_limit.gif', durations)
