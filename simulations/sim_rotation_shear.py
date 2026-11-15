# sim_rotation_shear.py
# Beyond the linear approximation in the resonant Jaynes-Cummings model.
#
# Expanding the Rabi frequencies 2 g sqrt(n+1) about the mean photon number,
#   order 1 (slope):      the field splits into two components that ROTATE by
#                         -/+ phi(t),  phi = g t / (2 sqrt(n_bar)): collapse and revival;
#   order 2 (curvature):  the angular velocity g/(2 sqrt(n)) depends on n, so each
#                         component is SHEARED by s = phi (differential rotation);
#   next order (1/n_bar): explicit corrections, e.g. the residual atom-field
#                         entanglement at t_r/2.
# Every Gaussian-overlap result below follows from one integral,
#   < exp(-i a x^2) >_{x ~ N(0, n_bar)} = (1 + 2 i a n_bar)^(-1/2),
# and is compared with the EXACT state, which for |e> x |alpha> is known in closed form:
#   |psi(t)> = sum_n c_n [ cos(g t sqrt(n+1)) |e,n> - i sin(g t sqrt(n+1)) |g,n+1> ].
#
# Generates (paper/figures/): fig_revival_envelopes, fig_wigner_shear,
#     fig_portrait_model, fig_convergence, fig_loss_scaling
# and prints the table of closed forms against exact values.
#
# Nguyen Khoi Nguyen (Alan), Boston University

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import poisson

from plotstyle import (use_style, save, panel_label, wigner_image, time_guides, sequential_colors,
                       GUIDE_HALF, GUIDE_FULL,
                       FULL_W, COLORS, PAPER_FIG, FIG_DIR)
from utils import wigner_negativity

use_style()
g = 1.0
PI = np.pi


# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------
def fock_cutoff(nb):
    """Poisson weight beyond this is < 1e-30 for nb >= 1."""
    return int(nb + 14 * np.sqrt(nb) + 40)


def amplitudes(nb, t):
    """c_n and the exact amplitudes of |e,n> and |g,n> at time t."""
    N = fock_cutoff(nb)
    n = np.arange(N)
    c = np.sqrt(poisson.pmf(n, nb))
    th = g * t * np.sqrt(n + 1)
    ce = c * np.cos(th)
    cg = np.zeros(N, complex)
    cg[1:] = -1j * c[:-1] * np.sin(th[:-1])
    return n, c, ce, cg


def purity_entropy(nb, t):
    """Field purity and atomic entropy (bits) from the rank-2 reduced state."""
    _, _, ce, cg = amplitudes(nb, t)
    pe, pg, ov = np.vdot(ce, ce).real, np.vdot(cg, cg).real, np.vdot(ce, cg)
    P = pe ** 2 + pg ** 2 + 2 * abs(ov) ** 2
    lam = 0.5 * (1 - np.sqrt(max(2 * P - 1, 0.0)))
    return P, binary_entropy(lam)


def binary_entropy(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-16)
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


def inversion_and_parity(nb, t):
    """<sigma_z>(t) = sum p_n cos(2 g t sqrt(n+1));  parity = sum p_n (-1)^n cos(2 g t sqrt(n+1))."""
    N = fock_cutoff(nb)
    n = np.arange(N)
    p = poisson.pmf(n, nb)
    ph = 2 * g * np.sqrt(n + 1)[:, None] * np.atleast_1d(t)[None, :]
    return (p[:, None] * np.cos(ph)).sum(0), ((p * (-1.0) ** n)[:, None] * np.cos(ph)).sum(0)


def Z_exact(k, nb):
    """Z_k = sum_n p_n exp(2 i g t_k sqrt(n+1) - 2 pi i k n) at t_k = 2 pi k sqrt(nb)/g.
    Re Z_k is the inversion at the k-th revival (k integer) or the photon-number
    parity at t = k t_r (k half-integer)."""
    N = fock_cutoff(nb)
    n = np.arange(N)
    return np.sum(poisson.pmf(n, nb) * np.exp(1j * (4 * PI * k * np.sqrt(nb * (n + 1)) - 2 * PI * k * n)))


def cat_fidelity_exact(nb, t):
    """max over beta, theta of <cat|rho_field|cat>, |cat> ~ |beta> + e^{i theta}|-beta>."""
    n, _, ce, cg = amplitudes(nb, t)
    lg = 0.5 * gammaln(n + 1)
    thetas = np.linspace(0, 2 * PI, 2001)

    def neg_f(z):
        beta = z[0] + 1j * z[1]
        kp = np.exp(-abs(beta) ** 2 / 2 + n * np.log(beta) - lg)
        km = kp * (-1.0) ** n
        num = 0.0
        for v in (ce, cg):
            a, b = np.vdot(kp, v), np.vdot(km, v)
            num = num + abs(a) ** 2 + abs(b) ** 2 + 2 * np.real(np.exp(1j * thetas) * np.conj(a) * b)
        den = 2 + 2 * np.real(np.exp(1j * thetas) * np.exp(-2 * abs(beta) ** 2))
        return -np.max(num / den)

    phi = g * t / (2 * np.sqrt(nb + 1))
    b0 = np.sqrt(nb) * np.exp(1j * phi)
    res = minimize(neg_f, [b0.real, b0.imag], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-12, 'maxiter': 4000})
    return -res.fun


def reduced_field(nb, t, N=None):
    n, _, ce, cg = amplitudes(nb, t)
    N = N or int(nb + 9 * np.sqrt(nb) + 25)
    return qt.Qobj(np.outer(ce[:N], ce[:N].conj()) + np.outer(cg[:N], cg[:N].conj()))


def hermite_functions(nmax, x):
    """Oscillator eigenfunctions phi_n(x), n = 0 .. nmax-1, by the stable three-term recurrence."""
    out = np.empty((nmax, len(x)))
    out[0] = PI ** -0.25 * np.exp(-x * x / 2)
    if nmax > 1:
        out[1] = np.sqrt(2.0) * x * out[0]
    for m in range(1, nmax - 1):
        out[m + 1] = np.sqrt(2.0 / (m + 1)) * x * out[m] - np.sqrt(m / (m + 1.0)) * out[m - 1]
    return out


def negativity_half_revival(nb, pts_per_fringe=31.7, X=16.0):
    """Wigner negativity delta = int |W| dx dp - 1 of the reduced field at t = pi sqrt(nb)/g.

    The reduced state has rank 2, rho = |A><A| + |B><B|, so W is the sum of two
    pure-state Wigner functions, each a Fourier transform in y of
    psi*(x + y) psi(x - y).  At t_r/2 the two lobes sit on the p axis, so psi(x)
    is confined to |x| < X for any nb and only the grid spacing has to follow the
    fringe period pi / sqrt(2 nb).  A non-integer number of points per fringe
    avoids locking the samples of |cos| to the fringes.  Cost ~ nb, against
    ~ nb^4 for a Laguerre evaluation on a square grid.
    """
    n, _, ce, cg = amplitudes(nb, PI * np.sqrt(nb) / g)
    dx = PI / np.sqrt(2 * nb) / pts_per_fringe
    M = int(2 * np.ceil(X / dx))
    x = (np.arange(M) - M // 2) * dx
    Hm = hermite_functions(len(n), x)
    psis = [ce @ Hm, cg @ Hm]
    J = np.arange(-M // 2, M // 2)
    dp = PI / (M * dx)
    total = 0.0
    for i in range(M):
        ip, im = i + J, i - J
        ok = (ip >= 0) & (ip < M) & (im >= 0) & (im < M)
        row = np.zeros(M, complex)
        for psi in psis:
            row[ok] += np.conj(psi[ip[ok]]) * psi[im[ok]]
        total += np.abs(np.fft.fft(np.fft.ifftshift(row)).real).sum()
    return total * dx / PI * dx * dp - 1


# ---------------------------------------------------------------------------
# Asymptotic results
# ---------------------------------------------------------------------------
def G(k):
    """Leading amplitude of Z_k: overlap of two Gaussians sheared by +s and -s, s = pi k."""
    return (1 + 1j * PI * k) ** -0.5


def j2(k):
    """Coefficient of 1/n_bar in Z_k = exp(2 pi i k (n_bar+1)) G(k) [1 + j2(k)/n_bar + ...]."""
    K = PI * k
    return -K * (8 * K ** 2 + 3) / (96 * (K - 1j) ** 3)


def envelope(t, nb, j):
    """|S_j(t)|: envelope of the j-th Poisson-summation term (inversion: j integer; parity: j half-integer)."""
    L1 = g * t / np.sqrt(nb + 1)
    q = g * t * nb / (2 * (nb + 1) ** 1.5)             # = -Lambda_2 n_bar: the shear of the overlap
    return (1 + q ** 2) ** -0.25 * np.exp(-nb * (L1 - 2 * PI * j) ** 2 / (2 * (1 + q ** 2)))


def purity_model(t, nb, order=1):
    """Field purity in the collapsed interval; order 0 is the n_bar -> infinity limit."""
    phi = g * t / (2 * np.sqrt(nb))
    P = 1 - 0.5 * np.cos(phi) ** 2
    if order:
        P = P - ((1 + phi ** 2) * np.sin(phi) ** 2 + phi * np.sin(phi) * np.cos(phi)) / (8 * nb)
    return P


F_INF = (1 + PI ** 2 / 16) ** -0.5                    # cat fidelity
A_PAR = (1 + PI ** 2 / 4) ** -0.25                    # parity amplitude at t_r/2
CHI = 0.5 * np.arctan(PI / 2)                         # and its phase
D_INF = (2 / PI) * (1 + PI ** 2 / 4) ** 0.25          # Wigner negativity
C_PUR = (4 + PI ** 2) / 32                            # n_bar (1 - purity)
SQ_DB = 10 * np.log10((1 + PI ** 2 / 8) - np.sqrt((1 + PI ** 2 / 8) ** 2 - 1))   # squeezing of one component


# ---------------------------------------------------------------------------
# Table: closed forms against exact values
# ---------------------------------------------------------------------------
print("Closed forms (n_bar -> infinity) against the exact state")
print(f"  revival height |G(1)| = {abs(G(1)):.4f}; parity amplitude |G(1/2)| = {A_PAR:.4f}, phase {CHI:.4f}; "
      f"cat fidelity |G(1/4)|^2 = {F_INF:.4f}; delta = {D_INF:.4f}; n(1-P) = {C_PUR:.4f}; squeezing {SQ_DB:.2f} dB")
print("   n_bar   |Z_1|   +1/n    |Z_1/2|  +1/n   parity(t_r/2)  formula   F_cat   n(1-P)   S_min    h(C/2n)")
for nb in (10, 25, 100, 400):
    t_half = PI * np.sqrt(nb) / g
    z1, zh = Z_exact(1, nb), Z_exact(0.5, nb)
    P, S = purity_entropy(nb, t_half)
    par_model = np.real(np.exp(1j * PI * (nb + 1)) * G(0.5) * (1 + j2(0.5) / nb))
    print(f"  {nb:5d}  {abs(z1):.4f}  {abs(G(1) * (1 + j2(1) / nb)):.4f}  {abs(zh):.4f}  "
          f"{abs(G(0.5) * (1 + j2(0.5) / nb)):.4f}    {zh.real:+.4f}     {par_model:+.4f}   "
          f"{cat_fidelity_exact(nb, t_half):.4f}  {nb * (1 - P):.4f}   {S:.4f}   {float(binary_entropy(C_PUR / (2 * nb))):.4f}")


# ---------------------------------------------------------------------------
# Fig A: inversion and parity with their analytic envelopes
# ---------------------------------------------------------------------------
nb = 10
t_r = 2 * PI * np.sqrt(nb) / g
t = np.linspace(0, 2.3 * t_r, 9000)
W, Par = inversion_and_parity(nb, t)

figA, (axa, axb) = plt.subplots(2, 1, figsize=(FULL_W, 3.9), sharex=True)
# every term of the Poisson sum that reaches into the plotted window (2.3 t_r) is included
for ax, y, js, color, ylabel in ((axa, W, (0, 1, 2, 3), COLORS['blue'], r'$\langle \sigma_z \rangle$'),
                                 (axb, Par, (0.5, 1.5, 2.5), COLORS['green'], r'$\langle (-1)^{\hat n} \rangle$')):
    ax.plot(t * g, y, color=color, lw=0.6)
    ax.axhline(0, color='0.6', lw=0.4)
    env = sum(envelope(t, nb, j) for j in js)
    ax.plot(t * g, env, color='k', lw=0.8, ls='--')
    ax.plot(t * g, -env, color='k', lw=0.8, ls='--')
    ax.set_ylabel(ylabel)
    ax.set_ylim(-1.08, 1.08)
    ax.set_yticks([-1, 0, 1])
    ax.set_xlim(0, t[-1] * g)
# The closed forms |G(k)| are the amplitudes at t = k t_r (the dotted/dashed guides); the bump
# itself peaks slightly later, at 2 pi k sqrt(n_bar + 1)/g, and a few per cent higher.
for k in (1, 2):
    axa.plot(k * t_r, abs(G(k)), 'o', ms=4, mfc='white', mec='k', mew=0.8, zorder=5)
    axa.text(k * t_r + (1.2 if k == 1 else -0.6), abs(G(k)) + 0.17, rf'$(1+{"" if k == 1 else "4"}\pi^2)^{{-1/4}}$',
             fontsize=8, ha='left' if k == 1 else 'right')
for k in (0.5, 1.5):
    axb.plot(k * t_r, abs(G(k)), 'o', ms=4, mfc='white', mec='k', mew=0.8, zorder=5)
axb.text(0.5 * t_r + 1.6, abs(G(0.5)) + 0.12, r'$(1+\pi^2/4)^{-1/4}$', fontsize=8, ha='left')
axb.text(1.5 * t_r + 1.6, abs(G(1.5)) + 0.22, r'$(1+9\pi^2/4)^{-1/4}$', fontsize=8, ha='left')
trans = axa.get_xaxis_transform()
for m, lab in ((0.5, r'$t_r/2$'), (1, r'$t_r$'), (1.5, r'$3t_r/2$'), (2, r'$2t_r$')):
    for ax in (axa, axb):          # same guide styles as in the other figures
        ax.axvline(m * t_r * g, **(GUIDE_FULL if float(m).is_integer() else GUIDE_HALF))
    axa.text(m * t_r * g, 1.03, lab, transform=trans, ha='center', va='bottom', fontsize=8)
axb.set_xlabel(r'$gt$')
panel_label(axa, 'a')
panel_label(axb, 'b')
save(figA, PAPER_FIG / 'fig_revival_envelopes')


# ---------------------------------------------------------------------------
# Fig B: Wigner function with the rotation-and-shear ellipses
# ---------------------------------------------------------------------------
def model_ellipses(nb, t, level=2.0):
    """1/e contours of the two Gaussian components: centre sqrt(2 nb) at angle -/+ phi, shear +/- s."""
    phi = g * t / (2 * np.sqrt(nb + 1))
    s = phi
    tau = np.linspace(0, 2 * PI, 200)
    out = []
    for sign in (+1, -1):                              # +1: component at angle -phi
        ang = -sign * phi
        V_loc = 0.5 * np.array([[1.0, sign * s], [sign * s, 1 + s * s]])     # (radial, tangential)
        R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        Lc = np.linalg.cholesky(R @ V_loc @ R.T)
        pts = np.sqrt(level) * Lc @ np.vstack([np.cos(tau), np.sin(tau)])
        centre = np.sqrt(2 * nb) * np.array([np.cos(ang), np.sin(ang)])
        out.append(pts + centre[:, None])
    return out


fracs = [(1 / 8, r'$t_r/8$'), (1 / 4, r'$t_r/4$'), (1 / 2, r'$t_r/2$'), (1, r'$t_r$')]
figB, axesB = plt.subplots(2, 4, figsize=(FULL_W, 3.45))
letters = iter('abcdefgh')
for row, nbB in enumerate((10, 40)):
    lim = np.sqrt(2 * nbB) + 3.6
    xv = np.linspace(-lim, lim, int(2 * lim / (PI / np.sqrt(2 * nbB) / 10)) | 1)
    for col, (fr, lab) in enumerate(fracs):
        ax = axesB[row, col]
        tt = fr * 2 * PI * np.sqrt(nbB) / g
        im = wigner_image(ax, qt.wigner(reduced_field(nbB, tt), xv, xv), xv)
        for pts, ls in zip(model_ellipses(nbB, tt), ('-', '--')):
            ax.plot(pts[0], pts[1], color='k', lw=0.7, ls=ls)
        panel_label(ax, next(letters), lab + (rf', $\bar{{n}} = {nbB}$' if col == 0 else ''), fontsize=8)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        tk = [-5, 0, 5] if nbB == 10 else [-10, 0, 10]
        ax.set_xticks(tk); ax.set_yticks(tk)
        ax.tick_params(labelleft=(col == 0))
        ax.set_xlabel(r'$x$')
    axesB[row, 0].set_ylabel(r'$p$')
cb = figB.colorbar(im, ax=axesB, shrink=0.85, aspect=30, pad=0.02, ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb.set_label(r'$W(x,p)$')
save(figB, PAPER_FIG / 'fig_wigner_shear')


# ---------------------------------------------------------------------------
# Fig C: inversion, entropy and purity with the analytic curves
# ---------------------------------------------------------------------------
tC = np.linspace(0, 1.35 * t_r, 2700)
WC, _ = inversion_and_parity(nb, tC)
PS = np.array([purity_entropy(nb, x) for x in tC])
env01 = envelope(tC, nb, 0) + envelope(tC, nb, 1)
valid = env01 < 0.02                                   # between collapse and revival
P0 = np.where(valid, purity_model(tC, nb, 0), np.nan)
P1 = np.where(valid, purity_model(tC, nb, 1), np.nan)
ent = lambda P: binary_entropy(0.5 * (1 - np.sqrt(np.clip(2 * P - 1, 0, 1))))

figC, (c1, c2, c3) = plt.subplots(3, 1, figsize=(FULL_W, 4.9), sharex=True)
c1.plot(tC * g, WC, color=COLORS['blue'], lw=0.7)
c1.plot(tC * g, env01, color='k', lw=0.8, ls='--', label='analytic envelope')
c1.plot(tC * g, -env01, color='k', lw=0.8, ls='--')
c1.axhline(0, color='0.6', lw=0.4)
c1.set_ylabel(r'$\langle \sigma_z \rangle$'); c1.set_ylim(-1.08, 1.08); c1.set_yticks([-1, 0, 1])
c1.legend(loc='lower left', bbox_to_anchor=(0.40, 0.0))      # between the two time guides
c2.plot(tC * g, PS[:, 1], color=COLORS['red'], lw=1.0, label='exact')
c2.plot(tC * g, ent(P0), color='k', lw=0.8, ls=':', label=r'$\bar{n} \to \infty$')
c2.plot(tC * g, ent(P1), color='k', lw=0.9, ls='--', label=r'with $1/\bar{n}$ term')
c2.axhline(1, color='0.6', lw=0.4)
c2.set_ylabel(r'$S(\rho_{\rm atom})$ (bit)'); c2.set_ylim(-0.05, 1.1); c2.set_yticks([0, 0.5, 1])
c2.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), borderaxespad=0.0, ncol=3, columnspacing=1.0,
          handlelength=1.8)                                         # above the axes, clear of the guides
c3.plot(tC * g, PS[:, 0], color=COLORS['green'], lw=1.0)
c3.plot(tC * g, P0, color='k', lw=0.8, ls=':')
c3.plot(tC * g, P1, color='k', lw=0.9, ls='--')
c3.axhline(0.5, color='0.6', lw=0.4)
c3.set_ylabel(r'$\mathcal{P} = \mathrm{Tr}\,\rho_{\rm field}^2$'); c3.set_ylim(0.42, 1.04)
c3.set_yticks([0.5, 0.75, 1.0]); c3.set_xlabel(r'$gt$')
for ax, letter in zip((c1, c2, c3), 'abc'):
    time_guides(ax, t_r * g, labels=(letter == 'a'))
    panel_label(ax, letter)
    ax.set_xlim(0, tC[-1] * g)
save(figC, PAPER_FIG / 'fig_portrait_model')


# ---------------------------------------------------------------------------
# Fig D: convergence of the asymptotic constants, and where the expansion fails
# ---------------------------------------------------------------------------
nbs = np.array([2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 64, 80, 100, 150, 200, 300, 400])
rows = []
for x in nbs:
    th = PI * np.sqrt(x) / g
    P, _ = purity_entropy(x, th)
    rows.append((abs(Z_exact(1, x)), abs(Z_exact(0.5, x)), cat_fidelity_exact(x, th), x * (1 - P)))
rows = np.array(rows)
nb_d = nbs[nbs >= 4]
delta_d = np.array([negativity_half_revival(x) for x in nb_d])
print(f"delta at n_bar = 1600, 6400: {negativity_half_revival(1600, 12.37):.5f}, {negativity_half_revival(6400, 12.37):.5f}")
print("delta at t_r/2: " + ", ".join(f"n={a}: {b:.4f}" for a, b in zip(nb_d, delta_d)) + f"   (limit {D_INF:.4f})")
assert abs(negativity_half_revival(10) - wigner_negativity(
    qt.wigner(reduced_field(10, PI * np.sqrt(10) / g), np.linspace(-7, 7, 401), np.linspace(-7, 7, 401)),
    np.linspace(-7, 7, 401))) < 5e-4, "FFT and Laguerre evaluations of delta disagree"

series = [  # label, exact values, leading constant, colour, marker
    (r'revival height $|Z_1|$', rows[:, 0], abs(G(1)), COLORS['blue'], 'o'),
    (r'parity amplitude $|Z_{1/2}|$', rows[:, 1], A_PAR, COLORS['green'], 's'),
    (r'cat fidelity $F_{\rm cat}$', rows[:, 2], F_INF, COLORS['orange'], 'D'),
    (r'impurity $\bar{n}(1-\mathcal{P})$', rows[:, 3], C_PUR, COLORS['red'], '^'),
]
figD, (d1, d2) = plt.subplots(1, 2, figsize=(FULL_W, 2.75))
for lab, y, c0, col, mk in series:
    d1.semilogx(nbs, y / c0, marker=mk, ms=3, lw=0.7, color=col, label=lab)
    d2.loglog(nbs, np.abs(y / c0 - 1), marker=mk, ms=3, lw=0.7, color=col)
d1.semilogx(nb_d, delta_d / D_INF, marker='v', ms=3, lw=0.7, color=COLORS['purple'], label=r'negativity $\delta$')
d2.loglog(nb_d, np.abs(delta_d / D_INF - 1), marker='v', ms=3, lw=0.7, color=COLORS['purple'])
for k, col, mk in ((1, COLORS['blue'], 'o'), (0.5, COLORS['green'], 's')):    # with the 1/n term
    z = np.array([abs(Z_exact(k, x)) for x in nbs])
    d2.loglog(nbs, np.abs(z / np.abs(G(k) * (1 + j2(k) / nbs)) - 1), marker=mk, ms=3, mfc='white', lw=0.7,
              ls='--', color=col)
xx = np.array([3.0, 400.0])
d2.loglog(xx, 0.5 / xx, color='0.5', lw=0.6); d2.text(120, 0.012, r'$\propto \bar{n}^{-1}$', fontsize=7.5, color='0.3')
d2.loglog(xx, 0.5 / xx ** 2, color='0.5', lw=0.6); d2.text(130, 4.5e-6, r'$\propto \bar{n}^{-2}$', fontsize=7.5, color='0.3')
d1.axhline(1, color='0.6', lw=0.4)
d1.set_ylim(0.45, 1.2)
d1.set_xlabel(r'$\bar{n}$'); d1.set_ylabel('exact / asymptotic')
d1.legend(fontsize=7, loc='lower right', labelspacing=0.25, borderaxespad=0.3)
d2.set_xlabel(r'$\bar{n}$'); d2.set_ylabel('relative error')
d2.set_ylim(2e-9, 1)          # keeps the corrected series in view down to n_bar = 400
panel_label(d1, 'a'); panel_label(d2, 'b')
save(figD, PAPER_FIG / 'fig_convergence')


# ---------------------------------------------------------------------------
# Fig E: cavity loss; needs figures/phase_diagram_data.npz (sim_phase_diagram.py)
# ---------------------------------------------------------------------------
data_file = FIG_DIR / 'phase_diagram_data.npz'
if data_file.exists():
    d = np.load(data_file)
    nbar_grid, kappa_grid, delta_map, kc = d['nbar'], d['kappa'], d['delta'], d['kappa_c_delta']
    pos = kappa_grid > 0
    figE, (e1, e2) = plt.subplots(1, 2, figsize=(FULL_W, 2.6))
    x = np.linspace(0, 6, 200)
    e1.semilogy(x, np.exp(-x), color='k', lw=0.9, ls='--', label=r'$e^{-\pi\kappa\bar{n}^{3/2}/g}$')
    ns = [4, 6, 10, 16, 25]
    for n_, c_, mk in zip(ns, sequential_colors(len(ns), 'viridis', 0.0, 0.85), 'osD^v'):
        j = int(np.argmin(np.abs(nbar_grid - n_)))
        e1.semilogy(PI * kappa_grid[pos] * n_ ** 1.5, delta_map[pos, j] / delta_map[0, j], marker=mk, ms=3,
                    lw=0.6, color=c_, label=rf'$\bar{{n}} = {n_}$')
    e1.set_xlim(0, 6); e1.set_ylim(3e-3, 1.3)
    e1.set_xlabel(r'$\pi\kappa\bar{n}^{3/2}/g$'); e1.set_ylabel(r'$\delta(\kappa)/\delta(0)$')
    e1.legend(fontsize=7, loc='upper right', labelspacing=0.2, handlelength=1.6)
    e2.loglog(nbar_grid, kc, 'o', ms=3.5, color='k', label=r'simulation, $\delta = 0.05$')
    e2.loglog(nbar_grid[2:], np.log(delta_map[0, 2:] / 0.05) / (PI * nbar_grid[2:] ** 1.5), color='k', lw=0.9,
              ls='--', label=r'$\ln(\delta_0/0.05)\,/\,\pi\bar{n}^{3/2}$')
    e2.set_xlim(1.5, 30)
    e2.set_xticks([2, 5, 10, 20]); e2.set_xticklabels(['2', '5', '10', '20'])
    e2.xaxis.set_minor_formatter(plt.NullFormatter())
    e2.set_xlabel(r'$\bar{n}$'); e2.set_ylabel(r'$\kappa_c/g$')
    e2.legend(fontsize=7, loc='lower left')
    panel_label(e1, 'a'); panel_label(e2, 'b')
    save(figE, PAPER_FIG / 'fig_loss_scaling')
else:
    print("fig_loss_scaling skipped: run sim_phase_diagram.py first")
