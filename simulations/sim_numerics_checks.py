# sim_numerics_checks.py
# How far can the numbers be trusted?  Four checks, each against the exact
# closed-form state of the resonant Jaynes-Cummings model,
#   |psi(t)> = sum_n c_n [ cos(g t sqrt(n+1)) |e,n> - i sin(g t sqrt(n+1)) |g,n+1> ],
# which needs no integrator and no truncation beyond the Poisson tail.
#
#   (a) Fock-space truncation: error of the entropy and purity at t_r/2 against N_cav,
#       and the photon-number weight left outside N_cav for thermal and squeezed fields;
#   (b) integrator tolerance: QuTiP mesolve against the exact inversion;
#   (c) phase-space grid: error of the Wigner negativity against points per fringe;
#   (d) time sampling: a trace sampled too coarsely shows a spurious envelope (aliasing).
#
# Generates (paper/figures/): fig_numerics_checks; prints the numbers quoted in the
# supplementary material.
#
# Nguyen Khoi Nguyen (Alan), Boston University

import time
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.special import gammaln
from scipy.stats import poisson

from plotstyle import use_style, save, panel_label, FULL_W, COLORS, PAPER_FIG
from utils import wigner_negativity

use_style()
g = 1.0
nb = 10.0
t_r = 2 * np.pi * np.sqrt(nb) / g
t_half = 0.5 * t_r


def exact_state(nb, t, N):
    n = np.arange(N)
    c = np.sqrt(poisson.pmf(n, nb))
    th = g * t * np.sqrt(n + 1)
    ce = c * np.cos(th)
    cg = np.zeros(N, complex)
    cg[1:] = -1j * c[:-1] * np.sin(th[:-1])
    return ce, cg


def purity_entropy(ce, cg):
    pe, pg, ov = np.vdot(ce, ce).real, np.vdot(cg, cg).real, np.vdot(ce, cg)
    norm = pe + pg
    P = (pe ** 2 + pg ** 2 + 2 * abs(ov) ** 2) / norm ** 2
    lam = 0.5 * (1 - np.sqrt(max(2 * P - 1, 0)))
    S = 0.0 if lam <= 0 else -lam * np.log2(lam) - (1 - lam) * np.log2(1 - lam)
    return P, S


def mesolve_state(nb, t, N, **opts):
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.sigmam())
    H = g * (a.dag() * sm + a * sm.dag())
    psi0 = qt.tensor(qt.coherent(N, np.sqrt(nb)), qt.basis(2, 0))
    return H, psi0, qt.tensor(qt.qeye(N), qt.sigmaz())


fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 4.5))
(ax_a, ax_b), (ax_c, ax_d) = axes

# ---------------------------------------------------------------------------
# (a) truncation
# ---------------------------------------------------------------------------
P_ref, S_ref = purity_entropy(*exact_state(nb, t_half, 200))
Ns = np.arange(14, 61, 2)
errS, errP = [], []
for N in Ns:
    # a simulation in a truncated space: coherent state built in N levels, evolved with the truncated H
    H, psi0, _ = mesolve_state(nb, t_half, int(N))
    st = qt.sesolve(H, psi0, [0, t_half], options={'atol': 1e-12, 'rtol': 1e-10, 'nsteps': 100000}).states[-1]
    rho_a = st.ptrace(1)
    errS.append(abs(qt.entropy_vn(rho_a, 2) - S_ref))
    errP.append(abs((st.ptrace(0) ** 2).tr().real - P_ref))
ax_a.semilogy(Ns, errS, 'o-', ms=3, lw=0.8, color=COLORS['red'], label=r'$|\Delta S|$ at $t_r/2$')
ax_a.semilogy(Ns, errP, 's--', ms=3, lw=0.8, color=COLORS['green'], label=r'$|\Delta\mathcal{P}|$ at $t_r/2$')
ax_a.semilogy(Ns, poisson.sf(Ns - 1, nb), color='k', lw=0.8, ls=':', label='Poisson weight beyond $N_{\\rm cav}$')
ax_a.set_xlabel(r'$N_{\rm cav}$'); ax_a.set_ylabel('error')
ax_a.set_ylim(1e-13, 1)
ax_a.legend(fontsize=7, loc='upper right')
panel_label(ax_a, 'a', rf'truncation, coherent field, $\bar{{n}} = {nb:g}$', fontsize=8)
print("(a) truncation, coherent n_bar = 10:")
for N in (20, 25, 30, 35, 50):
    k = int(np.argmin(np.abs(Ns - N)))
    print(f"    N_cav = {Ns[k]:3d}: |dS| = {errS[k]:.1e}, |dP| = {errP[k]:.1e}, Poisson tail = {poisson.sf(Ns[k] - 1, nb):.1e}")
q = nb / (nb + 1)
r = np.arcsinh(np.sqrt(nb))
m = np.arange(0, 3000)
p_sq = np.exp(gammaln(2 * m + 1) - 2 * (m * np.log(2) + gammaln(m + 1)) + 2 * m * np.log(np.tanh(r)) - np.log(np.cosh(r)))
print("    weight beyond N_cav:   thermal (n/(n+1))^N      squeezed vacuum")
for N in (50, 150, 300):
    print(f"      N_cav = {N:3d}:        {q ** N:.1e}                 {p_sq[2 * m >= N].sum():.1e}")

# ---------------------------------------------------------------------------
# (b) integrator tolerance
# ---------------------------------------------------------------------------
N = 50
H, psi0, sz = mesolve_state(nb, 0, N)
tl = np.linspace(0, 1.5 * t_r, 1501)
nn = np.arange(200); pp = poisson.pmf(nn, nb)
W_exact = (pp[:, None] * np.cos(2 * g * np.sqrt(nn + 1)[:, None] * tl[None, :])).sum(0)
print("(b) mesolve (Adams) against the exact inversion, 0 < t < 1.5 t_r, N_cav = 50:")
rtols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
errs, times = [], []
for rtol in rtols:
    t0 = time.time()
    res = qt.mesolve(H, psi0, tl, e_ops=[sz], options={'atol': rtol * 1e-2, 'rtol': rtol, 'nsteps': 200000})
    times.append(time.time() - t0)
    errs.append(np.abs(np.real(res.expect[0]) - W_exact).max())
    print(f"    rtol = {rtol:.0e} (atol = {rtol * 1e-2:.0e}): max |error| = {errs[-1]:.1e}   ({times[-1]:.2f} s)")
res = qt.mesolve(H, psi0, tl, e_ops=[sz])
err_default = np.abs(np.real(res.expect[0]) - W_exact).max()
print(f"    QuTiP defaults (atol 1e-8, rtol 1e-6): max |error| = {err_default:.1e}")
ax_b.loglog(rtols, errs, 'o-', ms=3, lw=0.8, color=COLORS['blue'])
ax_b.axvline(1e-6, color='0.5', lw=0.5, ls=':')
ax_b.text(1.3e-6, min(errs) * 1.5, 'QuTiP default', fontsize=7, color='0.3', rotation=90, va='bottom')
ax_b.set_xlabel('relative tolerance of the integrator'); ax_b.set_ylabel(r'max $|\Delta\langle\sigma_z\rangle|$')
panel_label(ax_b, 'b', 'integrator tolerance', fontsize=8)

# ---------------------------------------------------------------------------
# (c) phase-space grid for the Wigner negativity
# ---------------------------------------------------------------------------
ce, cg = exact_state(nb, t_half, 60)
rho = qt.Qobj(np.outer(ce, ce.conj()) + np.outer(cg, cg.conj()))
period = np.pi / np.sqrt(2 * nb)
L = 7.0


def delta_on_grid(dx):
    x = np.arange(-L, L + dx / 2, dx)
    return wigner_negativity(qt.wigner(rho, x, x), x)


d_ref = delta_on_grid(period / 40)
ppf = np.array([1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 14, 20])
d_err = np.array([abs(delta_on_grid(period / k) - d_ref) for k in ppf])
ax_c.loglog(ppf, d_err, 'o-', ms=3, lw=0.8, color=COLORS['purple'])
ax_c.axvline(2, color='0.5', lw=0.5, ls=':')
ax_c.text(2.08, d_err.min() * 1.3, 'Nyquist', fontsize=7, color='0.3', rotation=90, va='bottom')
_mm = np.array([3.0, 20.0])
ax_c.loglog(_mm, 0.13 * _mm ** -2.0, color='0.5', lw=0.6)            # the M^-2 law of a kinked integrand
ax_c.text(11, 0.13 * 11 ** -2.0 * 2.2, r'$\propto M^{-2}$', fontsize=7.5, color='0.3')
ax_c.set_xticks([2, 3, 5, 10, 20]); ax_c.set_xticklabels(['2', '3', '5', '10', '20'])
ax_c.xaxis.set_minor_formatter(plt.NullFormatter())
ax_c.set_xlabel(r'grid points per fringe, $M$'); ax_c.set_ylabel(r'$|\Delta\delta|$')
panel_label(ax_c, 'c', 'phase-space grid', fontsize=8)
print(f"(c) Wigner negativity at t_r/2 (reference {d_ref:.5f}, fringe period {period:.3f}):")
for k, e in zip(ppf, d_err):
    print(f"    {k:5.1f} points per fringe (dx = {period / k:.3f}): |error| = {e:.1e}")

# ---------------------------------------------------------------------------
# (d) time sampling
# ---------------------------------------------------------------------------
nb_d = 24.0
nn = np.arange(400); pp = poisson.pmf(nn, nb_d)
inv = lambda t: (pp[:, None] * np.cos(2 * g * np.sqrt(nn + 1)[:, None] * np.atleast_1d(t)[None, :])).sum(0)
t_fine = np.linspace(24, 38, 4000)
period_rabi = np.pi / (g * np.sqrt(nb_d + 1))
ax_d.plot(t_fine, inv(t_fine), color='0.35', lw=0.6, label='exact', zorder=2)
# the under-sampled trace is joined by lines (what a plotting routine would show); the adequately
# sampled one is drawn as markers only, so that it does not hide the exact curve
for pts, color, mk, lw_ in ((2.2, COLORS['red'], 'o', 0.9), (6.3, COLORS['blue'], 's', 0.0)):
    ts = np.arange(24, 38, period_rabi / pts)
    ax_d.plot(ts, inv(ts), marker=mk, ms=2.4, lw=lw_, color=color, label=f'{pts:g} samples per Rabi period',
              zorder=3 if lw_ else 4)
ax_d.set_xlabel(r'$gt$'); ax_d.set_ylabel(r'$\langle\sigma_z\rangle$')
ax_d.set_xlim(24, 38); ax_d.set_ylim(-0.75, 1.15)
ax_d.legend(fontsize=7, loc='upper right', ncol=1, labelspacing=0.2)
panel_label(ax_d, 'd', rf'time sampling, $\bar{{n}} = {nb_d:g}$', fontsize=8)
print(f"(d) Rabi period at n_bar = {nb_d:g}: {period_rabi:.3f}/g; the entropy oscillates with half that period")
save(fig, PAPER_FIG / 'fig_numerics_checks')
