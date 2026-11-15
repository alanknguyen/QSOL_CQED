#!/usr/bin/env python3
"""
Animation for "Rotation and shear in the Jaynes-Cummings model" (paper/ejp/jc_beyond_linear).

The exact Wigner function of the reduced cavity field for n_bar = 10 and n_bar = 40, at
equal rotation angle phi = g t / (2 sqrt(n_bar)), with the prediction of the
rotation-and-shear picture drawn on top and no adjustable parameter: each black curve is
the 1/e contour of a Gaussian centred at distance sqrt(2 n_bar) from the origin, at angle
-/+ phi' with phi' = g t / (2 sqrt(n_bar + 1)), sheared by s = phi (covariance matrix
V(+/-s) of the paper).  The third panel is the overlap |<Phi_-|Phi_+>| of the two
components, i.e. the envelope of the atomic inversion: it collapses when they separate
and returns only to (1 + pi^2)^(-1/4) = 0.55 at phi = pi, where they coincide in
position but are sheared in opposite senses.  The animation pauses at phi = pi/2 (the
cat state, shear pi/2, 6.3 dB of squeezing in each component) and at phi = pi.

Everything is computed from the exact amplitudes; no differential equation is solved.

Generates: animations/anim_rotation_shear.gif
"""

import sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plotstyle import (use_style, wigner_image, fig_to_frame, save_gif, COLORS, ANIM_DIR,
                       GUIDE_HALF, GUIDE_FULL)

g = 1.0
PI = np.pi
CASES = {10: dict(L=7.0, npts=281), 40: dict(L=12.6, npts=365)}       # window half-width, grid points
PHI = np.arange(0, 116) * PI / 100                                     # rotation angle; pi/2 and pi are frames 50, 100


def amplitudes(nb, t):
    """Exact amplitudes of |e,n> and |g,n> for |e> x |alpha = sqrt(nb)>."""
    n = np.arange(int(nb + 7 * np.sqrt(nb) + 14))
    c = np.sqrt(poisson.pmf(n, nb))
    th = g * t * np.sqrt(n + 1.0)
    cg = np.zeros(len(n), complex)
    cg[1:] = (-1j * c * np.sin(th))[:-1]
    return n, c * np.cos(th), cg


def wigner_frame(args):
    nb, phi = args
    _, ce, cg = amplitudes(nb, 2 * np.sqrt(nb) * phi / g)
    rho = qt.Qobj(np.outer(ce, ce.conj()) + np.outer(cg, cg.conj()))
    xv = np.linspace(-CASES[nb]['L'], CASES[nb]['L'], CASES[nb]['npts'])
    return qt.wigner(rho, xv, xv)


def overlap(nb, phi):
    """|<Phi_-|Phi_+>| = |sum_n p_n exp(2 i g t sqrt(n+1))|: the envelope of the inversion."""
    n = np.arange(int(nb + 9 * np.sqrt(nb) + 25))
    t = 2 * np.sqrt(nb) * np.atleast_1d(phi) / g
    return np.abs((poisson.pmf(n, nb)[:, None] * np.exp(2j * g * np.sqrt(n + 1.0)[:, None] * t[None, :])).sum(0))


def model_ellipses(nb, phi, level=2.0):
    """1/e contours of the two Gaussian components: angle -/+ phi', covariance V(+/-s) with s = phi."""
    s = phi
    phi_c = phi * np.sqrt(nb / (nb + 1.0))                 # the n+1 of the exact phases, kept in the position
    tau = np.linspace(0, 2 * PI, 200)
    out = []
    for sign in (+1, -1):
        ang = -sign * phi_c
        V = 0.5 * np.array([[1.0, sign * s], [sign * s, 1 + s * s]])           # (radial, tangential)
        R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        centre = np.sqrt(2 * nb) * np.array([np.cos(ang), np.sin(ang)])
        out.append(np.sqrt(level) * np.linalg.cholesky(R @ V @ R.T) @ np.vstack([np.cos(tau), np.sin(tau)]) + centre[:, None])
    return out


if __name__ == '__main__':
    use_style()
    print(f"Wigner functions: {len(PHI)} frames for n_bar = 10 and 40 ...")
    with Pool() as pool:
        W = {nb: pool.map(wigner_frame, [(nb, p) for p in PHI], chunksize=4) for nb in CASES}

    phi_line = np.linspace(0, PHI[-1], 1500)
    env = {nb: overlap(nb, phi_line) for nb in CASES}
    here = {nb: overlap(nb, PHI) for nb in CASES}
    colour = {10: COLORS['blue'], 40: COLORS['green']}
    print(f"overlap at phi = pi: {here[10][100]:.4f} (n_bar = 10), {here[40][100]:.4f} (n_bar = 40); "
          f"closed form {(1 + PI ** 2) ** -0.25:.4f}")

    frames, durations = [], []
    for k, phi in enumerate(PHI):
        fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.25), gridspec_kw={'width_ratios': [1, 1, 1.12]})
        for ax, nb in zip(axes[:2], CASES):
            xv = np.linspace(-CASES[nb]['L'], CASES[nb]['L'], CASES[nb]['npts'])
            wigner_image(ax, W[nb][k], xv)
            for pts, ls in zip(model_ellipses(nb, phi), ('-', '--')):
                ax.plot(pts[0], pts[1], color='k', lw=0.8, ls=ls)
            ax.set_xlim(xv[0], xv[-1]); ax.set_ylim(xv[0], xv[-1])
            ticks = [-5, 0, 5] if nb == 10 else [-10, 0, 10]
            ax.set_xticks(ticks); ax.set_yticks(ticks)
            ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$p$', labelpad=1)
            ax.set_title(rf'$\bar{{n}} = {nb}$: exact $W$, model contours', fontsize=9)
        ax = axes[2]
        for nb in CASES:
            ax.plot(phi_line / PI, env[nb], color='0.8', lw=0.6)
            sel = phi_line <= phi
            ax.plot(phi_line[sel] / PI, env[nb][sel], color=colour[nb], lw=1.0, label=rf'$\bar{{n}} = {nb}$')
            ax.plot(phi / PI, here[nb][k], 'o', color=colour[nb], ms=4.5, mec='k', mew=0.5, zorder=5)
        ax.axhline((1 + PI ** 2) ** -0.25, color='0.4', lw=0.6, ls=':')
        ax.text(0.53, (1 + PI ** 2) ** -0.25 + 0.025, r'$(1+\pi^2)^{-1/4}$', fontsize=7.5, color='0.3')
        ax.axvline(0.5, **GUIDE_HALF); ax.axvline(1.0, **GUIDE_FULL)
        ax.set_xlim(0, PHI[-1] / PI); ax.set_ylim(-0.03, 1.05)
        ax.set_xlabel(r'rotation angle $\varphi/\pi = gt/2\pi\sqrt{\bar{n}}$')
        ax.set_ylabel(r'overlap $|\langle\Phi_-|\Phi_+\rangle|$')
        ax.legend(loc='upper center', fontsize=7.5, ncol=2, handlelength=1.4, columnspacing=1.2)
        ax.set_title(rf'$\varphi = {phi / PI:.2f}\,\pi$, shear $s = \varphi = {phi:.2f}$', fontsize=9)
        frames.append(fig_to_frame(fig, dpi=100))
        durations.append({0: 1000, 50: 2000, 100: 2000}.get(k, 90))
    durations[-1] = 1200
    save_gif(frames, ANIM_DIR / 'anim_rotation_shear.gif', durations)
