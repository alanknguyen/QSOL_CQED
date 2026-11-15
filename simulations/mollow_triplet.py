# mollow_triplet.py
# Resonance fluorescence spectrum (Mollow triplet) of a driven two-level atom,
# computed with the quantum regression theorem in QuTiP.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585, taught by Prof. Luca Dal Negro
#
# Generates (paper/figures/): mollow_triplet_driving_strength
#
# Physics: a two-level atom driven on resonance (Delta = 0) fluoresces at
# three frequencies, omega_L and omega_L +/- Omega.  The central peak has
# HWHM = gamma/2 and the sidebands HWHM = 3 gamma/4; in the strong-driving
# limit the peak heights are 3:1 (centre:sideband) and the areas 2:1.  These
# follow from the eigenvalues of the optical Bloch equations,
#   lambda_0 = -gamma/2,  lambda_{+/-} = -3 gamma/4 +/- i sqrt(Omega^2 - gamma^2/16).

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from qutip import sigmam, sigmap, sigmaz, correlation_2op_1t, steadystate, expect

from plotstyle import use_style, save, sequential_colors, HALF_W, PAPER_FIG

use_style()


def mollow_spectrum(nu, Omega, gamma, Delta=0.0, tau_max=100.0, n_tau=20000):
    """Incoherent fluorescence spectrum S(nu), nu = omega - omega_L, from
    <sigma+(tau) sigma-(0)>_ss:  S(nu) = (1/pi) Re int_0^inf G(tau) e^{-i nu tau} dtau."""
    H = 0.5 * Delta * sigmaz() + 0.5 * Omega * (sigmap() + sigmam())   # frame rotating at omega_L
    c_ops = [np.sqrt(gamma) * sigmam()]
    tau = np.linspace(0, tau_max, n_tau)
    corr = correlation_2op_1t(H, None, tau, c_ops, sigmap(), sigmam())
    rho_ss = steadystate(H, c_ops)
    corr_inc = corr - expect(sigmap(), rho_ss) * expect(sigmam(), rho_ss)   # drop the elastic part
    return np.array([trapezoid(corr_inc * np.exp(-1j * w * tau), tau).real for w in nu]) / np.pi


gamma = 1.0
Omega_values = [0.5, 2.0, 5.0, 10.0, 15.0]
nu = np.linspace(-20, 20, 2001)
colors = sequential_colors(len(Omega_values), 'viridis', 0.0, 0.8)
OFFSET = 1.15                       # vertical spacing of the stacked, peak-normalised spectra

fig, ax = plt.subplots(figsize=(HALF_W, 2.9))      # one column of the two-column layout
for i, (Om, color) in enumerate(zip(Omega_values, colors)):
    S = mollow_spectrum(nu, Om, gamma)
    assert S.min() > -1e-6 * S.max(), "spectrum should be non-negative"
    S_norm = S / S.max()
    base = i * OFFSET
    ax.plot(nu, S_norm + base, color=color, lw=1.1)
    ax.axhline(base, color='0.8', lw=0.4, zorder=0)
    ax.text(19.3, base + 0.42, rf'$\Omega = {Om:g}\,\gamma$', ha='right', va='bottom', fontsize=8)
    if Om > gamma:                  # expected sideband positions
        split = np.sqrt(Om ** 2 - gamma ** 2 / 16)
        for sgn in (-1, 1):
            ax.plot([sgn * split] * 2, [base - 0.10, base - 0.02], color='k', lw=0.7)
    centre = S_norm[len(nu) // 2]
    side = S_norm[np.argmin(np.abs(nu - Om))]
    print(f"Omega/gamma = {Om:5.1f}: sideband/centre height = {side / centre:.3f}")

ax.set_xlabel(r'$(\omega - \omega_L)/\gamma$')
ax.set_ylabel('spectrum (normalized, offset)')
ax.set_xlim(-20, 20)
ax.set_ylim(-0.2, len(Omega_values) * OFFSET + 0.05)
ax.set_yticks([])
save(fig, PAPER_FIG / 'mollow_triplet_driving_strength')
