# mollow_triplet.py
# Resonance fluorescence spectrum (Mollow triplet) of a driven two-level atom.
# Computed via the quantum regression theorem using QuTiP.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: mollow_triplet_driving_strength.png
#
# Physics: a two-level atom driven on resonance (Delta = 0) fluoresces at
# three frequencies: omega_L and omega_L +/- Omega. The central peak has
# HWHM = gamma/2; the sidebands have HWHM = 3*gamma/4. In the strong-driving
# limit the peak height ratio is 3:1 (center:side). These widths and positions
# follow from the eigenvalues of the optical Bloch equation matrix:
#   lambda_0 = -gamma/2,  lambda_{+/-} = -3*gamma/4 +/- i*Omega.

import numpy as np
import matplotlib.pyplot as plt
from qutip import (basis, sigmam, sigmap, sigmaz,
                    correlation_2op_1t, steadystate, expect)


def mollow_spectrum(nu_array, Omega, gamma, Delta=0.0,
                    tau_max=100, n_tau=20000):
    """
    Incoherent fluorescence spectrum S(nu) where nu = omega - omega_L.
    Computed from <sigma+(tau) sigma-(0)>_ss via the quantum regression theorem.
    """
    # Hamiltonian in the rotating frame at omega_L
    H = (Delta / 2) * sigmaz() + (Omega / 2) * (sigmap() + sigmam())
    c_ops = [np.sqrt(gamma) * sigmam()]

    # Two-time correlation <sigma+(tau) sigma-(0)>_ss
    taulist = np.linspace(0, tau_max, n_tau)
    corr = correlation_2op_1t(H, None, taulist, c_ops, sigmap(), sigmam())

    # Subtract elastic (coherent) component
    rho_ss = steadystate(H, c_ops)
    corr_inc = corr - expect(sigmap(), rho_ss) * expect(sigmam(), rho_ss)

    # Fourier transform: S(nu) = (1/pi) Re int_0^inf G(tau) exp(-i nu tau) dtau
    dtau = taulist[1] - taulist[0]
    S = np.zeros(len(nu_array))
    for k, nu in enumerate(nu_array):
        S[k] = (1 / np.pi) * np.real(
            np.trapezoid(corr_inc * np.exp(-1j * nu * taulist), dx=dtau))

    return np.maximum(S, 0)


# -- Spectrum vs driving strength (Fig for paper) --
gamma = 1.0
Omega_values = [0.5, 2.0, 5.0, 10.0, 15.0]
nu = np.linspace(-20, 20, 2000)
colors = plt.cm.viridis(np.linspace(0, 1, len(Omega_values)))

fig, ax = plt.subplots(figsize=(12, 7))
for i, Om in enumerate(Omega_values):
    print(f"Computing: Omega/gamma = {Om}")
    S = mollow_spectrum(nu, Om, gamma)
    S_norm = S / np.max(S) if np.max(S) > 0 else S
    offset = i * 0.3
    ax.plot(nu, S_norm + offset, color=colors[i], linewidth=2.5,
            label=rf'$\Omega = {Om}\gamma$')
    if Om > gamma:
        ax.axvline(Om, color=colors[i], ls='--', alpha=0.3)
        ax.axvline(-Om, color=colors[i], ls='--', alpha=0.3)

ax.set_xlabel(r'$\omega - \omega_L$ [units of $\gamma$]', fontsize=14)
ax.set_ylabel('Normalized emission intensity (offset)', fontsize=14)
ax.set_title('Mollow Triplet: Dependence on Driving Strength', fontsize=15)
ax.legend(fontsize=12)
ax.set_xlim(-20, 20)
fig.tight_layout()
fig.savefig('mollow_triplet_driving_strength.png', dpi=300, bbox_inches='tight')
fig.savefig('mollow_triplet_driving_strength.pdf', bbox_inches='tight')
print("Saved: mollow_triplet_driving_strength")
