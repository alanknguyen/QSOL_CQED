# utils.py
# Nguyen Khoi Nguyen (Alan), Boston University
# Began as a course project taught by Prof. Luca Dal Negro

import numpy as np
from qutip import (
    basis, tensor, destroy, sigmam, sigmap, sigmaz,
    coherent, thermal_dm, squeeze, ket2dm, Qobj
)

DEFAULT_PARAMS = {
    "g": 1.0,        # vacuum Rabi coupling [sets time and energy units]
    "n_bar": 10.0,    # mean photon number
    "N_cav": 35,      # Fock space truncation
    "kappa": 0.0,     # cavity decay rate
    "gamma": 0.0,     # spontaneous emission rate
}


def jc_operators(N_cav):
    """Joint atom x cavity operators for the Jaynes-Cummings model."""
    I_cav = Qobj(np.eye(N_cav))
    I_atom = Qobj(np.eye(2))

    a  = tensor(destroy(N_cav), I_atom)       # cavity annihilation
    sm = tensor(I_cav, sigmam())               # sigma^-
    sp = tensor(I_cav, sigmap())               # sigma^+
    sz = tensor(I_cav, sigmaz())               # sigma_z

    H_int = a.dag() * sm + a * sp              # a^dag sigma^- + a sigma^+

    return {"a": a, "sm": sm, "sp": sp, "sz": sz, "H_int": H_int}


def jc_hamiltonian(N_cav, g=1.0, omega_c=0.0, omega_a=0.0):
    """H_JC = omega_c a^dag a + (omega_a/2) sigma_z + g (a^dag sigma^- + a sigma^+)"""
    ops = jc_operators(N_cav)
    H = (omega_c * ops["a"].dag() * ops["a"]
         + 0.5 * omega_a * ops["sz"]
         + g * ops["H_int"])
    return H, ops


def initial_state(N_cav, field_type="coherent", n_bar=10.0, atom="excited"):
    """
    Initial density matrix rho_0 = rho_field x rho_atom.
    field_type: 'coherent', 'thermal', 'squeezed', 'fock'
    """
    rho_atom = ket2dm(basis(2, 0 if atom == "excited" else 1))

    if field_type == "coherent":
        rho_field = ket2dm(coherent(N_cav, np.sqrt(n_bar)))
    elif field_type == "thermal":
        rho_field = thermal_dm(N_cav, n_bar)
    elif field_type == "squeezed":
        # <n> = sinh^2(r)
        r = np.arcsinh(np.sqrt(n_bar))
        rho_field = ket2dm(squeeze(N_cav, r) * basis(N_cav, 0))
    elif field_type == "fock":
        rho_field = ket2dm(basis(N_cav, int(round(n_bar))))
    else:
        raise ValueError(f"Unknown field_type: {field_type}")

    return tensor(rho_field, rho_atom)


def collapse_time(g=1.0):
    """Operational collapse time t_c = 1/g (the dephasing estimate gives O(1/g))."""
    return 1.0 / g


def wigner_negativity(W, xvec):
    """Wigner negativity volume delta = int |W| dx dp - 1 (Kenfack & Zyczkowski 2004).

    Zero for any state with non-negative W; equals twice the integrated
    negative volume, since int W = 1.
    """
    dx = xvec[1] - xvec[0]
    return float(np.sum(np.abs(W)) * dx**2 - 1.0)


def revival_time(n_bar, g=1.0):
    """t_r = 2 pi sqrt(n_bar) / g"""
    return 2.0 * np.pi * np.sqrt(n_bar) / g


def squeezed_vacuum(N_cav, n_bar):
    """Squeezed vacuum with <n> = sinh^2(r) = n_bar, built from its exact Fock
    amplitudes  c_{2m} = (-tanh r)^m sqrt((2m)!) / (2^m m! sqrt(cosh r)).

    This avoids exponentiating a truncated generator.  The photon-number tail is
    heavy, p(2m) ~ (n_bar/(n_bar+1))^m / sqrt(pi m): for n_bar = 10 the weight
    beyond n = 150 is 1.6e-4 and N_cav = 300 is needed to push it below 1e-6.
    """
    from scipy.special import gammaln
    r = np.arcsinh(np.sqrt(n_bar))
    m = np.arange(0, (N_cav + 1) // 2)
    log_amp = (0.5 * gammaln(2 * m + 1) - m * np.log(2.0) - gammaln(m + 1)
               + m * np.log(np.tanh(r)) - 0.5 * np.log(np.cosh(r)))
    c = np.zeros(N_cav, dtype=complex)
    c[2 * m] = (-1.0) ** m * np.exp(log_amp)
    return Qobj(c.reshape(-1, 1)).unit()


def cat_fidelity(rho_field, n_bar, return_params=False):
    """Fidelity of a field state with the best-fitting two-component cat,

        F = max_{beta, theta} <cat| rho |cat>,   |cat> ~ |beta> + e^{i theta} |-beta>,

    with beta scanned around i sqrt(n_bar) (the half-revival cat of an initially
    excited atom and real alpha sits on the p axis) and theta over [0, 2 pi).
    With return_params=True also returns the optimal (beta, theta).
    """
    N = rho_field.shape[0]
    thetas = np.linspace(0, 2 * np.pi, 721)
    best = (0.0, None, None)
    for scale in np.linspace(0.85, 1.10, 11):
        for dphi in np.linspace(-0.2, 0.2, 9):
            beta = scale * np.sqrt(n_bar) * np.exp(1j * (np.pi / 2 + dphi))
            kp = coherent(N, beta, method='analytic').full().ravel()
            km = coherent(N, -beta, method='analytic').full().ravel()
            R = rho_field.full()
            rpp = np.real(np.vdot(kp, R @ kp))
            rmm = np.real(np.vdot(km, R @ km))
            c = np.vdot(kp, R @ km)            # <beta| rho |-beta>
            ov = np.vdot(kp, km)               # <beta|-beta> = exp(-2|beta|^2)
            F = ((rpp + rmm + 2 * np.real(np.exp(1j * thetas) * c))
                 / (2 + 2 * np.real(np.exp(1j * thetas) * ov)))
            k = int(np.argmax(F))
            if F[k] > best[0]:
                best = (float(F[k]), beta, float(thetas[k]))
    return best if return_params else best[0]


def photon_parity(rho_field):
    """<(-1)^n> of a field state."""
    p = np.real(rho_field.diag())
    return float(np.sum(p[0::2]) - np.sum(p[1::2]))