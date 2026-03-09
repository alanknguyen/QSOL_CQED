# utils.py
# Nguyen Khoi Nguyen (Alan), Boston University
# Advised by Prof. Luca Dal Negro

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
    """t_c ~ 1/g"""
    return 1.0 / g


def revival_time(n_bar, g=1.0):
    """t_r = 2 pi sqrt(n_bar) / g"""
    return 2.0 * np.pi * np.sqrt(n_bar) / g