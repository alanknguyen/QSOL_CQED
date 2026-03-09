"""
Shared utilities for Jaynes-Cummings simulations.

Provides operator constructors, default parameters, and helper functions
used across all simulation scripts.
"""

import numpy as np
from qutip import (
    basis, tensor, destroy, sigmam, sigmap, sigmaz,
    coherent, thermal_dm, squeeze, ket2dm, Qobj
)


# ── Default physical parameters ──────────────────────────────────

DEFAULT_PARAMS = {
    "g": 1.0,           # Vacuum Rabi coupling (sets time/energy unit)
    "n_bar": 10.0,       # Mean photon number
    "N_cav": 35,         # Fock space truncation
    "kappa": 0.0,        # Cavity decay rate
    "gamma": 0.0,        # Spontaneous emission rate
}


# ── Operator constructors ────────────────────────────────────────

def jc_operators(N_cav):
    """
    Construct Jaynes-Cummings operators in the joint
    atom ⊗ cavity Hilbert space.

    Parameters
    ----------
    N_cav : int
        Fock space truncation for the cavity.

    Returns
    -------
    dict with keys: 'a', 'sm', 'sp', 'sz', 'H_int'
        a     : cavity annihilation operator
        sm    : atomic lowering operator σ⁻
        sp    : atomic raising operator σ⁺
        sz    : Pauli-Z operator σ_z
        H_int : interaction Hamiltonian a†σ⁻ + a σ⁺ (without g prefactor)
    """
    I_cav = Qobj(np.eye(N_cav))
    I_atom = Qobj(np.eye(2))

    a  = tensor(destroy(N_cav), I_atom)
    sm = tensor(I_cav, sigmam())
    sp = tensor(I_cav, sigmap())
    sz = tensor(I_cav, sigmaz())

    H_int = a.dag() * sm + a * sp

    return {"a": a, "sm": sm, "sp": sp, "sz": sz, "H_int": H_int}


def jc_hamiltonian(N_cav, g=1.0, omega_c=0.0, omega_a=0.0):
    """
    Full Jaynes-Cummings Hamiltonian in the rotating frame.

    H = ω_c a†a + (ω_a/2) σ_z + g (a†σ⁻ + a σ⁺)
    """
    ops = jc_operators(N_cav)
    H = (omega_c * ops["a"].dag() * ops["a"]
         + 0.5 * omega_a * ops["sz"]
         + g * ops["H_int"])
    return H, ops


# ── Initial state constructors ───────────────────────────────────

def initial_state(N_cav, field_type="coherent", n_bar=10.0, atom="excited"):
    """
    Construct the initial atom ⊗ field state.

    Parameters
    ----------
    N_cav : int
        Fock space truncation.
    field_type : str
        One of 'coherent', 'thermal', 'squeezed', 'fock'.
    n_bar : float
        Mean photon number.
    atom : str
        'excited' or 'ground'.

    Returns
    -------
    rho0 : Qobj
        Initial density matrix in the joint Hilbert space.
    """
    # Atom state
    if atom == "excited":
        rho_atom = ket2dm(basis(2, 0))
    else:
        rho_atom = ket2dm(basis(2, 1))

    # Field state
    if field_type == "coherent":
        alpha = np.sqrt(n_bar)
        rho_field = ket2dm(coherent(N_cav, alpha))
    elif field_type == "thermal":
        rho_field = thermal_dm(N_cav, n_bar)
    elif field_type == "squeezed":
        r = np.arcsinh(np.sqrt(n_bar))
        psi = squeeze(N_cav, r) * basis(N_cav, 0)
        rho_field = ket2dm(psi)
    elif field_type == "fock":
        n = int(round(n_bar))
        rho_field = ket2dm(basis(N_cav, n))
    else:
        raise ValueError(f"Unknown field_type: {field_type}")

    return tensor(rho_field, rho_atom)


# ── Characteristic timescales ────────────────────────────────────

def collapse_time(g=1.0):
    """Approximate collapse time t_c ≈ 1/g."""
    return 1.0 / g


def revival_time(n_bar, g=1.0):
    """First revival time t_r = 2π√n̄ / g."""
    return 2.0 * np.pi * np.sqrt(n_bar) / g
