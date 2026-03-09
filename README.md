# Quantum States of Light in Cavity QED

**Wigner Function Dynamics and Atom-Field Entanglement in the Jaynes-Cummings Model**

<p align="center">
  <img src="paper/figures/banner_qsol.png" width="100%">
</p>

[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![QuTiP 5.2](https://img.shields.io/badge/QuTiP-5.2-green.svg)](https://qutip.org/)

Nguyen Khoi Nguyen (Alan), Boston University  
Advised by Prof. Luca Dal Negro, EC 585 / EC 777

---

## Overview

Computational study of quantum light-matter interaction in the Jaynes-Cummings (JC) model. All dynamics are computed from first principles using QuTiP, with no phenomenological approximations.

**Three principal results:**

1. Time-resolved Wigner function snapshots showing Schrodinger cat-state formation at `t = t_r/2` with Wigner negativity `delta = 0.43`
2. Systematic comparison of atom-field entanglement entropy across coherent, thermal, squeezed, and Fock initial field states
3. Quantitative decoherence study: cavity decay `kappa/g = 0.02` reduces cat-state Wigner negativity by over 80%

---

## Animations

### Wigner Function Evolution

A coherent state splits into a Schrodinger cat state during the collapse of Rabi oscillations, then partially re-localizes at the first revival.

<p align="center">
  <img src="animations/anim_wigner_evolution.gif" width="750">
</p>

Left: Wigner function `W(x,p)` of the reduced cavity field. Right: atomic inversion `<sigma_z>(t)` with a moving time marker. The animation pauses at the cat-state time (`t = t_r/2`) and at the first revival (`t = t_r`).

### Entanglement, Inversion, and Purity

Entanglement entropy saturates at 1 bit during collapse (atom and field maximally entangled, cat state forms) and dips at revival (system approximately refactorizes).

<p align="center">
  <img src="animations/anim_entanglement.gif" width="850">
</p>

Left: Wigner function. Center: inversion. Right: entropy (red) and field purity (green). Purity drops to 0.5 during collapse, consistent with a two-component superposition.

### Decoherence Destroys the Cat State

Cavity photon loss erases the interference fringes on a timescale `~ 1/(kappa * n_bar)`, far shorter than the bare cavity lifetime.

<p align="center">
  <img src="animations/anim_decoherence.gif" width="750">
</p>

Left: Wigner function at `t = t_r/2` as `kappa/g` increases from 0 to 0.15. Right: Wigner negativity volume `delta` tracking the loss of non-classicality.

---

## Static Figures

### Cat-State Detail

<p align="center">
  <img src="paper/figures/fig_cat_state_detail.png" width="750">
</p>

Cross-section `W(x, p=0)` confirms deep negativity (`W ~ -0.22`) between the two coherent components. Fringe spacing `~ pi / sqrt(2 n_bar) ~ 0.7`.

### Entanglement Across Field States

<p align="center">
  <img src="paper/figures/fig_entanglement_comparison.png" width="800">
</p>

Only the coherent state (blue) shows a clean entropy dip at `t_r`. Fock oscillates periodically. Thermal dephases permanently. Squeezed is intermediate.

### Entropy Scaling with Photon Number

<p align="center">
  <img src="paper/figures/fig_entropy_nbar_scaling.png" width="800">
</p>

`t_c ~ 1/g` (independent of `n_bar`). `t_r = 2 pi sqrt(n_bar) / g`. Clean revivals emerge for `n_bar >= 9`.

### Dissipative Entanglement

<p align="center">
  <img src="paper/figures/fig_dissipative_entanglement.png" width="800">
</p>

### Decoherence Table

| `kappa/g` | Wigner negativity `delta` | Field purity | Cat state visible? |
|-----------|--------------------------|--------------|-------------------|
| 0.00      | 0.425                    | 0.960        | Yes               |
| 0.02      | 0.071                    | 0.460        | Marginal          |
| 0.05      | 0.011                    | 0.418        | No                |
| 0.10      | < 0.001                  | 0.386        | No                |

### Coherent vs Thermal Inversion

<p align="center">
  <img src="paper/figures/jaynes_cummings_comparison.png" width="800">
</p>

`W(t) = sum_n P(n) cos(2 sqrt(n+1) t)`. Left: Poisson weights (coherent) produce collapse-revival. Right: Bose-Einstein weights (thermal) produce permanent dephasing.

### Mollow Triplet

<p align="center">
  <img src="paper/figures/mollow_triplet_driving_strength.png" width="800">
</p>

Computed via the quantum regression theorem (QuTiP `correlation_2op_1t`). Sidebands at `+/- Omega` with HWHM = `3 gamma/4`, peak ratio 3:1.

### Static Wigner Functions for Fock States

<p align="center">
  <img src="paper/figures/wigner_fock_combined.png" width="800">
</p>

`W_n(x,p) = (-1)^n / pi * L_n(2 r^2) * exp(-r^2)`. Ring structure and negativity grow with `n`.

---

## Repository Structure

```
.
├── paper/
│   ├── merged_paper.tex              # Full manuscript (REVTeX 4.2)
│   └── figures/                      # All static figures (PDF + PNG) + banner
│
├── simulations/
│   ├── utils.py                      # JC operators, initial states, timescales
│   ├── sim_wigner_evolution.py       # Wigner W(x,p) dynamics + decoherence
│   ├── sim_entanglement_dynamics.py  # Von Neumann entropy + purity + scaling
│   ├── jaynes_cummings_comparison.py # Coherent vs thermal inversion
│   ├── mollow_triplet.py             # Fluorescence spectrum via regression theorem
│   └── wigner_fock_states.py         # Static Wigner functions for Fock states
│
├── animations/
│   ├── anim_wigner_evolution.py      # Wigner evolution GIF generator
│   ├── anim_entanglement.py          # Entropy + Wigner GIF generator
│   ├── anim_decoherence.py           # Decoherence sweep GIF generator
│   ├── thumbnail_banner.py           # Dark banner image generator
│   ├── anim_wigner_evolution.gif     # 120 frames, 2.6 MB
│   ├── anim_entanglement.gif         # 120 frames, 2.3 MB
│   └── anim_decoherence.gif          # 60 frames, 750 KB
│
├── requirements.txt
├── LICENSE
└── README.md
```

## Simulation Details

### Jaynes-Cummings Hamiltonian

```
H_JC = omega_c a^dag a + (omega_a / 2) sigma_z + g (a^dag sigma^- + a sigma^+)
```

On resonance: `H = g (a^dag sigma^- + a sigma^+)`.

### Lindblad Master Equation

```
d rho / dt = -i [H, rho] + kappa D[a] rho + gamma D[sigma^-] rho
```

### Wigner Negativity Volume

```
delta = integral |W(x,p)| dx dp - 1
```

`delta = 0` for classical states. `delta > 0` for non-classical states.

### Von Neumann Entropy

```
S(rho_atom) = -Tr[rho_atom log_2 rho_atom]
```

Ranges from 0 (product state) to 1 bit (maximally entangled).

### Mollow Spectrum

Eigenvalues of the Bloch matrix on resonance:

```
lambda_0 = -gamma/2           (central peak, HWHM = gamma/2)
lambda_+/- = -3gamma/4 +/- i Omega   (sidebands, HWHM = 3gamma/4)
```

### Parameters

| Parameter | Symbol | Default |
|-----------|--------|---------|
| Vacuum Rabi coupling | `g` | 1.0 |
| Mean photon number | `n_bar` | 10 |
| Fock truncation | `N_cav` | 35-50 |
| Cavity decay | `kappa/g` | 0-0.2 |
| Collapse time | `t_c` | `~ 1/g` |
| Revival time | `t_r` | `2 pi sqrt(n_bar) / g` |

### Figure-to-Script Map

| Output | Script |
|--------|--------|
| `fig_inversion_snapshots` | `sim_wigner_evolution.py` |
| `fig_wigner_evolution` | `sim_wigner_evolution.py` |
| `fig_cat_state_detail` | `sim_wigner_evolution.py` |
| `fig_wigner_decoherence` | `sim_wigner_evolution.py` |
| `fig_entanglement_comparison` | `sim_entanglement_dynamics.py` |
| `fig_coherent_entropy_purity` | `sim_entanglement_dynamics.py` |
| `fig_entropy_nbar_scaling` | `sim_entanglement_dynamics.py` |
| `fig_dissipative_entanglement` | `sim_entanglement_dynamics.py` |
| `jaynes_cummings_comparison` | `jaynes_cummings_comparison.py` |
| `mollow_triplet_driving_strength` | `mollow_triplet.py` |
| `wigner_fock_combined` | `wigner_fock_states.py` |
| `anim_wigner_evolution.gif` | `anim_wigner_evolution.py` |
| `anim_entanglement.gif` | `anim_entanglement.py` |
| `anim_decoherence.gif` | `anim_decoherence.py` |
| `banner_qsol.png` | `thumbnail_banner.py` |

---

## Quick Start

```bash
git clone https://github.com/alanknguyen/QSOL_CQED.git
cd QSOL_CQED
pip install -r requirements.txt
```

Static figures:

```bash
cd simulations
python sim_wigner_evolution.py            # ~2 min
python sim_entanglement_dynamics.py       # ~5 min
python jaynes_cummings_comparison.py      # ~10 sec
python mollow_triplet.py                  # ~1 min
python wigner_fock_states.py              # ~5 sec
```

Animations:

```bash
cd animations
python anim_wigner_evolution.py           # ~5 min
python anim_entanglement.py              # ~5 min
python anim_decoherence.py              # ~3 min
python thumbnail_banner.py              # ~30 sec
```

---

## Citation

```bibtex
@article{nguyen2025quantum,
  author  = {Nguyen, Nguyen Khoi},
  title   = {Quantum States of Light in Cavity {QED}: A Computational Study
             of {Wigner} Function Dynamics and Atom-Field Entanglement
             in the {Jaynes-Cummings} Model},
  journal = {arXiv preprint arXiv:XXXX.XXXXX},
  year    = {2025},
}
```

## License

MIT. See [LICENSE](LICENSE).

## Acknowledgments

Prepared under the guidance of Prof. Luca Dal Negro at Boston University (EC 585, EC 777). Simulations use [QuTiP](https://qutip.org/) by J. R. Johansson, P. D. Nation, and F. Nori.
