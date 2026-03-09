# Quantum States of Light in Cavity QED

**Wigner Function Dynamics and Atom–Field Entanglement in the Jaynes–Cummings Model**

[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![QuTiP](https://img.shields.io/badge/QuTiP-5.2-green.svg)](https://qutip.org/)

> Nguyen Khoi Nguyen (Alan), Boston University
> Prepared under the guidance of Prof. Luca Dal Negro

---

## Overview

This repository contains the simulation code, figures, and manuscript source for a computational study of quantum light–matter interaction in the Jaynes–Cummings model. The paper combines a pedagogical treatment of quantum optical states with original numerical results on:

1. **Wigner function evolution** — time-resolved phase-space snapshots showing the dynamical formation of a Schrödinger cat state during the collapse of Rabi oscillations
2. **Entanglement dynamics** — systematic comparison of the von Neumann entropy across coherent, thermal, squeezed, and Fock initial field states
3. **Decoherence** — quantitative study of how cavity decay destroys non-classical features, with Wigner negativity and purity as metrics

All simulations use [QuTiP](https://qutip.org/) (Quantum Toolbox in Python).

## Repository Structure

```
.
├── paper/
│   ├── merged_paper.tex          # Full manuscript (REVTeX 4.2)
│   ├── figures/                   # All publication figures (PDF + PNG)
│   │   ├── fig_wigner_evolution.pdf
│   │   ├── fig_cat_state_detail.pdf
│   │   ├── fig_wigner_decoherence.pdf
│   │   ├── fig_entanglement_comparison.pdf
│   │   ├── fig_coherent_entropy_purity.pdf
│   │   ├── fig_entropy_nbar_scaling.pdf
│   │   ├── fig_dissipative_entanglement.pdf
│   │   ├── fig_inversion_snapshots.pdf
│   │   └── ...
│   └── supplementary/             # OBE derivations, additional figures
│
├── simulations/
│   ├── sim_wigner_evolution.py    # Wigner function time evolution + decoherence
│   ├── sim_entanglement_dynamics.py  # Von Neumann entropy + purity + scaling
│   └── utils.py                   # Shared helpers (operators, parameters)
│
├── notebooks/
│   └── full_analysis.ipynb        # Interactive walkthrough of all results
│
├── requirements.txt
├── LICENSE
└── README.md
```

## Quick Start

### Requirements

- Python ≥ 3.10
- QuTiP ≥ 5.0
- NumPy, SciPy, Matplotlib

### Installation

```bash
git clone https://github.com/YOUR_USERNAME/jc-quantum-light.git
cd jc-quantum-light
pip install -r requirements.txt
```

### Reproduce all figures

```bash
cd simulations
python sim_wigner_evolution.py          # Figures 1–4 (Wigner dynamics)
python sim_entanglement_dynamics.py     # Figures 5–8 (entanglement)
```

Output figures are saved to the working directory as both PDF (for LaTeX) and PNG (for preview).

### Interactive exploration

```bash
jupyter notebook notebooks/full_analysis.ipynb
```

The notebook walks through each simulation with inline commentary, parameter sweeps, and additional visualizations not included in the paper.

## Key Results

| Figure | Description | Script |
|--------|-------------|--------|
| Fig. 1 | Atomic inversion with snapshot markers | `sim_wigner_evolution.py` |
| Fig. 2 | Wigner function evolution (6 panels) | `sim_wigner_evolution.py` |
| Fig. 3 | Cat state detail + W(x,0) cross-section | `sim_wigner_evolution.py` |
| Fig. 4 | Wigner decoherence vs. κ/g | `sim_wigner_evolution.py` |
| Fig. 5 | Entanglement entropy: 4 field states | `sim_entanglement_dynamics.py` |
| Fig. 6 | Coherent state: inversion + entropy + purity | `sim_entanglement_dynamics.py` |
| Fig. 7 | Entropy scaling with n̄ | `sim_entanglement_dynamics.py` |
| Fig. 8 | Dissipative entanglement dynamics | `sim_entanglement_dynamics.py` |

### Highlight: Cat-state formation at t = t_r/2

<p align="center">
  <img src="paper/figures/fig_cat_state_detail.png" width="700">
</p>

The Wigner function at half the revival time shows two coherent-state lobes connected by interference fringes with Wigner negativity volume δ = 0.43, a definitive signature of non-classicality.

## Simulation Parameters

All simulations use dimensionless units with vacuum Rabi coupling g = 1.

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Mean photon number | n̄ | 10 (default) |
| Fock space truncation | N_cav | 35–50 |
| Cavity decay rate | κ/g | 0–0.2 |
| Spontaneous emission | γ/g | 0 |
| Wigner grid | — | 200 × 200 over [−7, 7]² |

## Citation

If you use this code or build on these results, please cite:

```bibtex
@article{nguyen2025quantum,
  author  = {Nguyen, Nguyen Khoi},
  title   = {Quantum States of Light in Cavity {QED}: A Computational Study
             of {Wigner} Function Dynamics and Atom--Field Entanglement
             in the {Jaynes--Cummings} Model},
  journal = {arXiv preprint arXiv:XXXX.XXXXX},
  year    = {2025},
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

This work was prepared under the guidance of Prof. Luca Dal Negro at Boston University (EC 585). Simulations were performed using [QuTiP](https://qutip.org/), developed by J. R. Johansson, P. D. Nation, and F. Nori.
