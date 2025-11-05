# Flavoured Leptogenesis Simulator

# Project Description
This project simulates the evolution of matter-antimatter asymmetry in the early Universe, focusing on:
- CP violation effects
- Flavour-dependent effects
- Baryon asymmetry generation via solving Boltzmann equations for heavy neutrinos and lepton flavours

# Requirements
- Python >= 3.8
- Libraries: `numpy`, `scipy`, `matplotlib`, `sympy`

Install dependencies:
```bash
pip install numpy scipy matplotlib sympy
# Usage

1-Open the script leptogenesis.py.

2-Select the mode
3-Adjust parameters (M_i, Yukawas, CP asymmetries) as needed.
4-Run the script
5-Outputs:
Flavour-specific lepton asymmetries Y_ΔL
Final baryon asymmetry Y_B
Plots of the evolution of densities
# Notes

Default values are illustrative; check physical realism.

Washout terms and ΔL=2 processes are approximated.

Symbolic mode outputs the equations without solving numerically.
