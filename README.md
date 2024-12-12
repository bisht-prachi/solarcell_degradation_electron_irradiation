## Solar Cell Degradation Under Electron Irradiation

This repository contains tools to simulate the degradation of GaAs solar cells under 1 MeV electron irradiation. The key file, `iv_kphi_gaas.cpp`, generates the degraded IV characteristics of a GaAs solar cell, capturing the effects of irradiation-induced defect creation.

---

## Overview

The primary degradation mechanism arises from defect creation at deep levels within the solar cell material. These defects act as recombination centers for minority charge carriers, significantly reducing the minority carrier diffusion length L. This reduction is quantified by the diffusion length damage coefficient K_L, which is energy-dependent and determines the impact of radiation on L. The relation governing the degraded diffusion length is:

1/L^2 = 1/Lo^2 + K_L*phi 

Where:  
- Lo: Initial diffusion length (pre-irradiation).  
- K_L: Diffusion length damage coefficient (for electrons and holes).  
- phi: Total electron fluence (number of electrons per unit area).  

The provided code simulates the degradation under 1 MeV electron incidence. The parameter \(K_L\) is specific to this energy but can be analytically calculated for other incident energies <= 10 MeV using the `defectIntroductionRate.py` script.

---

## Repository Contents

1. **`iv_kphi_gaas.cpp`**  
   A C++ program that computes the degraded IV characteristics of a GaAs solar cell under irradiation. The program uses fixed K_L values for 1 MeV electrons.

2. **`defectIntroductionRate.py`**  
   A Python script to calculate \(K_L\) analytically for any incident energy <= 10 MeV.

---

## References

1. **Primary Research**  
   *"Irradiation-induced degradation in solar cells: Characterization of recombination centers"*  
   J. C. Bourgoin and M. Zazoui

2. **Associated Publication**  
   *"Modeling diffusion length damage coefficient in GaAs and InGaP solar cells under electron irradiation"*  
   J. Appl. Phys. **131**, 104503 (2022).  
   DOI: [10.1063/5.0079456](https://doi.org/10.1063/5.0079456)  
   *Author*: P. Bisht

---

## How to Use

1. **Compile and Run `iv_kphi_gaas.cpp`**  
   - Use a C++ compiler (e.g., `g++`) to compile the program:
     ```bash
     g++ iv_kphi_gaas.cpp -o iv_kphi_gaas
     ./iv_kphi_gaas
     ```
   - Adjust input parameters to simulate different fluence levels.

2. **Calculate \(K_L\) for Other Energies**  
   - Run the `defectIntroductionRate.py` script to compute \(K_L\) for specific electron energies:
     ```bash
     python defectIntroductionRate.py
     ```
   - Integrate the calculated \(K_L\) values into the C++ program as needed.

---

This repository is a useful tool for researchers studying solar cell performance under space radiation environments, with a focus on quantifying degradation due to electron irradiation.
