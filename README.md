# solarcell_degradation_electron_irradiation
A simple code ("iv_kphi_gaas.cpp") to produce the degraded IV of GaAs solar cell under 1 MeV electron irradiation.

The primary mechanism of degradation is via defect creation in deep levels that act as recombination centers for minority charge carriers manifesting in a severely reducing the minority carrier diffusion length (Lo), the most important parameter for generating photocurrent. This effect is quantified by the parameter "K_L" i.e the diffusion length damage coefficient (for electrons and holes), a parameter dependent on incident energy only. The reduced L is given by the relation (approximated to first order effects only):

1/L^2 = 1/Lo^2 + KL*phi 

where phi is the total number of electrons falling on the solar cell. Here (in "iv_kphi_gaas.cpp") KL value is taken only for 1 MeV incidence but can be calculated analytically for any incident energy <= 10 MeV (see "defectIntroductionRate.py" script).

reference: Irradiation-induced degradation in solar cell: characterization of recombination centres - J C Bourgoin and M Zazoui

Publication for this work: Modeling diffusion length damage coefficient in GaAs and InGaP solar cells under electron irradiation, 
J. Appl. Phys. 131, 104503 (2022); https://doi.org/10.1063/5.0079456,
P. Bisht
