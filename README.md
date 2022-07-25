# solarcell_degradation_electron_irradiation
A simple code to produce the degraded IV of GaAs solar cell under 1 MeV electron irradiation.

The primary mechanism of degradation is via defectcreation in deep levels that act as recombination centers for minority charge carriers manifesting in a severely reducing the minority carrier diffusion length (L), the most import parameter for generating photocurrent. This effect is quantified by the parameter "K_L" i.e the diffusion lengthdamage coefficient (for electrons and holes), a parameter dependent on incident energy only.

1/L^2 = 1/L_o^22 + K_L*\phi 

where \phi is the total number of electrons falling on the solar cell. Here KL value is taken only for 1 MeV incidence but can be calculated analytically (different code).

see: Irradiation-induced degradation in solar cell: characterization ofrecombination centres - J C Bourgoin and M Zazoui
