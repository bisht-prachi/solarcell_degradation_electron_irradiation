#            29/07/2021. Prachi Bisht
#
#            Script for calculating defect introduction rate in GaAs under
#            electron irradiation (relativistic) at room temperature.
#            See reference "Irradiation Damage in Germanium and Silicon 
#            due to Electrons and Gamma Rays"- Julius H. Cahn, J. Appl. 
#            Phys. 30, 1310 (1959)            
         

import numpy as np
import matplotlib.pyplot as plt

Z2 = 31                          #Ga atomic number
M2 = 72                          #Ga atomic weight
na = 4.42e22                     #no. of atoms per unit volume (per cc) of GaAs
m = 1/1823                       #e- mass in amu
Ed_Ga_As = 10                    #minimum displacement energy ~ 10 eV for Ga and As
mc_squared = 0.511               #electron rest mass energy in MeV


#           see p. 1311 of reference
def fcahn(Z2,M2,E,Ed):
    alpha = Z2/137
    gamma = (E/mc_squared) + 1
    beta = np.sqrt(1 - (1/gamma**2) )                     #v/c
    Tm = 2*m*(E/M2)*(2 + E/mc_squared)*1e6                #maximum transfer energy in eV
    coeff = (2.5e-25*(Z2**2))/((beta**4)*(gamma**2))
    
    kappa = 0.8                                               #displacement efficiency for primary knock-on atom
    y = kappa*Tm/(2*Ed)                                       #available energy fraction for causing displacement
    if(y < 0.5):
        k = 0                                                 #no displacement if availabele energy < 2Ed/kappa
    
    elif (y > 0.5 and y < 1):
        term1 = 2*y - 1
        term2 = beta*(beta + np.pi*alpha)*np.log(2*y)
        term3 = 2*np.pi*alpha*beta*(1 - (np.sqrt(2*y)))
        sigma_d = coeff*( term1 - term2 - term3 )    
        k = na*sigma_d
    
    elif (y >= 1):
        term1 = np.log(y)
        term2 = beta*(beta + np.pi*alpha)*(1 - (0.3069/y))
        term3 = 2*np.pi*alpha*beta*(1 - (0.5858/np.sqrt(y)))
        sigma_d = coeff*y*( 1 + term1 - term2 + term3 )    
        k = na*sigma_d
        
        
    return k



E = 1                                                   #electron irradiation energy in MeV (don't go less than rest mass energy))

print("Incident energy, DIR(Ed = 10eV),  DIR(Ed = 21eV),  DIR(Ed = 25eV)")
while (e <= 10.0):
    k1 = fcahn(Z2,M2,E,10)                              #DIR for minimum displacement energy 10 eV
    k2 = fcahn(Z2,M2,E,21)                              #DIR for minimum displacement energy 21 eV                  
    k3 = fcahn(Z2,M2,E,25)                              #DIR for minimum displacement energy 25 eV
    print( E, k1, k2, k3)
    E += 0.5

#           Author Note: Although the minimum threshold energy, Ed, is 10 eV for Ga and As sublattice, the DIR values match the experimental
#           values for Ed = 25 eV. This suggests defect recombination processes occuring at room temperature. 
#           Author Note 2: to find damage coefficient for a property like diffusion length i.e. KL at some incident energy E > 1 MeV (provided
#           KL at 1 MeV is known), multiply the factor (k(E)/k(1 MeV)) with KL(1 MeV). Works for intermediate energy range (1-20 MeV)
#           Author note 3: Here, the rate of fall of energy of incident electrons (bit mouthful, sorry) in the sample is ignored.
#           Works if the sample is thin enough and incident energy isn't too high.
