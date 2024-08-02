import numpy as np
from scipy.optimize import fsolve

Ncores = 16

#Cutoff frequency
step = 1e-2
Omega_c = 10;
#omega_arr = np.arange(-1.2*Omega_c,1.2*Omega_c,step)
omega_arr_fine = np.arange(-Omega_c,Omega_c,step)

#Initial guess for Delta
initD = [1,-1]

#Now extend to two gaps, s+-
max_it = 10000
beta = 1.74

Vhe=200;
Vee=1;
Vhh=1;
Ne=1e-3;
Nh=2.6e-3;
Ntot=Nh+Ne;

Tmin = 2e-3
Tmax = 1

#n_imp = 6e-4;
#Gamma = n_imp/(np.pi * Ntot);
#Gamma = 0;

#c=cot(delta_0) is a measure of scattering strength.
#c=0: Unitary limit
#c>1: Born limit
c = 0; #Unitary limit

#Partial DoSes:
nh = Nh/Ntot;
ne = Ne/Ntot;

#Initial guess for Delta
def Eqn(Delta):
    ChiN = lambda Delta: -1/2*Delta*(np.log(1-Omega_c/np.sqrt(Delta**2+Omega_c**2))-np.log(1+Omega_c/np.sqrt(Delta**2+Omega_c**2)))
    return np.array([Delta[0]+(Vhh*Nh*ChiN(Delta[0]) + Vhe*Ne*ChiN(Delta[1])),Delta[1]+(Vee*Ne*ChiN(Delta[1]) + Vhe*Nh*ChiN(Delta[0]))]);
initD = fsolve(lambda D: Eqn(D),[1,-1],xtol=1e-13)

print('Delta at T=0: ',initD) # returns the e and h gap size for zero impurity level at T = 0?
