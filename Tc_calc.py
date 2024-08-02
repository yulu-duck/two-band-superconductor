import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
from const import *

curr_dir = os.path.dirname(os.path.realpath(__file__))

def M(T,Gamma):
    def eta(omega):
        return Gamma/((1+c**2)*abs(omega))

    #def varDelta(omega,Delta):
        #return Gamma/((1+c**2)*abs(omegatilde(omega)))*(nh*Delta_h+ne*Delta_e)/Delta

    def omegatilde(omega):
        return omega*(1+eta(omega))

    #def Deltatilde(omega, Delta):
        #return Delta*(1+varDelta(omega,Delta))

    #def ChiSubInt(omega,Delta):
        #return np.pi*Deltatilde(omega, Delta)/abs(omegatilde(omega))

    def omegaMatsu(n,temperature):
        omega = (2*n+1)*np.pi*temperature;
        return omega


    a_n = lambda w: np.arctan(Omega_c/w)
    a_n_w = 0
    a_n_w_2 = 0
    for n in range(-max_it,max_it):
        w_tild = omegatilde(omegaMatsu(n, T))
        a_n_w = a_n_w + a_n(w_tild)/w_tild
        a_n_w_2 = a_n_w_2 + a_n(w_tild)/w_tild**2

    M_1 = np.array([[Vhh*Nh,Vhe*Ne],[Vhe*Nh,Vee*Ne]])
    M_2 = np.array([[nh*(Vhh*Nh+Vhe*Ne),ne*(Vhh*Nh+Vhe*Ne)],[nh*(Vhe*Nh+Vee*Ne),ne*(Vhe*Nh+Vee*Ne)]])


    return -2*T*(a_n_w*M_1 + Gamma/(1+c**2)*a_n_w_2*M_2)


G = 0; G_tot = []; Tc = []
while True:
    G_tot.append(G)
    Tc.append(fsolve(lambda T: la.det(np.identity(2)-M(T,G)), 0.1)[0])
    G = G + 0.01
    if Tc[-1]<1e-3: del G_tot[-1]; del Tc[-1]; break

np.savetxt(curr_dir+'/Tc.txt',np.array([G_tot,Tc]))

plt.plot(G_tot,Tc)
plt.xlabel('Gamma');plt.ylabel('Tc')
plt.savefig(curr_dir+'/Tc.pdf'); plt.close()
























print(' ')
