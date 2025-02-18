import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import os
from const import *
from gapeqn_library import sum_FPade
curr_dir = os.path.dirname(os.path.realpath(__file__))

def M(T,Gamma):
    # define the summands
    eta = lambda omega: Gamma/((1+c**2)*abs(omega))
    omegatilde = lambda omega: omega*(1+eta(omega))
    a_n = lambda w: np.arctan(Omega_c/w)
    func_a_n_w = lambda w: a_n(omegatilde(w)) / omegatilde(w) # summand for a_n_w
    func_a_n_w_2 = lambda w: a_n(omegatilde(w)) / omegatilde(w)**2
    
    # sum by Pade Decomposition
    N_FPade = 50
    a_n_w = np.real(sum_FPade(func_a_n_w,T,N_FPade))
    a_n_w_2 = np.real(sum_FPade(func_a_n_w_2,T,N_FPade))

    M_1 = np.array([[Vhh*Nh,Vhe*Ne],[Vhe*Nh,Vee*Ne]])
    M_2 = np.array([[nh*(Vhh*Nh+Vhe*Ne),ne*(Vhh*Nh+Vhe*Ne)],[nh*(Vhe*Nh+Vee*Ne),ne*(Vhe*Nh+Vee*Ne)]])


    return -2*T*(a_n_w*M_1 + Gamma/(1+c**2)*a_n_w_2*M_2)

G_max = 0.5
G_tot = np.arange(0,G_max,0.001)
Tc = np.zeros(G_tot.shape)
for i in range(len(G_tot)):
    G = G_tot[i]
    Tc[i] = fsolve(lambda T: la.det(np.identity(2)-M(T,G)), 0.1)[0]
    if i%50 == 0:
        print(f'{i} out of {len(G_tot)} computation completed' )

"""
G = 0; G_tot = []; Tc = []
while True:
    G_tot.append(G)
    Tc.append(fsolve(lambda T: la.det(np.identity(2)-M(T,G)), 0.1)[0])
    G = G + 0.01
    if Tc[-1]<1e-3: del G_tot[-1]; del Tc[-1]; break
"""
np.savetxt(curr_dir+'/Tc.txt',np.array([G_tot,Tc]))

plt.plot(G_tot,Tc)
plt.xlabel('Gamma');plt.ylabel('Tc')
plt.savefig(curr_dir+'/Tc.pdf'); plt.close()
























print(' ')
