import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from multiprocessing import Pool, freeze_support
from scipy.interpolate import interp1d
from gapeqn_library import *
from const import *
import time
import os

curr_dir = os.path.dirname(os.path.realpath(__file__))
startT = time.time()
w_imag = -1j*omega_arr_fine-10**(-1.1)


def SigmaCorrections(omega,Delta_h,Delta_e,Gamma):
    # Bang equations (11) and (12)
    omegatilde_cons = lambda omega,Tildes: omega + sigma_h_0(Tildes,Gamma) + sigma_e_0(Tildes,Gamma);
    Deltatilde_h_cons = lambda Delta_h,Tildes: Delta_h + sigma_h_1(Tildes,Gamma) + sigma_e_1(Tildes,Gamma);
    Deltatilde_e_cons = lambda Delta_e,Tildes: Delta_e + sigma_h_1(Tildes,Gamma) + sigma_e_1(Tildes,Gamma);

    #Self-consistent solution for omegatilde, Deltatilde_h, Deltatilde_e
    Eqn = lambda Tildes: np.array([Tildes[0]-omegatilde_cons(omega,Tildes),Tildes[1]-Deltatilde_h_cons(Delta_h,Tildes),Tildes[2]-Deltatilde_e_cons(Delta_e,Tildes)])

    TildesOut = fsolve(Eqn,[omega,Delta_h,Delta_e],xtol=1e-12);
    G = np.array([[G_h_0(TildesOut), G_h_1(TildesOut)],[ G_e_0(TildesOut), G_e_1(TildesOut)]])
    Sigma = np.array([[sigma_h_0(TildesOut,Gamma), sigma_h_1(TildesOut,Gamma)],[sigma_e_0(TildesOut,Gamma), sigma_e_1(TildesOut,Gamma)]]);

    return TildesOut.tolist()


def Sigma_Corrections_wrapper(Omega_Tilde,Delta_h,Delta_e,Gamma):
    TildesOut = np.array([SigmaCorrections(Omega_Tilde[n],Delta_h[n],Delta_e[n],Gamma) for n in range(len(Omega_Tilde))])
    G_h_imag = pade(w_imag,TildesOut[:,0],[G_h_0(t) for t in TildesOut])
    G_e_imag = pade(w_imag,TildesOut[:,0],[G_e_0(t) for t in TildesOut])
    return -Ntot*np.real(G_h_imag+G_e_imag)/(np.pi*Ntot)


def G_and_Sigma_wrapper(Omega_Tilde,Delta_h,Delta_e,Gamma):
    G_arr, Sigma_arr, Omega_Tilde, Delta_Tilde = G_and_Sigma(Omega_Tilde,Delta_h,Delta_e,Gamma)

    G_h_imag = pade(w_imag,Omega_Tilde,G_arr[:,0,0])
    G_e_imag = pade(w_imag,Omega_Tilde,G_arr[:,1,0])

    return -Ntot*np.real(G_h_imag+G_e_imag)/(np.pi*Ntot)



if __name__ == "__main__":
    freeze_support()

    # Delta unrenormalized
    Trange,Drange_unren = Delta(Tmax)
    print('There are ',len(Trange), 'temperature steps')
    if la.norm(Drange_unren[-1])<=1e-3: Tc = Trange[-1]; print('Tc: ',Tc)

    # Plot Delta unrenormalized
    plt.plot(Trange,Drange_unren[:,0],label=r'$\Delta_h$'); plt.plot(Trange,Drange_unren[:,1],label=r'$\Delta_e$')
    plt.axvline(x=Tc,color='k',linestyle='--'); plt.text(1.01*Tc,min(Drange_unren[:,1]),'Tc')
    plt.xlim(0,1.1*Tc); plt.ylim(1.1*min(Drange_unren[:,1]),1.1*max(Drange_unren[:,0]))
    plt.legend(); plt.savefig(curr_dir+'/Deltas_unrenormalized.pdf'); plt.close()


    p = Pool(Ncores)
    Tc_imp_arr = np.loadtxt(curr_dir+'/Tc.txt')

    Gamma = Tc_imp_arr[0,10]
    Tc_imp = Tc_imp_arr[1,10]
    print('Gamma: ',Gamma)
    print('Tc_imp: ',Tc_imp)

    temps = np.linspace(0.01,Tc_imp-0.01,Ncores)


    if True:
        DeltaOut, Tildes_total = SolveGap(Drange_unren[0],Trange[0],Gamma)
        #np.savetxt(curr_dir+'/DeltaOut.txt',DeltaOut)
        np.savetxt(curr_dir+'/Tildes_total.txt',Tildes_total)
    else:
        Tildes_total = np.loadtxt(curr_dir+'/Tildes_total.txt')


    Omega_Tilde = np.array([Tildes_total[i][0] for i in range(len(Tildes_total))])
    Delta_h = np.array([Tildes_total[i][1] for i in range(len(Tildes_total))])
    Delta_e = np.array([Tildes_total[i][2] for i in range(len(Tildes_total))])


    print('Length of Omega_Tilde: ',np.shape(Omega_Tilde))
    print('Min and Max of Omega_Tilde: ')
    print(min(Omega_Tilde))
    print(max(Omega_Tilde))

    Omega_Tilde_T = [Omega_Tilde for i in range(len(temps))]
    Delta_h_T = [Delta_h*T_dep(T,Tc_imp) for T in temps]
    Delta_e_T = [Delta_e*T_dep(T,Tc_imp) for T in temps]

#################################

    RESULT = p.starmap(G_and_Sigma_wrapper,zip(Omega_Tilde_T,Delta_h_T,Delta_e_T,Gamma*np.ones(len(temps))))
    np.savetxt(curr_dir+'/RESULT.txt',RESULT)




    endT = time.time()
    print('Time: '+'{:.1f}'.format(endT-startT))




'''
# Find Delta_pade through G
to = np.column_stack((Omega_Tilde,Delta_h,Delta_e))

if True:
    G_h_0_imag = pade(w_imag,Omega_Tilde,[G_h_0(x) for x in to])
    G_e_0_imag = pade(w_imag,Omega_Tilde,[G_e_0(x) for x in to])
    np.savetxt(curr_dir+'/G_h_0_imag.txt',G_h_0_imag)
    np.savetxt(curr_dir+'/G_e_0_imag.txt',G_e_0_imag)
else:
    G_h_0_imag = np.loadtxt(curr_dir+'/G_h_0_imag.txt',dtype=complex)
    G_e_0_imag = np.loadtxt(curr_dir+'/G_e_0_imag.txt',dtype=complex)


delta_h_pade = lambda w, G: (np.sqrt(nh**2 - G**2)*w)/G
delta_h_pade_arr = delta_h_pade(w_imag,G_h_0_imag)

delta_e_pade = lambda w, G: (np.sqrt(ne**2 - G**2)*w)/G
delta_e_pade_arr = delta_e_pade(w_imag,G_e_0_imag)
'''
