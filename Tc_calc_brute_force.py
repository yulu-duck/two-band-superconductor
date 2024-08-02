import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from multiprocessing import Pool
from gapeqn_library import *
from const import *
import time
import os
import sys

curr_dir = os.path.dirname(os.path.realpath(__file__))
startT = time.time()


def SolveGap_helper(D,T,g):
    DeltaOut, Tildes_total = SolveGap(D,T,g)
    np.savetxt(curr_dir+"/results/Drange_ren_Gamma_"+"{:.5f}".format(g)+'/DeltaOut_T_'+"{:.5f}".format(T),DeltaOut)
    return DeltaOut


if __name__ == "__main__":
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
    for i in range(len(Tc_imp_arr)):
        Gamma = Tc_imp_arr[0,i]
        Tc_imp = Tc_imp_arr[1,i]
        if not os.path.exists(curr_dir+"/results/Drange_ren_Gamma_"+"{:.5f}".format(Gamma)):
            os.makedirs(curr_dir+"/results/Drange_ren_Gamma_"+"{:.5f}".format(Gamma))

        Drange_ren = p.starmap(SolveGap_helper,zip(Drange_unren,Trange,Gamma*np.ones(len(Trange))))
        np.savetxt(curr_dir+"/results/Drange_ren_Gamma_"+"{:.5f}".format(Gamma)+'/Drange_ren.txt',Drange_ren)

        Drange_ren = np.loadtxt(curr_dir+"/results/Drange_ren_Gamma_"+"{:.5f}".format(Gamma)+'/Drange_ren.txt')


        # Plot empirical formula
        t_dep = np.array([np.tanh(beta*np.sqrt(Tc_imp/t-1)) for t in Trange])
        plt.plot(Trange,Drange_ren[0,0]*t_dep,label=r'$\Delta_h\cdot \tanh{\left(\beta\sqrt{T_c/T-1}\right)}$'); plt.plot(Trange,Drange_ren[0,1]*t_dep,label=r'$\Delta_e\cdot \tanh{\left(\beta\sqrt{T_c/T-1}\right)}$')

        # Plot Delta renormalized
        plt.plot(Trange,Drange_ren[:,0],label=r'$\Delta_h$'); plt.plot(Trange,Drange_ren[:,1],label=r'$\Delta_e$')

        # Plot difference
        plt.plot(Trange,Drange_ren[:,0]-Drange_ren[0,0]*t_dep,color='k'); plt.plot(Trange,Drange_ren[:,1]-Drange_ren[0,1]*t_dep,color='k')


        plt.axvline(x=Tc,color='k',linestyle='--'); plt.text(1.01*Tc,min(min(Drange_unren[:,1]),min(Drange_ren[:,1])),'Tc')
        plt.axvline(x=Tc_imp,color='k',linestyle='--'); plt.text(0.88*Tc_imp,min(min(Drange_unren[:,1]),min(Drange_ren[:,1])),'Tc imp')
        plt.xlim(0,1.1*Tc); plt.ylim(1.1*min(min(Drange_unren[:,1]),min(Drange_ren[:,1])),1.1*max(max(Drange_unren[:,0]),max(Drange_ren[:,0])))
        plt.legend(); plt.savefig(curr_dir+'/results/Deltas_renormalized_'+"{:.5f}".format(Gamma)+'.pdf'); plt.close()











    endT = time.time()
    print('Time: '+'{:.1f}'.format(endT-startT))
