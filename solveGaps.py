import numpy as np
import matplotlib.pyplot as plt
from gapeqn_library import *
from const import *
import time
from datetime import timedelta
import pandas as pd
import os
curr_dir = os.path.dirname(os.path.realpath(__file__))
"""
The script calculates and stores results of the gap equations with impurity.
"""
savePlot = True
Gamma_arr = np.linspace(0.0,0.20, 10)
Tc_imp_arr = np.loadtxt(curr_dir+'/Tc.txt')
Tmin = 0.01
T_len = 80
start_time = time.time()
DeltaOut_list = []
save_list = [] # lists of [Gamma, T/Tc, DeltaOut, TildesTotal], with TildesTotal being N*3 arrays of omegaTilda, Delta_h_tilde, Delta_e_tilde
save_df = pd.DataFrame(columns=['Gamma', 'T/Tc', 'Delta_h','Delta_e', 'TildesTotal'])
for i in range(len(Gamma_arr)): # smart guess
    Gamma = Gamma_arr[i]
    Tc = closest_values(Tc_imp_arr[0], Tc_imp_arr[1],Gamma)[0]
    T_arr = np.linspace(0.3*Tc, Tmin, T_len)
    for j in range(len(T_arr)):
        T_eval = T_arr[j]
        normalised_temperature = T_eval / Tc
        if j == 0: # start from Tc, work its way backwards
            gapGuess = initD
        elif np.abs(DeltaOut_list[-1][0]) < 0.1: # use the last value as an initial guess
            gapGuess = initD
        else:
            gapGuess = DeltaOut_list[-1]
        
        DeltaOut, TildesTotal = SolveGap_FPade(initD,T_eval,Gamma)
        DeltaOut_list.append(DeltaOut)
        # save the results
        Delta_e, Delta_h = DeltaOut
        TildesTotal_str = np.array2string(TildesTotal, separator=',') 

        save_df.loc[i*T_len + j] = [Gamma, normalised_temperature, Delta_h, Delta_e,TildesTotal_str]
    if savePlot:
        Delta_h_list = np.array(save_df.iloc[i*T_len:(i+1)*T_len]['Delta_h'])
        Delta_e_list = np.array(save_df.iloc[i*T_len:(i+1)*T_len]['Delta_e'])
        T_normalised_list = np.array(save_df.iloc[i*T_len:(i+1)*T_len]['T/Tc'])
        plt.plot(T_normalised_list,Delta_h_list,'+',label = 'holes')
        plt.plot(T_normalised_list,Delta_e_list,'+',label = 'electrons')
        plt.xlabel('T/Tc')
        plt.ylabel('Delta')
        plt.title(f'Gap sizes for Gamma = {Gamma:.2f}')
        plt.savefig(curr_dir+f'/Gaps_Long_0911/gaps_short_zoomed_in_{i}.pdf'); plt.close()
    if i%5 == 0:
        end_time = time.time()
        formatted_time = str(timedelta(seconds= (end_time - start_time)))
        print(f'{i} out of {len(Gamma_arr)} Gamma values completed. Cumulative time: {formatted_time}.')


# Save the DataFrame to a text file
output_file = curr_dir+'/Gaps_Long_0911/gaps_short_zoomed_in.txt'
save_df.to_csv(output_file, sep='\t', index=False, header=True)