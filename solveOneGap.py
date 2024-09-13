import numpy as np
import matplotlib.pyplot as plt
from Onegapeqn_library import *
from const import *
import time
from datetime import timedelta
from scipy.optimize import fsolve, root_scalar
import pandas as pd
import os
curr_dir = os.path.dirname(os.path.realpath(__file__))

DeltaOut = 0.99; Gamma = 0.4

T_arr = np.linspace(0.01,0.4,60)
Delta_arr = np.zeros(len(T_arr))

for i in range(len(T_arr)):
    temperature = T_arr[i]
    Delta_val = SolveGap_FPade(DeltaOut, temperature, Gamma, N_FPade=200)[0][0]
    Delta_arr[i] = Delta_val

plt.title(f'One-band Gamma = {Gamma}')
plt.plot(T_arr, Delta_arr, '+')
plt.show()