import numpy as np
import numpy.linalg as la
from numpy import dot
from numpy.linalg import inv
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root
from scipy.integrate import quad, simpson
from const import *
from scipy.interpolate import interp1d, interp2d, UnivariateSpline
from cfr import CFR # module to compute the Pade Decomposition of FD function
import os 
curr_dir = os.path.dirname(os.path.realpath(__file__))
def closest_values(arr1, arr2, val):
    # Compute the absolute difference for each value in val against all elements in arr1
    diff = np.abs(arr1[:, np.newaxis] - val)
    # Find the indices of the minimal differences for each value in val
    indices = diff.argmin(axis=0)
    # Use the indices to retrieve corresponding elements from arr2
    return arr2[indices]
def f_FD(T,w):
    if w/T > 700: return 0 #This is to prevent overflows of np.exp(). Not necessary, but avoids RuntimeWarnings.
    else: return 1/(1+np.exp(w/T))

def d_dw_f_FD(T,w):
    if T<Tmin: return 0
    elif w/T > 700: return 0 #This is to prevent overflows of np.cosh(). Not necessary, but avoids RuntimeWarnings.
    else: return -1/(2*T*(np.cosh(w/T)+1))

def T_dep(T,Tc):
    return np.tanh(beta*np.sqrt(Tc/T-1))

def d_dT_T_dep(T,Tc):
    return -(beta*Tc)/(np.sqrt((Tc - T)/T)*T**2*(1 + np.cosh(2*beta*np.sqrt(Tc/T-1))))

def omegaMatsu(n,temperature):
    omega = (2*n+1)*np.pi*temperature;
    return omega

#def Knight_shift_alt(Tildes,T,Tc):
#    G_h_0_T_imag = np.array([G_h_0_T(x,T,Tc) for x in Tildes])
#    G_e_0_T_imag = np.array([G_e_0_T(x,T,Tc) for x in Tildes])
#
#    N_T = -np.real(G_h_0_T_imag+G_e_0_T_imag)/np.pi
#
#    N = interp1d(omega_arr_fine,N_T,kind='cubic')
#
#    return quad(lambda w: -d_dw_f_FD(T,w)*N(-w),0,Omega_c-0.1) # The term -0.1 is there because quad sometimes gets outside the interpolation range


def Knight_shift(Tildes,T,Tc):
    #N = interp1d(omega_arr_fine,Ntot,kind='cubic')
    #return quad(lambda w: -d_dw_f_FD(T,w)*N(w),0,Omega_c-0.1) # The term -0.1 is there because quad sometimes gets outside the interpolation range
    bb = omega_arr_fine>0

    G_h_0_T_imag = np.array([G_h_0_T(x,T,Tc) for x in Tildes])
    G_e_0_T_imag = np.array([G_e_0_T(x,T,Tc) for x in Tildes])

    N_T = -np.real(G_h_0_T_imag+G_e_0_T_imag)/np.pi

    d_dw_f_FD_arr = np.array([d_dw_f_FD(T,w) for w in omega_arr_fine])

    return simpson(-d_dw_f_FD_arr[bb]*N_T[bb],omega_arr_fine[bb])


def specific_heat(Tildes,T,Tc):
    w_range = np.array([-np.imag(x[0]) for x in Tildes])
    bb = w_range>0

    T_arr = T*np.ones(len(w_range))

    G_h_0_T_imag = np.array([G_h_0_T(x,T,Tc) for x in Tildes])
    G_e_0_T_imag = np.array([G_e_0_T(x,T,Tc) for x in Tildes])

    d_dT_G_h_0_T_imag = np.array([d_dT_G_h_0_T(x,T,Tc) for x in Tildes])
    d_dT_G_e_0_T_imag = np.array([d_dT_G_e_0_T(x,T,Tc) for x in Tildes])

    N_T = -np.real(G_h_0_T_imag+G_e_0_T_imag)/np.pi
    d_dT_N_T = -np.real(d_dT_G_h_0_T_imag+d_dT_G_e_0_T_imag)/np.pi

    f_FD_arr = np.array([f_FD(T,w) for w in w_range])
    d_dw_f_FD_arr = np.array([d_dw_f_FD(T,w) for w in w_range])

    S_arr = (1-f_FD_arr[bb])*np.log(1-f_FD_arr[bb])+f_FD_arr[bb]*np.log(f_FD_arr[bb])

    return simpson(-d_dw_f_FD_arr[bb]*(w_range[bb]/T_arr[bb])**2*N_T[bb]-S_arr*d_dT_N_T[bb], w_range[bb])

#Solve for Delta without Impurities
def Eqn(Delta,temperature):
    # equations for a clean 2-band superconductor
    energy = lambda epsi,Delta: (Delta**2 + epsi**2)**(1/2);
    GapIntegr = lambda epsi,temperature,Delta: Delta/(2*energy(epsi,Delta))*np.tanh(energy(epsi,Delta)/(2*temperature));
    ChiN = lambda temperature, Delta: 2*quad(lambda epsi: GapIntegr(epsi,temperature,Delta),0,Omega_c)[0];

    return np.array([Delta[0]+(Vhh*Nh*ChiN(temperature,Delta[0]) + Vhe*Ne*ChiN(temperature,Delta[1])),Delta[1]+(Vee*Ne*ChiN(temperature,Delta[1]) + Vhe*Nh*ChiN(temperature,Delta[0]))]);

def Delta(T):
    # solves for the superconducting gap size (h and e) for a range of temperatures from T_min (const) to T in steps of 0.01*T_min. It returns out[0] = temperature array and out[1] = array of solution tuple (Delta_h, Delta_e). The solver starts from T_min using the gap at T = 0 as the initial guess there and works its way upwards till T or until it hits the phase transition back to normal.
    Delta = []; t=[];
    t.append(Tmin)
    Delta.append(fsolve(lambda D: Eqn(D,0.01),initD,xtol=1e-13))
    while t[-1]<T and la.norm(Delta[-1])>1e-3:
        t.append(t[-1]+0.01*t[-1]*la.norm(Delta[-1]))
        Delta.append(fsolve(lambda D: Eqn(D,t[-1]),Delta[-1],xtol=1e-13))
    Delta = np.array([i.tolist() for i in Delta])
    return t,Delta


# Bang equations (13)-(17)
omegatilde = lambda Tildes: Tildes[0];
Deltatilde_h = lambda Tildes: Tildes[1];
Deltatilde_e = lambda Tildes: Tildes[2];

G_h_0 = lambda Tildes: nh * omegatilde(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_h(Tildes)**2);
G_e_0 = lambda Tildes: ne * omegatilde(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_e(Tildes)**2);
G_h_1 = lambda Tildes: nh * Deltatilde_h(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_h(Tildes)**2);
G_e_1 = lambda Tildes: ne * Deltatilde_e(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_e(Tildes)**2);

# The following two blocks should only take Tildes at T=0 as an input
G_h_0_T = lambda Tildes, T, Tc: nh * omegatilde(Tildes)/np.sqrt(omegatilde(Tildes)**2 + (Deltatilde_h(Tildes)*T_dep(T,Tc))**2);
G_e_0_T = lambda Tildes, T, Tc: ne * omegatilde(Tildes)/np.sqrt(omegatilde(Tildes)**2 + (Deltatilde_e(Tildes)*T_dep(T,Tc))**2);
G_h_1_T = lambda Tildes, T, Tc: nh * Deltatilde_h(Tildes)*T_dep(T,Tc)/np.sqrt(omegatilde(Tildes)**2 + (Deltatilde_h(Tildes)*T_dep(T,Tc))**2);
G_e_1_T = lambda Tildes, T, Tc: ne * Deltatilde_e(Tildes)*T_dep(T,Tc)/np.sqrt(omegatilde(Tildes)**2 + (Deltatilde_e(Tildes)*T_dep(T,Tc))**2);

d_dT_G_h_0_T = lambda Tildes, T, Tc: -nh * (omegatilde(Tildes)*Deltatilde_h(Tildes)**2*T_dep(T,Tc))/(omegatilde(Tildes)**2 + (Deltatilde_h(Tildes)*T_dep(T,Tc))**2)**(3/2)*d_dT_T_dep(T,Tc);
d_dT_G_e_0_T = lambda Tildes, T, Tc: -ne * (omegatilde(Tildes)*Deltatilde_e(Tildes)**2*T_dep(T,Tc))/(omegatilde(Tildes)**2 + (Deltatilde_e(Tildes)*T_dep(T,Tc))**2)**(3/2)*d_dT_T_dep(T,Tc);
d_dT_G_h_1_T = lambda Tildes, T, Tc: nh * (omegatilde(Tildes)**2*Deltatilde_h(Tildes))/(omegatilde(Tildes)**2 + (Deltatilde_h(Tildes)*T_dep(T,Tc))**2)**(3/2)*d_dT_T_dep(T,Tc);
d_dT_G_e_1_T = lambda Tildes, T, Tc: ne * (omegatilde(Tildes)**2*Deltatilde_e(Tildes))/(omegatilde(Tildes)**2 + (Deltatilde_e(Tildes)*T_dep(T,Tc))**2)**(3/2)*d_dT_T_dep(T,Tc);

D =  lambda Tildes: c**2 + (G_h_0(Tildes) + G_e_0(Tildes))**2 + (G_h_1(Tildes) + G_e_1(Tildes))**2;

Tmat_h_0 = lambda Tildes: G_h_0(Tildes)/D(Tildes);
Tmat_e_0 = lambda Tildes: G_e_0(Tildes)/D(Tildes);
Tmat_h_1 = lambda Tildes: G_h_1(Tildes)/D(Tildes);
Tmat_e_1 = lambda Tildes: G_e_1(Tildes)/D(Tildes);

sigma_h_0  = lambda Tildes,Gamma: Gamma*Tmat_h_0(Tildes);
sigma_e_0  = lambda Tildes,Gamma: Gamma*Tmat_e_0(Tildes);
sigma_h_1  = lambda Tildes,Gamma: Gamma*Tmat_h_1(Tildes);
sigma_e_1  = lambda Tildes,Gamma: Gamma*Tmat_e_1(Tildes);


#Solve gap equation
def SolveGap(DeltaOut,temperature,Gamma):
    # returns the unnormalised gaps DeltaOut and renormalised gaps Tildes given an initial guess of unnormalised gaps
    # Bang equations (7)
    def SelfConsEqn(Delta): # self-consistent condition: Delta == Deltacalc
        Delta_h1 = Delta[0];
        Delta_e1 = Delta[1];
        ChiN_h, ChiN_e, Tildes_total = ChiMatsu(temperature,Delta_h1,Delta_e1,Gamma);

        Deltacalc = np.array([Delta_h1+(Vhh*Nh*ChiN_h + Vhe*Ne*ChiN_e),Delta_e1+(Vee*Ne*ChiN_e + Vhe*Nh*ChiN_h)]);
        return Deltacalc, Tildes_total

    wrapper_SelfConsEqn = lambda Delta: SelfConsEqn(Delta)[0]

    #Self-consistently solve the gap equation, starting with estimates Delta_h and Delta_e.
    Delta0 = np.array([DeltaOut[0], DeltaOut[1]]);
    DeltaOut = fsolve(wrapper_SelfConsEqn,Delta0,xtol=1e-13);
    DeltaOut = DeltaOut.tolist()

    Tildes_total = SelfConsEqn(DeltaOut)[1]
    return DeltaOut, Tildes_total

def ChiMatsu(temperature,Delta_h,Delta_e,Gamma):
    # IMPORTANT: the ordering of Matsubara frequencies runs from [-omega_max, ..., -omega_min, omega_min, ..., omega_max]. This is reflected in the Tildes_total generated
    # Calculate the pair susceptibility given the unnormalised order parameters from explicit summation over ALL matsubara frequencies, likely the rate-limiting step. Returns the renormalised omega_n and SC order-parameters at each Matsubara frequencies. can be improved using the Pade summation scheme. 
    
    #For each term in the sum, an integral over energy has to be evaluated:
    #ChiSubInt = 2*integral(Delta/(omega**2+epsi**2+Delta**2)) from 0 to Omega_c;

    #This is actually analytically done, and it is the same as:
    ChiSubInt = lambda omega,Delta: 2*np.arctan(Omega_c/np.sqrt(omega**2+Delta**2))*Delta/np.sqrt(omega**2+Delta**2);



    ChiN_h =0; ChiN_e = 0;
    Chi_h = []; Chi_e=[]; Tildes_total = [];
    n=0;
    #The sum should go from -infty to infty

    for n in range(-max_it(temperature),max_it(temperature)):
        omega = omegaMatsu(n,temperature);
        #Find renormalised omega and Delta_h, Delta_e
        Tildes, _, _ = SigmaCorrections(omega,Delta_h,Delta_e,Gamma);
        Tildes_total.append(Tildes.tolist())

        # Bang equation (10)
        Chi_h.append(temperature*ChiSubInt(Tildes[0],Tildes[1]));
        ChiN_h = ChiN_h + Chi_h[-1];

        Chi_e.append(temperature*ChiSubInt(Tildes[0],Tildes[2]));
        ChiN_e = ChiN_e + Chi_e[-1];

        n = n+1;



    return ChiN_h, ChiN_e, Tildes_total

def find_best_Pade_params(N_FPade, FPade_dict: dict):
    # Find the most suitable Pade decomposition roots and weights based on the cutoff number given
    # Convert the keys to integers. 
    N_arr = sorted(int(key) for key in FPade_dict.keys())

    # Find the closest key that is greater than or equal to N_FPade
    for key in N_arr:
        if key >= N_FPade:
            return FPade_dict[str(key)]
    
    # If no such key exists, print an error message
    print(f"A larger Pade decomposition cutoff {N_FPade} is needed, currently only {N_arr[-1]} is available")

def sum_FPade(func, temperature, N_FPade,args = ()):
    # WARNING: The result is generally complex. Need to recast it into real by np.real() in many cases.
    # summation using Pade decomposition for an arbitrary function func whose first argument is the Matsubara frequency, the rest of the arguments goes to args
    try: # if running in IPython environment
        FPade_dict = np.load('FPade_Poles_and_Weights.npy',allow_pickle=True).item()
    except:
        try: # if running in Python environment
            FPade_dict = np.load(curr_dir+'/FPade_Poles_and_Weights.npy',allow_pickle=True).item()
        except: # if the lookup table for weights and poles hasn't been generated before. This only works in IPython environments.
            print('Could not find pole and weights for Pade decomposition. Generating weights and poles...')
            FPade_maxN_list = np.array([5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,150,200,500,1000])
            FPade_dict = generate_FPade_poles_and_weights(FPade_maxN_list)
            print('Weights and Poles generated and saved.')
    positive_poles, weights, = find_best_Pade_params(N_FPade, FPade_dict)

    omega = positive_poles*temperature
    try: # if the function is vectorisable
        return np.sum(func(omega, *args) * weights)
    except: # else, do a for loop to sum over each term
        summed = 0
        for i in range(len(positive_poles)):
            pole = positive_poles[i]
            weight = weights[i]
            omega = pole * temperature

            summand = func(omega,*args) * weight
            summed += summand

def SolveGap_FPade(DeltaOut,temperature,Gamma,smartTruncate = True, N_FPade = None):
    # WARNING: convergence only guaranteed for T > 1mK and Gamma < 0.16
    """
    if temperature < 1e-3 or Gamma > 0.16:
        print('Convergence to 1e-6 is not guaranteed')
    """ 
    if N_FPade == None:
        N_FPade_safe = 1000 # a safe value which gurentees 1e-6 relative precision for T > 1mK and Gamma < 0.16. 
        if smartTruncate: # smartTruncate allows the Pade decomposition to truncate at earlier terms for elevated termperatures
            if temperature > 0.1:
                N_FPade = 50
            else: # a linear scaling in beta: N = 50 when T = 0.1 and N = 400 when T = 0.1
                N_FPade = np.minimum(50 + (0.10/temperature - 1) * 7,400)
        else:
            N_FPade = N_FPade_safe
    # returns the unnormalised gaps DeltaOut and renormalised gaps Tildes given an initial guess of unnormalised gaps; with a cutoff N_FPade in the Pade-decomposed Matsubara summation

    try: # if running in IPython environment
        FPade_dict = np.load('FPade_Poles_and_Weights.npy',allow_pickle=True).item()
    except:
        try: # if running in Python environment
            FPade_dict = np.load(curr_dir+'/FPade_Poles_and_Weights.npy',allow_pickle=True).item()
        except:
            print('Could not find pole and weights for Pade decomposition. Generating weights and poles...')
            FPade_maxN_list = np.array([5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,150,200,500,1000])
            FPade_dict = generate_FPade_poles_and_weights(FPade_maxN_list)
            print('Weights and Poles generated and saved.')
    positive_poles, weights, = find_best_Pade_params(N_FPade, FPade_dict)
    
    # Bang equations (7)
    def SelfConsEqn(Delta): # self-consistent condition: Delta == Deltacalc
        Delta_h1 = Delta[0];
        Delta_e1 = Delta[1];
        ChiN_h, ChiN_e, Tildes_total = ChiMatsu_FPade(temperature,Delta_h1,Delta_e1,Gamma,positive_poles,weights);

        Deltacalc = np.array([Delta_h1+(Vhh*Nh*ChiN_h + Vhe*Ne*ChiN_e),Delta_e1+(Vee*Ne*ChiN_e + Vhe*Nh*ChiN_h)]);
        return Deltacalc, Tildes_total

    wrapper_SelfConsEqn = lambda Delta: SelfConsEqn(Delta)[0]

    #Self-consistently solve the gap equation, starting with estimates Delta_h and Delta_e.
    Delta0 = np.array([DeltaOut[0], DeltaOut[1]]);
    DeltaOut = fsolve(wrapper_SelfConsEqn,Delta0,xtol=1e-13);
    DeltaOut = DeltaOut.tolist()

    Tildes_total = SelfConsEqn(DeltaOut)[1]
    return DeltaOut, Tildes_total

def ChiMatsu_FPade(temperature,Delta_h,Delta_e,Gamma,positive_poles, weights):
    # Calculate the S.C. susceptibility by doing the Matsubara sum using Pade decomposition of Fermi-Dirac from the given positive poles and residues of the Fermi-Diract function
    ChiSubInt = lambda omega,Delta: 2*np.arctan(Omega_c/np.sqrt(omega**2+Delta**2))*Delta/np.sqrt(omega**2+Delta**2);

    Chi_h_arr = np.zeros(positive_poles.shape,dtype = np.float64)
    Chi_e_arr = np.zeros(positive_poles.shape, dtype = np.float64)
    Tildes_total = np.zeros((np.size(positive_poles),3))
    for i in range(len(positive_poles)):
        pole = positive_poles[i]
        omega = np.real(pole * temperature) # the frequency to be evaluated

        #Find renormalised omega and Delta_h, Delta_e
        Tildes, _, _ = SigmaCorrections(omega,Delta_h,Delta_e,Gamma);
        Tildes_total[i] = Tildes

        #Bang Eq(10)
        Chi_h_arr[i] = temperature * ChiSubInt(Tildes[0],Tildes[1])
        Chi_e_arr[i] = temperature * ChiSubInt(Tildes[0],Tildes[2])
    ChiN_h = np.real(np.sum(weights * Chi_h_arr))
    ChiN_e = np.real(np.sum(weights * Chi_e_arr))
    return ChiN_h, ChiN_e, Tildes_total


def SigmaCorrections(omega,Delta_h,Delta_e,Gamma):
    # Bang equations (11) and (12)
    omegatilde_cons = lambda omega,Tildes: omega + sigma_h_0(Tildes,Gamma) + sigma_e_0(Tildes,Gamma);
    Deltatilde_h_cons = lambda Delta_h,Tildes: Delta_h + sigma_h_1(Tildes,Gamma) + sigma_e_1(Tildes,Gamma);
    Deltatilde_e_cons = lambda Delta_e,Tildes: Delta_e + sigma_h_1(Tildes,Gamma) + sigma_e_1(Tildes,Gamma);

    #Self-consistent solution for omegatilde, Deltatilde_h, Deltatilde_e given unnormalised omega and Deltas
    Eqn = lambda Tildes: np.array([Tildes[0]-omegatilde_cons(omega,Tildes),Tildes[1]-Deltatilde_h_cons(Delta_h,Tildes),Tildes[2]-Deltatilde_e_cons(Delta_e,Tildes)])

    TildesOut = fsolve(Eqn,[omega,Delta_h,Delta_e],xtol=1e-12);
    G_res = np.array([[G_h_0(TildesOut), G_h_1(TildesOut)],[ G_e_0(TildesOut), G_e_1(TildesOut)]])
    Sigma = np.array([[sigma_h_0(TildesOut,Gamma), sigma_h_1(TildesOut,Gamma)],[sigma_e_0(TildesOut,Gamma), sigma_e_1(TildesOut,Gamma)]]);

    return TildesOut, G_res, Sigma

def generate_FPade_poles_and_weights(FPade_maxN_list):
    """
    Generates poles and weights for the Pade decomposition of the Fermi-Dirac distribution's continued fraction representation.

    ### Description:
    The function converts a Matsubara sum over a function `f(i omega_n)` into a sum over positive poles using a Pade decomposition approach:
    
    \[
    \sum_{-\infty}^{\infty} f(i \omega_n) \approx \sum_{j=1}^{N} b_j \cdot f\left(\frac{a_j}{k_B T}\right)
    \]
    
    where:
    - \( \omega_n = (2n+1)\pi k_B T \) are Matsubara frequencies.
    - \( a_j \) are the poles, and \( b_j \) are the corresponding weights.
    
    The function returns a dictionary containing these poles and weights for each value of `N` in the input list.

    ### Parameters:
    - `FPade_maxN_list`: `list[int]`
        - A list of integers specifying the maximum number of poles (denoted by `N`) for which poles and weights should be generated.
        - Example: `[5, 10, 15, 20, 25, 30, 35, 40, 45, 50]`.

    ### Returns:
    - `FPade_params_dict`: `dict[str, np.ndarray]`
        - A dictionary where each key is a string representation of an integer `N` (from `FPade_maxN_list`), and the corresponding value is a `2xN` numpy array:
            - The first row contains the `N` positive poles.
            - The second row contains the corresponding `N` weights.
        - Example:
          ```python
          {
              "5": np.array([
                      [pole1, pole2, pole3, pole4, pole5],  # First row: Poles
                      [weight1, weight2, weight3, weight4, weight5]  # Second row: Weights
                  ]),
              ...
          }
          ```

    ### Procedure:
    1. For each `N` in `FPade_maxN_list`, calculate the poles and residues for the continued fraction representation (CFR) with a cutoff at `2*N`.
    2. Filter the calculated poles to keep only the positive values.
    3. Take the reciprocal of the positive poles.
    4. Double the corresponding residues to account for symmetry (since only positive poles are kept).
    5. Store the poles and weights in the dictionary as a `2xN` numpy array.
    6. Save the dictionary as a `.npy` file.
    
    ### Warnings:
    - The poles returned are not sorted, so all of them must be included when using the result.
    
    ### Example Usage:
    ```python
    FPade_maxN_list = [5, 10, 15, 20]
    FPade_params = generate_FPade_poles_and_weights(FPade_maxN_list)
    ```
    """
    
    FPade_params_dict = {}
    
    for N in FPade_maxN_list:
        # Set cutoff for continued fraction representation (CFR)
        cfr_cutoff = 2 * N
        
        # Generate poles and residues using CFR (assuming CFR is a predefined function)
        poles, residues = CFR(cfr_cutoff).get_poles_and_residues()
        
        # Filter positive poles and corresponding residues
        mask = poles >= 0
        filtered_poles = 1 / poles[mask]    # Take reciprocal of positive poles
        filtered_residues = 2 * residues[mask] # Double the weights
        
        # Store results in the dictionary
        result = np.array([filtered_poles, filtered_residues])
        FPade_params_dict[str(N)] = result
        print(f'Completed generating {N}')
    
    # Save the dictionary to a file
    np.save('FPade_Poles_and_Weights.npy', FPade_params_dict)
    
    return FPade_params_dict

def unnormalised_Delta_from_tildes(tildes, Gamma):
    # return Delta_e and Delta_h from a 1*3 array of [omega_tilde, Deltatilde_h, Deltatilde_e] under impurity level Gamma
    omega_tilde, Deltatilde_h, Deltatilde_e = tildes
    Delta_h = Deltatilde_h - sigma_h_1(tildes, Gamma) - sigma_e_1(tildes, Gamma)
    Delta_e = Deltatilde_e - sigma_h_1(tildes, Gamma) - sigma_e_1(tildes, Gamma) 
    omega_n = omega_tilde - sigma_h_0(tildes, Gamma) - sigma_h_0(tildes, Gamma)
    return omega_n, Delta_h, Delta_e

def sqrtM(M):
    A = M[0,0]; B = M[0,1]; C = M[1,0]; D = M[1,1]
    assert (np.sqrt(A**2+B**2) > -A),"This will result in a singular square root!"
    tau = A + D
    delta = A*D-B*C
    s = np.sqrt(delta)
    t = np.sqrt(tau + 2*s)
    return 1/t * np.array([[A+s,B],[C,D+s]])

def M(re,im):
    return np.array([[re,-im],[im,re]])

def SigmaCorrections_complex(omega,Delta_h,Delta_e,Gamma):
    #Work out the renormalisation shifts of omega and Delta
    #These are the modified frequencies and gap values

    G_h_0 = lambda v: nh * dot(M(v[0],v[1]),inv(sqrtM(dot(M(v[0],v[1]),M(v[0],v[1])) + dot(M(v[2],v[3]),M(v[2],v[3])))));
    G_e_0 = lambda v: ne * dot(M(v[0],v[1]),inv(sqrtM(dot(M(v[0],v[1]),M(v[0],v[1])) + dot(M(v[4],v[5]),M(v[4],v[5])))));
    G_h_1 = lambda v: nh * dot(M(v[2],v[3]),inv(sqrtM(dot(M(v[0],v[1]),M(v[0],v[1])) + dot(M(v[2],v[3]),M(v[2],v[3])))));
    G_e_1 = lambda v: ne * dot(M(v[4],v[5]),inv(sqrtM(dot(M(v[0],v[1]),M(v[0],v[1])) + dot(M(v[4],v[5]),M(v[4],v[5])))));

    D =  lambda v: c**2 + dot((G_h_0(v) + G_e_0(v)),(G_h_0(v) + G_e_0(v))) + dot((G_h_1(v) + G_e_1(v)),(G_h_1(v) + G_e_1(v)));

    Tmat_h_0 = lambda v: dot(G_h_0(v),inv(D(v)));
    Tmat_e_0 = lambda v: dot(G_e_0(v),inv(D(v)));
    Tmat_h_1 = lambda v: dot(G_h_1(v),inv(D(v)));
    Tmat_e_1 = lambda v: dot(G_e_1(v),inv(D(v)));

    sigma_h_0  = lambda v: Gamma*Tmat_h_0(v);
    sigma_e_0  = lambda v: Gamma*Tmat_e_0(v);
    sigma_h_1  = lambda v: Gamma*Tmat_h_1(v);
    sigma_e_1  = lambda v: Gamma*Tmat_e_1(v);

    #This is the feedback - how to make it self-consistent, given we have the Tildes
    omegatilde_cons = lambda v: M(np.real(omega),np.imag(omega)) + sigma_h_0(v) + sigma_e_0(v)
    Deltatilde_h_cons = lambda v: M(np.real(Delta_h),np.imag(Delta_h)) + sigma_h_1(v) + sigma_e_1(v);
    Deltatilde_e_cons = lambda v: M(np.real(Delta_e),np.imag(Delta_e)) + sigma_h_1(v) + sigma_e_1(v);

    #Self-consistent solution for omegatilde, Deltatilde_h, Deltatilde_e
    z1 = lambda v: M(v[0],v[1])-omegatilde_cons(v)
    z2 = lambda v: M(v[2],v[3])-Deltatilde_h_cons(v)
    z3 = lambda v: M(v[4],v[5])-Deltatilde_e_cons(v)
    Eqn = lambda v: np.array([z1(v)[0][0],z1(v)[1][0],z2(v)[0][0],z2(v)[1][0],z3(v)[0][0],z3(v)[1][0]])

    # def Eqn(v):
    #    vv = np.array([np.real(omega),v[0],v[1],v[2],v[3],v[4]])
    #    return np.array([z1(vv)[1][0],z2(vv)[0][0],z2(vv)[1][0],z3(vv)[0][0],z3(vv)[1][0]])
    #TildesOut = np.insert(TildesOut,0,np.real(omega))

    #TildesOut = fsolve(Eqn,[np.real(omega),np.imag(omega),np.real(Delta_h),np.imag(Delta_h),np.real(Delta_e),np.imag(Delta_e)],xtol=1e-12)
    TildesOut = root(Eqn,[np.real(omega),np.imag(omega),np.real(Delta_h),np.imag(Delta_h),np.real(Delta_e),np.imag(Delta_e)],method='lm').x

    G_res = np.array([[G_h_0(TildesOut)[0][0] + 1j*G_h_0(TildesOut)[1][0], G_h_1(TildesOut)[0][0] + 1j*G_h_1(TildesOut)[1][0]],[ G_e_0(TildesOut)[0][0] + 1j*G_e_0(TildesOut)[1][0], G_e_1(TildesOut)[0][0] + 1j*G_e_1(TildesOut)[1][0]]])
    Sigma = np.array([[sigma_h_0(TildesOut)[0][0] + 1j*sigma_h_0(TildesOut)[1][0], sigma_h_1(TildesOut)[0][0] + 1j*sigma_h_1(TildesOut)[1][0]],[sigma_e_0(TildesOut)[0][0] + 1j*sigma_e_0(TildesOut)[1][0], sigma_e_1(TildesOut)[0][0] + 1j*sigma_e_1(TildesOut)[1][0]]]);

    return TildesOut, G_res, Sigma

def G_and_Sigma(omega,Delta_h,Delta_e,Gamma):
    # omega = 1D array; Delta_h, Delta_e, Gamma = float
    # returns the Green's function and self-energy, and normalised order parameters given frequencies to evaluate the green's functions at and unnormalised gap size. omega can be either real (Matsubara frequencies) or complex (general frequencies)
    # returns G_arr and Sigma_arr are N*2*2 arrays, Omega_Tilde and Delta_Tilde are size N 1D array
    use_complex_frequency = np.any(np.iscomplex(omega)) # determine whether complex or real frequencies are used

    G_arr0 = np.zeros((len(omega), 2, 2),dtype = np.complex128)
    Sigma_arr0 = np.zeros((len(omega), 2, 2),dtype = np.complex128)
    Omega_Tilde = np.zeros(len(omega),dtype = np.complex128)
    Delta_Tilde = np.zeros((len(omega), 2),dtype = np.complex128)

    for i in range(len(omega)):
        if use_complex_frequency:
            Tildesout,G_res,Sigma = SigmaCorrections_complex(omega[i],Delta_h,Delta_e,Gamma)
        else:
            Tildesout,G_res,Sigma = SigmaCorrections(omega[i],Delta_h,Delta_e,Gamma)
        
        Omega_Tilde[i] = (Tildesout[0]) # store the 
        Delta_Tilde[i] = ([Tildesout[1],Tildesout[2]])
        G_arr0[i] = G_res
        Sigma_arr0[i] = Sigma

    return G_arr0, Sigma_arr0, Omega_Tilde, Delta_Tilde

def padecoeffs(z,u):
    #Calculate Pade approximant coefficients, following the procedure given in the appendix of Vidberg 1977

    #g(p,i) is g_p (z_i)
    #a_i = g_i(z_i), diagonal elements of g
    #number of coefficients N = number of points z and values u(z)

    N = len(z);
    a = np.zeros(N,dtype=np.complex128); g = np.zeros((N,N),dtype=np.complex128);

    g[0,:]=u; a[0]=g[0,0];

    for p in range(1,N):
        for i in range(p,N):
            #if np.abs(g[p-1,p-1]-g[p-1,i])<1e-10: g[p,i]=0;break
            g[p,i]=(g[p-1,p-1]-g[p-1,i])/((z[i]-z[p-1])*g[p-1,i]);
        a[p]=g[p,p];
    return a

def padeeval(z_vec,z_arr,a):
    #Calculate Pade approximation at point z_vec, given the data points z_arr and the coefficients a

    if not isinstance(z_vec, (list, tuple, np.ndarray)): z_vec = np.array([z_vec])
    N = len(a);
    c = np.zeros(len(z_vec))+0j
    f = lambda n: (z-z_arr[n])*a[n+1]

    for m in range(len(z_vec)):
        Qpoly = []
        Ppoly = []

        z = z_vec[m];
        Qpoly.append(1/(1+f(0)));
        Ppoly.append(1);

        for n in range(1,N-1):
            Qpoly.append(1/(1+f(n)*Qpoly[-1]));
            Ppoly.append(1/(1+f(n)*Ppoly[-1]));

        c[m]=a[0]*np.prod(np.array(Qpoly)/np.array(Ppoly));

    return c

def pade(x_fine,x,y):
    # Evalaute the Pade approximant at x_fine, given values at (x,y) in the complex plane
    a = padecoeffs(x, y)
    return padeeval(x_fine,x,a)

# def padeeval(z_vec,z_arr,a):
#     #Calculate Pade approximation at point z, given the data points z_arr and the coefficients a
#
#     if not isinstance(z_vec, (list, tuple, np.ndarray)): z_vec = np.array([z_vec])
#     N = len(a);
#     c = np.zeros(len(z_vec))+0j
#
#     for m in range(len(z_vec)):
#         apoly = []
#         bpoly = []
#
#         apoly.append(0); apoly.append(a[0]);
#         bpoly.append(1); bpoly.append(1);
#
#         z = z_vec[m];
#         for n in range(1,N-1):
#             apoly.append(apoly[n]+(z-z_arr[n-1])*a[n]*apoly[n-1]);
#             bpoly.append(bpoly[n]+(z-z_arr[n-1])*a[n]*bpoly[n-1]);
#             if np.abs(apoly[-2]/bpoly[-2]-apoly[-1]/bpoly[-1]) < 1e-10: break
#
#         if bpoly[-1]==0: print('denominator = 0 in Pade approximant'); c[m]=0;
#         else: c[m]=apoly[-1]/bpoly[-1];
#
#     return c