import numpy as np
import numpy.linalg as la
from numpy import dot
from numpy.linalg import inv
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root
from scipy.integrate import quad, simpson
from const_oneBand import *
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

#Solve for Delta without Impurities
def Eqn(Delta,temperature):
    # equations for a clean 2-band superconductor
    energy = lambda epsi,Delta: (Delta**2 + epsi**2)**(1/2);
    GapIntegr = lambda epsi,temperature,Delta: Delta/(2*energy(epsi,Delta))*np.tanh(energy(epsi,Delta)/(2*temperature));
    ChiN = lambda temperature, Delta: 2*quad(lambda epsi: GapIntegr(epsi,temperature,Delta),0,Omega_c)[0];

    return Delta + Vhh*Nh*ChiN(temperature,Delta)
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


G_h_0 = lambda Tildes: nh * omegatilde(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_h(Tildes)**2);
G_h_1 = lambda Tildes: nh * Deltatilde_h(Tildes)/np.sqrt(omegatilde(Tildes)**2 + Deltatilde_h(Tildes)**2);


D =  lambda Tildes: c**2 + (G_h_0(Tildes))**2 + (G_h_1(Tildes))**2;

Tmat_h_0 = lambda Tildes: G_h_0(Tildes)/D(Tildes);
Tmat_h_1 = lambda Tildes: G_h_1(Tildes)/D(Tildes);

sigma_h_0  = lambda Tildes,Gamma: Gamma*Tmat_h_0(Tildes);
sigma_h_1  = lambda Tildes,Gamma: Gamma*Tmat_h_1(Tildes);

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
            print('N_FPade')
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
        if not np.isscalar(Delta):
            Delta = Delta[0] # fsolve gives an input [val], a list of length 1; but this doesn't make the SelfConsEqn vectorizable
        Delta_h1 = Delta
        ChiN_h, Tildes_total = ChiMatsu_FPade(temperature,Delta_h1,Gamma,positive_poles,weights);

        Deltacalc = Delta_h1+(Vhh*Nh*ChiN_h)
        # print(f'Evaluating SelfConsEqn, the difference is {Deltacalc}')
        return Deltacalc, Tildes_total

    wrapper_SelfConsEqn = lambda Delta: SelfConsEqn(Delta)[0]

    #Self-consistently solve the gap equation, starting with estimates Delta_h and Delta_e.
    Delta0 = DeltaOut;
    DeltaOut = fsolve(wrapper_SelfConsEqn,Delta0,xtol=1e-13);
    DeltaOut = DeltaOut.tolist()

    Tildes_total = SelfConsEqn(DeltaOut)[1]
    return DeltaOut, Tildes_total

def ChiMatsu_FPade(temperature,Delta_h,Gamma,positive_poles, weights):
    # Calculate the S.C. susceptibility by doing the Matsubara sum using Pade decomposition of Fermi-Dirac from the given positive poles and residues of the Fermi-Diract function
    ChiSubInt = lambda omega,Delta: 2*np.arctan(Omega_c/np.sqrt(omega**2+Delta**2))*Delta/np.sqrt(omega**2+Delta**2);

    Chi_h_arr = np.zeros(positive_poles.shape,dtype = np.float64)
    Tildes_total = np.zeros((np.size(positive_poles),2))
    for i in range(len(positive_poles)):
        pole = positive_poles[i]
        omega = np.real(pole * temperature) # the frequency to be evaluated

        #Find renormalised omega and Delta_h, Delta_e
        Tildes, _, _ = SigmaCorrections(omega,Delta_h,Gamma);
        Tildes_total[i,:] = Tildes

        #Bang Eq(10)
        Chi_h_arr[i] = temperature * ChiSubInt(Tildes[0],Tildes[1])
    ChiN_h = np.real(np.sum(weights * Chi_h_arr))
    # print(f'ChiN_h: {ChiN_h}')
    # print(f'Tildes_total: {Tildes_total}')
    return ChiN_h, Tildes_total


def SigmaCorrections(omega,Delta_h,Gamma):
    # Bang equations (11) and (12)
    omegatilde_cons = lambda omega,Tildes: omega + sigma_h_0(Tildes,Gamma) 
    Deltatilde_h_cons = lambda Delta_h,Tildes: Delta_h + sigma_h_1(Tildes,Gamma) 

    #Self-consistent solution for omegatilde, Deltatilde_h, Deltatilde_e given unnormalised omega and Deltas
    Eqn = lambda Tildes: np.array([Tildes[0]-omegatilde_cons(omega,Tildes),Tildes[1]-Deltatilde_h_cons(Delta_h,Tildes)])

    TildesOut = fsolve(Eqn,[omega,Delta_h],xtol=1e-12);
    G_res = np.array([[G_h_0(TildesOut), G_h_1(TildesOut)]])
    Sigma = np.array([[sigma_h_0(TildesOut,Gamma), sigma_h_1(TildesOut,Gamma)]]);

    # print(f'TildesOut: {TildesOut}')
    return TildesOut, G_res, Sigma


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