To get started (after checking dependencies etc,),
1. run the two examples in solveGaps.py, explore different temperature and Gammas. For large enough Gammas and low temperatures, jumps in Delta(temp) should be observed. I got stuck on understanding those jumps. Some record can be found at the end of my diary, note that was in chronological order.

2. Go through section V of test.ipynb for calculation of heat capacities. Some data are included in this folder, just to save a bit of time in generating them. These easily take 3-6 hours (on my laptop). As a general rule of thumb, a single data point (temp, Gamma) takes about 2 seconds on my machine, though it might take longer for lower temperature. More data can be found in the old folder.

(the one band tests might be useful sanity checks. If interested, a detailed breakdown of the Pade decomposition algorithm is provided in comments of some function in gapeqn_library.py )
For B not equal to 0, a similar mean-field expansion should be possible. And a similar scheme can be set up to calculate the heat capacity. The bottleneck really are the jumps in order parameter at low temperatures and high impurity levels.


cft.py is a utility function for Pade decomposition. It allows calculation for positions and weights of poles in continued fraction representation of the Fermi-Dirac function.

const.py stores values of physical constants and phenomenological parameters (specifically, interaction strength V_ee etc.) At some point, it might be worth trying varying these parameters. 

FPade_Poles_and_Weights.npy stores information about poles of the Fermi-Dirac distribution (universal). In its absence, this will be generated automatically when calling the function to solve the gap equation. But that can take quite a while (fast for N<100, but can take hours for N ~ 1000), so better just to store it.

gapeqn_library.py contains function necessary to solve the gap equations and is imported by calculation scripts.

solveGaps.py gives two simple examples of calculating the gaps.

early_work_two_band_SC.ipynb  archived my earlier attacks at the problem, including testing the formalism on a single-band BCS superconductor, a clean twob-band superconductor, and comparing different schemes for accelerating the Matsubara sum. It should be largely irrelevant for future development, but might help illustrate my thoughts. Some theories are also sketched there.  

test.ipynb documents my progress of calculating the heat capacities. There are some detours at the beginning, but they might serve as sanity checks and useful illustration for how the code works. Each section should run independently of each other, after some global setup at the beginning.

A more complete log of results can be found in the notes subfolder, specifically the Diary markdown file. There is also a more detailed account for my derivation of the mean-field expansion in Twoband_theory.pdf, though it is not completely up to date (and far from complete).

The One_band folder applies this formalism to a single-band BCS superconductor, just as a sanity check. 

The codes are built on legacy of Maximilian Daschner, who solved the gap equations near the real(dynamical) axis and calculated thermodynamic quantities from density of states. cfr.py is from https://github.com/freude/Continued-fraction-representation . 
Yu Lu yl837, Feb 2025