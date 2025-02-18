# Superconducting Gap Equation Solver

This repository contains a collection of utilities and scripts for solving the gap equations and calculating thermodynamic quantities in superconducting systems. The code is built on the legacy of Maximilian Daschner, who solved the gap equations near the real (dynamical) axis and calculated thermodynamic quantities from the density of states. It also includes tools for Pade decomposition of the Fermi-Dirac function.

## Getting Started

1. **Run Examples in `solveGaps.py`**:
   After ensuring that all dependencies are properly installed, run the two examples in `solveGaps.py`. Explore different temperature and Gamma values. For large enough Gamma and low temperatures, you should observe jumps in `Delta(temp)`. This behavior is still under investigation and some records related to this can be found at the end of the diary (chronological order).

2. **Calculate Heat Capacities in `test.ipynb`**:
   Go through Section V of `test.ipynb` to calculate heat capacities. Some pre-generated data are included in this folder to save time, as these calculations can take several hours (on a typical laptop). As a rough guideline, a single data point (temp, Gamma) takes about 2 seconds on my machine, although this might take longer for lower temperatures. More data can be found in the old folder.

   *Note: The "one band" tests may serve as useful sanity checks.*

3. **Pade Decomposition Breakdown**:
   If you are interested, a detailed breakdown of the Pade decomposition algorithm is provided in comments of certain functions within `gapeqn_library.py`.

## File Descriptions

- **`cft.py`**: Utility for Pade decomposition. It calculates the positions and weights of poles in the continued fraction representation of the Fermi-Dirac function.
  
- **`const.py`**: Stores values for physical constants and phenomenological parameters (e.g., interaction strength `V_ee`). You might want to experiment with varying these parameters at some point.

- **`FPade_Poles_and_Weights.npy`**: Contains precomputed poles of the Fermi-Dirac distribution (universal). If this file is absent, it will be generated automatically when solving the gap equation. However, this process can take a long time (fast for `N < 100`, but can take hours for `N ~ 1000`), so it’s better to store it once generated.

- **`gapeqn_library.py`**: Contains the necessary functions to solve the gap equations and is imported by the calculation scripts.

- **`solveGaps.py`**: Contains two simple examples of calculating gaps, to help get you started.

- **`early_work_two_band_SC.ipynb`**: Archives earlier attempts at solving the problem, including testing formalism on single-band BCS superconductors, clean two-band superconductors, and comparing different Matsubara sum acceleration schemes. While largely irrelevant for future development, it may be helpful in illustrating early thought processes and theories.

- **`test.ipynb`**: Documents progress on calculating heat capacities. Some detours are included at the beginning, but they serve as sanity checks and useful illustrations of how the code works. Each section should run independently after some global setup at the beginning.

- **`Diary.md`**: Contains a more complete log of results and observations. A more detailed account of my derivation of the mean-field expansion can also be found in **`Twoband_theory.pdf`**, though it is not up-to-date or complete.

- **`One_band/`**: Applies the formalism to a single-band BCS superconductor as a sanity check.

## Additional Notes

- **For `B ≠ 0`**: A similar mean-field expansion can be used. This can also be applied to calculate the heat capacity. The main bottleneck is the jumps in the order parameter at low temperatures and high impurity levels.

- **Legacy of Maximilian Daschner**: The code is built upon the work of Maximilian Daschner, who solved the gap equations near the real (dynamical) axis and calculated thermodynamic quantities from the density of states. The `cfr.py` file comes from [freude/Continued-fraction-representation](https://github.com/freude/Continued-fraction-representation).

## Contact Information

Yu Lu, yl837, February 2025
