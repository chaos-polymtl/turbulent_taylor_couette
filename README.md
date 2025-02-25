# Turbulent Taylor Couette Benchmark

This respository contains all the data and post-processing scripts that lead to the results in the publication titled "An implicit large-eddy simulation study of the turbulent Taylor-Couette flow with an inner rotating cylinder".

It contains: 
1. The results for all the simulations presented in the paper (nine cases plus two additional cases with different CFL numbers). The results include:
 * enstrophy
 * kinetic energy
 * kinetic energy rate
 * viscous dissipation
 * torque
2. Post-processing scripts to calculate the kinetic energy rate and generate each of the figures that can be found in the plots folder.
3. Reference data from Wang & Jourdan (2021). 

**Note:** all the results in this repository were obtained using Lethe v1.0 (branch: master; short: 8186f39) and deal.II v9.7.0-pre (branch: master; short: e7135bd5d1).