# political-constraints
Fortran and Matlab codes for the project "Political Constraints and Sovereign Default Premia"

Fortran codes:
 1. main_linear_g.f90 is the main fortran code that computes and simulates the model.
 2. helper_new.f90 contains the non-calibrated parameters and the functions used in main_linear_g.f90.
 3. new_loss.f90 uses NLOPT to calibrate the remaining parameters of the model.
 
 Matlab Codes:
 1. gen_new.m uses tauchen1d.m to create a discrete state space for the aggregate TFP shock along with other discrete grid points.
 2. z_simulation_new.m creates a simulated series for the aggregate TFP shock.
 3. counterfactual.m and event_analysis.m codes the event study
 4. load_linear.m, load_simul_linear.m loads the fortran results while simulation_plots.m plots the simulated series.
