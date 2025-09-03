# hetpore1D

Simulates uniaxial cyclic loading of a heterogeneous nonlinear poroelastic material

This repository is the supporting code for the paper Cyclic loading of a heterogeneous nonlinear poroelastic material. The solvers take a 1D heterogeneous poroelastic material, with free flux and an applied stress (AS) or applied displacement (AD) at Z=0, and no flux or displacement at Z=1. Heterogeneity is imposed as a local decrease in stiffness or permeability. Other forms of heterogeneity are also accommodated for within the options of the code, which are encoded in the "params" of the solver functions.

This repository is organised in the following way:

Parameters
- Function to assign default parameters for AS or AD
- Function to create custom parameters

Solvers
4 functions which solve for Phi (Lagrangian porosity) or phi (Eulerian porosity) and take params as input:
- 2 solvers with the AS boundary condition at Z=0 (Lagrangian and Eulerian)
- 2 solvers with the AD boundary condition at Z=0 (Lagrangian and Eulerian)

Varying_location
- Vary_l runs through AS and AD for homogeneous stiffness / permeability and fixed stiffness / permeability damage magnitude d for 3 damage locations l, and plots figures from the paper and additional ones
- cumulative_strain_n_flux_vs_Z: function to plot cumulative strain and flux against Z (Figure 3)
- strain_n_flux_vs_Z: function to plot strain and flux over one steady cycle against Z (Figure 5)

Varying_d_and_omega
- calc_metric: function to compute net strain and flux (m=1 corresponds to net strain and flux, m=2 is an alternative metric)
- plot_stiff_metrics : function to generate Figure 6
- plot_perm_metrics : function to generate Figure 7
- plot_stiff_metrics_AD_appendix: function to generate Figure 12
- plot_stiff_metrics_AD_appendix: function to generate Figure 13
- plot_freq_metrics_extra: function to generate extra figure
- Vary_d_n_l: function to run through AS and AD for homogeneous stiffness / permeability and heterogeneous stiffness /permeability for 3 damage locations and d =  [0:0.02:0.9]
- Vary_freq: function to run through AS and AD for homogeneous stiffness / permeability and heterogeneous stiffness /permeability for 3 damage locations and omega = [0.5:0.5:50]
- Vary_n_plot: runs through Vary_d_n_l and Vary_freq and generates figures from paper and additional ones

Additional figures
- relative_error calculates and plots the relative error as shown in Figure 9 
- plot_transient_phase: generates Figure 10

slanCM is a package used for the colour schemes in the plots, downloaded here: https://uk.mathworks.com/matlabcentral/fileexchange/120088-200-colormap