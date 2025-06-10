% Plots for paper

% Decreased stiffness plot
plot_strain_n_flux_vs_Z(Stiff_D_Ts,Stiff_D_dUdZ,Stiff_D_Qs,Zs,Stiff_D_params,-0.3,0.3,-1,1)
print('Stiff-damage-profiles-VS-Z','-dpdf','-r0')

% Decreased permeability plot
plot_strain_n_flux_vs_Z(Perm_D_Ts,Perm_D_dUdZ,Perm_D_Qs,Zs,Perm_D_params,-0.3,0.3,-1,1)
print('Perm-damage-profiles-VS-Z','-dpdf','-r0')
