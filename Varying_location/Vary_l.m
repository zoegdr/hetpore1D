% This code take a 1D heterogeneous poroelastic material, and applies a cyclic load or a
% cyclic displacement at Z = 0, and no flux through Z = 1.
% Heterogeneity is imposed as a local decrease in stiffness of magnitude ds
% and location ls, or in permeability of magnitude dp and location lp.
% The code runs the solution for fixed omega and ds/dp, and three locations ls/lp 
% and plots Figure 3, 5, 11, 14 and 2 other non featured figures

% ------ NOTE BEFORE RUNNING -----
% Additional cell defined for storing increase in stiffness damage (Inc_) -
% uncomment / comment out as appropriate for this.
% Can also modify params.damage for other forms of heterogeneity.
% Can also modify other parameters.
% Figures 3,5,11 set d = [0 0.35 0.35 0.35]
% Figure 14 set d = [0 0.8 0.8 0.8]


% Create parameters
params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;
params.omega = 10; % choose fixed omega
params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = 'local-'; % 'local+' for Inc
params.lp = 0;
params.dp = 0;
params.ls = 0;
params.ds = 0;
params.v = 1/10;

params1 = params; params1.Astar = 0.2; % AS
params2 = params; params2.Astar = 0.1; % AD

d = [0 0.35 0.35 0.35]; % Figs 3,5,11
% d = [0 0.8 0.8 0.8]; % Fig 14
locs = [0 0.25 0.5 0.75];

% Cells to store solutions
Stiff_D_Ts = cell(2,4); Stiff_D_dUdZ = cell(2,4); Stiff_D_Qs = cell(2,4); Stiff_D_params = cell(2,4); % Stiff decrease, 1 = AS, 2 = AD
Perm_D_Ts = cell(2,4); Perm_D_dUdZ = cell(2,4); Perm_D_Qs = cell(2,4); Perm_D_params = cell(2,4); % Perm decrease, 1 = AS, 2 = AD
% --- Alternative stiff solutions:
% Inc_T = cell(2,4); Inc_dUdZ = cell(2,4); Inc_Q = cell(2,4); Inc_params = cell(2,4); % Stiff increase, 1 = AS, 2 = AD 


for i = 1:4
    
    % --------- Stiff ----------
    params1.ds = d(i); params2.ds = d(i);
    params1.ls = locs(i); params2.ls = locs(i);

    [params1,Ts,Zs,Phis,~,~,~,Qs,~,~] = cyclic_uniaxial_stress(params1); % AS
    Stiff_D_Ts{1,i} = Ts; Stiff_D_dUdZ{1,i} = Phis; Stiff_D_Qs{1,i} = Qs; Stiff_D_params{1,i} = params1;
    % --- Alternative:
    % Inc_T{1,i} = Ts; Inc_dUdZ{1,i} = Phis; Inc_Q{1,i} = Qs; Inc_params{1,i} = params1;

    try 
    [params2,Ts,~,Phis,~,~,~,Qs,~,~] = cyclic_uniaxial_disp(params2); % AD
    catch ME
    end
    if Ts(end) < params2.p*2*pi/params2.omega % assign NaN if breaks
       Stiff_D_Ts{2,i} = NaN; Stiff_D_dUdZ{2,i} = NaN; Stiff_D_Qs{2,i} = NaN; Stiff_D_params{2,i} = params2;
       % --- Alternative:
       % Inc_T{2,i} = NaN; Inc_dUdZ{2,i} = NaN; Inc_Q{2,i} = NaN; Inc_params{2,i} = params2;
    else  
        Stiff_D_Ts{2,i} = Ts; Stiff_D_dUdZ{2,i} = Phis; Stiff_D_Qs{2,i} = Qs; Stiff_D_params{2,i} = params2;
        % --- Alternative:
        % Inc_T{2,i} = Ts; Inc_dUdZ{2,i} = Phis; Inc_Q{2,i} = Qs; Inc_params{2,i} = params2;

    end

    % --------- Perm ----------
    params.ds = 0; params.ls = 0; % reset stiff damage to 0
    params1.dp = d(i); params2.dp = d(i);
    params1.lp = locs(i); params2.lp = locs(i);


    [params1,Ts,Zs,Phis,~,~,~,Qs,~,~] = cyclic_uniaxial_stress(params1); % AS
    Perm_D_Ts{1,i} = Ts; Perm_D_dUdZ{1,i} = Phis; Perm_D_Qs{1,i} = Qs; Perm_D_params{1,i} = params1;


    try
    [params2,Ts,~,Phis,~,~,~,Qs,~,~] = cyclic_uniaxial_disp(params2); % AD
    catch ME
    end
    if Ts(end) < params2.p*2*pi/params2.omega
       Perm_D_Ts{2,i} = NaN; Perm_D_dUdZ{2,i} = NaN; Perm_D_Qs{2,i} = NaN; Perm_D_params{2,i} = params2;
    else  
        Perm_D_Ts{2,i} = Ts; Perm_D_dUdZ{2,i} = Phis; Perm_D_Qs{2,i} = Qs; Perm_D_params{2,i} = params2;

    end

end

strain_n_flux_vs_Z(Stiff_D_Ts,Stiff_D_dUdZ,Stiff_D_Qs,Zs,Stiff_D_params,-0.3,0.3,-1,1) % Figure 3
% print('Stiff-damage-profiles-VS-Z','-dpdf','-r0')

cumulative_strain_n_flux_vs_Z(Stiff_D_Ts,Stiff_D_dUdZ,Stiff_D_Qs,Zs) % Figure 5
% print('Compare-abs-Stiff-damage-Locations','-dpdf','-r0')

strain_n_flux_vs_Z(Perm_D_Ts,Perm_D_dUdZ,Perm_D_Qs,Zs,Perm_D_params,-0.3,0.3,-1,1) % Figure 14
% print('Perm-damage-profiles-VS-Z','-dpdf','-r0')

cumulative_strain_n_flux_vs_Z(Perm_D_Ts,Perm_D_dUdZ,Perm_D_Qs,Zs) % Does not feature
% print('Compare-abs-Perm-damage-Locations','-dpdf','-r0')

% strain_n_flux_vs_Z(Inc_T,Inc_dUdZ,Inc_Q,Zs,Inc_params,-0.3,0.3,-1,1) % Does not feature
% % print('Stiff-INC-damage-profiles-VS-Z','-dpdf','-r0')
% 
% cumulative_strain_n_flux_vs_Z(Inc_T,Inc_dUdZ,Inc_Q,Zs) % Figure 11
% print('Compare-abs-Stiff-INC-damage-Locations','-dpdf','-r0')