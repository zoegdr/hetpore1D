% Initiate parameters
params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;
params.omega = 10;
params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = 'local-';
params.lp = 0;
params.dp = 0;
params.ls = 0;
params.ds = 0;
params.v = 1/10;

% Range of amplitudes to test
As = [0.05 0.1 0.2 0.3];

% Cells to store solutions
T = cell(2,4); dUdZ = cell(2,4); Q = cell(2,4); paramss = cell(2,4);

for i = 1:4

    params.Astar = As(i);

    [params,Ts,Zs,Phis,~,~,~,Qs,~,~,Fluxes,~,~] = tendon_uniaxial_cyclic_load_stressdiff(params); % AS
    T{1,i} = Ts; dUdZ{1,i} = Phis/params.Astar; Q{1,i} = Qs/params.Astar; paramss{1,i} = params; % Store and normalise

    try
    [params,Ts,Zs,Phis,~,~,~,Qs,~,~,~] = cylic_uniaxial_Lag_disp(params); % AD
    catch ME
    end
    if Ts(end) < params.p*2*pi/params.omega
       T{2,i} = NaN; dUdZ{2,i} = NaN; Q{2,i} = NaN; paramss{2,i} = params;
    else  
       T{2,i} = Ts; dUdZ{2,i} = Phis/params.Astar; Q{2,i} = Qs/params.Astar; paramss{2,i} = params;
    end

end

% plot_strain_n_flux_vs_Z(T,dUdZ,Q,Zs,paramss,-3,3,-6,6)