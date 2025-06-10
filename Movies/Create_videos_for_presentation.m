% Create videos for presentation - add gitignore for this file

% Choose loading type
loading = 'stress';

% Initialise parameters
params = default_parameters(loading);

% Run simulations

    % Default
    params.omega = 10;
    [params,T1,Zs,Phi1,Uss1,S1,P1,Q1,dUdZ1,dP1,k1] = run_simulation(loading,params);

    % Damaged
    params.dp = 0.35;
    params.lp = 0.5;
    [params,T2,~,Phi2,Uss2,S2,P2,Q2,dUdZ2,dP2,k2] = run_simulation(loading,params);


% Make strain video
% format_params.ylimit = [0 0.3];
% format_params.ytick = [0 0.15 0.3];
% title = 'Strain-Mid-perm-damage-mod-loading';
% make_video_pres_format(title,Zs,Phi1,T1,Phi2,T2,format_params)

% Make flux video
format_params.ylimit = [-0.5 0.5];
format_params.ytick = [-0.5 0 0.5];
title = 'Flux-Mid-perm-damage-mod-loading';
make_video_pres_format(title,Zs,Q1,T1,Q2,T2,format_params)