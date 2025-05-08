% Metrics Uniform 

% In this code:
% - set parameters 
% - choose metrics to calculate (comment out those not needed)
% - metrics are saved into vectors of length omega

%% Set parameters

params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;

params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = '';
params.lp = 0;
params.dp = 0;
params.ds = 0;
params.ls = 0;

omegas = 1:1:30;

% ......................................................................
%% Empty structures for saving

min_total_porosity = zeros(length(omegas),1);
max_total_porosity = zeros(length(omegas),1);
delta_phi = zeros(length(omegas),1);

max_flow = zeros(length(omegas),1);
min_flow = zeros(length(omegas),1);
delta_Q = zeros(length(omegas),1);

max_flow_trans = zeros(length(omegas),1);
min_flow_trans = zeros(length(omegas),1);
delta_Q_trans = zeros(length(omegas),1);

% .....................................................................

%% Calculate metrics

for o = 1:length(omegas)
    params.omega = omegas(o);
    % Restrict data to one cycle
    t1 = 18*pi/params.omega;
    t2 = 20*pi/params.omega;
            
    % Applied load ------------------------
    params.Astar = 0.2;
    [~,Ts,Zs,Phis,~,Uss,~,Qs,~,dUdZ,~] = tendon_uniaxial_cyclic_load_stressdiff(params); % applied load

    [~,T1] = min(abs(Ts-t1));
    [~,T2] = min(abs(Ts-t2));
    T = Ts(T1:T2);

    % Maximum and minimum total true porosity
    phi_true_tot = sum((Phis+params.Phi0)./(1+Phis),2)/params.N;
    max_total_porosity(o) = max(phi_true_tot(T1:T2));
    min_total_porosity(o) = min(phi_true_tot(T1:T2));
    delta_phi(o) = max_total_porosity(o) - min_total_porosity(o)+params.Phi0;

    % Max and min Q
    max_flow(o) = max(Qs(T1:T2));
    min_flow(o) = min(Qs(T1:T2));
    delta_Q(o) = max_flow(o)+min_flow(o); 

    % Transient region
    t3 = 0*pi/params.omega;
    t4 = 4*pi/params.omega;
    [~,T3] = min(abs(Ts-t3));
    [~,T4] = min(abs(Ts-t4));

    % Max and min Q transient
    max_flow_trans(o) = max(Qs(T3:T4));
    min_flow_trans(o) = min(Qs(T3:T4));
    delta_Q_trans(o) = max_flow_trans(o)+min_flow_trans(o); 


end