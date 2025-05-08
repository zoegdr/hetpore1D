% Metrics Local Stiffness 

% In this code:
% - set parameters 
% - choose metrics to calculate (comment out those not needed)
% - metrics are saved into mass matrix [magnitude, location, omega]

%% Set parameters

params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;
params.c=1/16;

params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = 'local-';
params.lp = 0;
params.dp = 0;

omegas = [1 10 30];
locations = [0.2 0.5 0.8];
magnitudes = [0 0.2 0.5 0.8];


% ......................................................................
%% Empty structures for saving

% max_strain_at_damage = zeros(length(magnitudes),length(locations),length(omegas));
% T_max_strain_at_damage = zeros(length(magnitudes),length(locations),length(omegas));
max_nominal_strain = zeros(length(magnitudes),length(locations),length(omegas));

% min_total_porosity = zeros(length(magnitudes),length(locations),length(omegas));
% max_total_porosity = zeros(length(magnitudes),length(locations),length(omegas));
% 
% max_flow_into_damage = zeros(length(magnitudes),length(locations),length(omegas));
% T_max_flow_into_damage = zeros(length(magnitudes),length(locations),length(omegas));

% .....................................................................

%% Calculate metrics

for o = 1:length(omegas)
    params.omega = omegas(o);
    % Restrict data to one cycle
    t1 = 18*pi/params.omega;
    t2 = 20*pi/params.omega;
    for i = 1:length(magnitudes)
        params.ds = magnitudes(i);
        for j = 1:length(locations)
            params.ls = locations(j);
            
            % Applied load ------------------------
            params.Astar = 0.2;
            [~,Ts,Zs,Phis,~,Uss,~,~,~,dUdZ,~] = tendon_uniaxial_cyclic_load_stressdiff(params); % applied load

            [~,T1] = min(abs(Ts-t1));
            [~,T2] = min(abs(Ts-t2));
            T = Ts(T1:T2);

            % Maximum Strain
            % [~,I1] = min(abs(Zs-params.ls));
            % [max_strain_at_damage(i,j,o),I2] = max(dUdZ(T1:T2,I1));
            % T_max_strain_at_damage(i,j,o) = T(I2);

            % Maximum nominal strain
            max_nominal_strain(i,j,o) = max(Uss(T1:T2,1));

            % Maximum and minimum total true porosity
            % phi_true_tot = sum((Phis+params.Phi0)./(1+Phis),2)/params.N;
            % max_total_porosity(i,j,o) = max(phi_true_tot(T1:T2));
            % min_total_porosity(i,j,o) = min(phi_true_tot(T1:T2));


            % % Applied displacement ----------------
            % params.Astar = 0.1;
            % [~,Ts,Zs,~,~,~,~,Qs,~] = cylic_uniaxial_Lag_disp_2(params);
            % 
            % [~,T3] = min(abs(Ts-t1));
            % [~,T4] = min(abs(Ts-t2));
            % T = Ts(T3:T4);
            % 
            % Max flow into damage from Z = ls-0.1 and Z=ls+0.1
            % [~,I3] = min(abs(Zs-(params.ls-0.1)));
            % [~,I4] = min(abs(Zs-(params.ls+0.1)));
            % [max_flow_into_damage(i,j,o),I5] = max(Qs(T3:T4,I3)-Qs(T3:T4,I4));
            % T_max_flow_into_damage(i,j,o) = T(I5);

        end 
    end
end