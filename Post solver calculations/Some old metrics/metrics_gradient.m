% Metrics Gradient Stiffness

% In this code:
% - set parameters 
% - choose metrics to calculate (comment out those not needed)
% - metrics are saved into mass matrix [magnitude, location, omega]

params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;

params.perm_law = 'KC';
params.stress_law = 'neo';
params.lp = 0;
params.dp = 0;
params.ls = 0;

omegas = [1 10 30];
magnitudes = [0 0.2 0.5 0.8];

% ......................................................................
% Empty structures for saving


% max_strain_negative_gradient = zeros(length(magnitudes),length(omegas)); 
% max_strain_positive_gradient = zeros(length(magnitudes),length(omegas));
% max_strain_negative_gradient_disp = zeros(length(magnitudes),length(omegas)); 
max_strain_positive_gradient_disp = zeros(length(magnitudes),length(omegas));


% .....................................................................

for o = 1:length(omegas)
    params.omega = omegas(o);
    % Restrict data to one cycle
    t1 = 18*pi/params.omega;
    t2 = 20*pi/params.omega;
    for i = 1:length(magnitudes)
        params.ds = magnitudes(i);
            % 
            % % Applied load ------------------------
            % params.Astar = 0.2;
            % 
            % % Negative Gradient
            % params.damage = 'dec1';
            % [~,Ts,Zs,Phis,~,Uss,~,~,~,dUdZ,~] = tendon_uniaxial_cyclic_load_stressdiff(params); % applied load
            % [~,T1] = min(abs(Ts-t1));
            % [~,T2] = min(abs(Ts-t2));
            % max_strain_negative_gradient(i,o) = max(dUdZ(T1:T2,end));
            % 
            % % Positive Gradient
            % params.damage='dec2';
            % [~,Ts,Zs,Phis,~,Uss,~,~,~,dUdZ,~] = tendon_uniaxial_cyclic_load_stressdiff(params); % applied load
            % [~,T3] = min(abs(Ts-t1));
            % [~,T4] = min(abs(Ts-t2));
            % max_strain_positive_gradient(i,o) = max(dUdZ(T3:T4,1));

            % Applied displacement ----------------
            params.Astar = 0.1;

            % % Negative Gradient
            % params.damage = 'dec1';
            % [~,Ts,Zs,~,~,~,~,Qs,dUdZ] = cylic_uniaxial_Lag_disp_2(params); 
            % [~,T1] = min(abs(Ts-t1));
            % [~,T2] = min(abs(Ts-t2));
            % max_strain_negative_gradient_disp(i,o) = max(dUdZ(T1:T2,end));

            % Positive Gradient
            params.damage='dec2';

            [~,Ts,Zs,~,~,~,~,Qs,dUdZ] = cylic_uniaxial_Lag_disp_2(params);
            if Ts(end)<params.p*2*pi/params.omega
                warning('Simulation Failed')
                max_strain_positive_gradient_disp(i,o) = 0;
            else
                [~,T3] = min(abs(Ts-t1));
                [~,T4] = min(abs(Ts-t2));
                max_strain_positive_gradient_disp(i,o) = min(dUdZ(T3:T4,1));
            end

    end
end