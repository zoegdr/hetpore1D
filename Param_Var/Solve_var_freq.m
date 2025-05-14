% This code take a 1D heterogeneous poroelastic material, and applies a cyclic load or a
% cyclic displacement at Z = 0, and no flux through Z = 1.
% Heterogeneity is imposed as a local decrease in stiffness of magnitude ds
% and location ls, or in permeability of magnitude dp and location lp.
% The code runs the solution for various loading frequencies omega, and
% three locations ls or lp

% ------ NOTE BEFORE RUNNING -----
% In for loops:
% - For stiffness heterogeneity:
%   * i loop: params.ds = 0
%   * j loop: params.ds = [0,1) and params.ls = locs(j)
% - For permeability heterogeneity:
%   * i loop: params.dp = 0
%   * j loop: params.dp = [0,1) and params.lp = locs(j)
% Can also modify params.damage for other forms of heterogeneity.
% Can also modify other parameters.


% Create default parameters
params.N = 400;
params.p = 20;
params.Phi0 = 0.55;
params.nu = 0.3;
params.omega = 1;
params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = 'local-';
params.lp = 0;
params.dp = 0;
params.ls = 0;
params.ds = 0;
params.v = 1/10;

% Parameters to vary
omegas = 0.5:0.5:50;
d = [0 0.35 0.35 0.35];
l = [0 0.25 0.5 0.75];
c = length(omegas);

% Cells to store solutions {n,n}:
% n \ m =   1 AL   |   2 AD
% ----------------------------------
% 1     |         Stiff
% ----------------------------------
% 2     |         Perm
T{1,1} = cell(c,4); Phi{1,1} = cell(c,4); dP{1,1} = cell(c,4); k{1,1} = cell(c,4); S{1,1} = cell(c,4); U{1,1} = cell(c,4); Q{1,1} = cell(c,4); % Stiff, AL
T{1,2} = cell(c,4); Phi{1,2} = cell(c,4); dP{1,2} = cell(c,4); k{1,2} = cell(c,4); S{1,2} = cell(c,4); U{1,2} = cell(c,4); Q{1,2} = cell(c,4); % Stiff, AD
T{2,1} = cell(c,4); Phi{2,1} = cell(c,4); dP{2,1} = cell(c,4); k{2,1} = cell(c,4); S{2,1} = cell(c,4); U{2,1} = cell(c,4); Q{2,1} = cell(c,4); % Perm, AL
T{2,2} = cell(c,4); Phi{2,2} = cell(c,4); dP{2,2} = cell(c,4); k{2,2} = cell(c,4); S{2,2} = cell(c,4); U{2,2} = cell(c,4); Q{2,2} = cell(c,4); % Perm AD

% Cell of matrices to store metrics {p,q}:
% p \ q =   1 dU_int   |   2 Q_int
% ----------------------------------
% 1     |         Stiff AL
% ----------------------------------
% 2     |         Stiff AD
% ----------------------------------
% 3     |         Perm AL 
% ----------------------------------
% 4     |         Perm AD
metric_var_freq = cell(4,2);
for p = 1:4 
    for q = 1:2
        metric_var_freq{p,q} = zeros(c,4); %(omegas,locations)
    end
end

% Solve varying omegas and l
for i = 1:5 %length(omegas)

    disp(['i = ' num2str(i)]); % track progress
    params.omega = omegas(i);

    % % First calculate for no damage (i.e. uniform)
    % params.ds = 0;
    %     % -- Applied Load
    %     params.Astar = 0.2;
    %     [params,Ts,~,Phis,U_wall,Uss,~,Qs,Ss,~,Fluxes,dPs,ks] = tendon_uniaxial_cyclic_load_stressdiff(params);
    %     T{1,1}{1,1} = Ts; Phi{1,1}{1,1} = Phis; dP{1,1}{1,1} = dPs; k{1,1}{1,1} = ks; S{1,1}{1,1} = Ss; U{1,1}{1,1} = Uss; Q{1,1}{1,1} = Qs;
    %     % -- Applied Displacement
    %     params.Astar = 0.1;
    %     [params,Ts,Zs,Phis,Uss,Ss,Ps,Qs,dUdZ,dPs,ks] = cylic_uniaxial_Lag_disp_2(params);
    %     T{1,2}{1,1} = Ts; Phi{1,2}{1,1} = Phis; dP{1,2}{1,1} = dPs; k{1,2}{1,1} = ks; S{1,2}{1,1} = Ss; U{1,2}{1,1} = Uss; Q{1,2}{1,1} = Qs;
    % 
    %     % --
   
    for j = 1:4 % j = 1 is no damage, j = 2,3,4 is damage at different locations
        
        %  ----------------------------------------------------------- STIFF -----------------------------------------------------------
        % Vary stiffness damage parameters
        params.ds = d(j); params.ls = l(j);
            
            % % -- Applied Load
            % params.Astar = 0.2;
            % [params,Ts,Zs,Phis,~,Uss,~,Qs,Ss,~,~,dPs,ks] = tendon_uniaxial_cyclic_load_stressdiff(params);
            % T{1,1}{i,j} = Ts; Phi{1,1}{i,j} = Phis; dP{1,1}{i,j} = dPs; k{1,1}{i,j} = ks; S{1,1}{i,j} = Ss; U{1,1}{i,j} = Uss; Q{1,1}{i,j} = Qs; % store solution
            % metric_var_freq{1,1}(i,j) = calc_metric(Phis,Zs,Ts,1); % calculate strain metric
            % metric_var_freq{1,2}(i,j) = calc_metric(Qs,Zs,Ts,1); % calculate flux metric
            
            % -- Applied Displacement
            try
            params.Astar = 0.1;
            [params,Ts,~,Phis,Uss,Ss,~,Qs,~,dPs,ks] = cylic_uniaxial_Lag_disp(params);
            catch ME
            end
            if Ts(end) < params.p*2*pi/params.omega % if code breaks assign NaN
               T{1,2}{i,j} = NaN; Phi{1,2}{i,j} = NaN; dP{1,2}{i,j} = NaN; k{1,2}{i,j} = NaN; S{1,2}{i,j} = NaN; U{1,2}{i,j} = NaN; Q{1,2}{i,j} = NaN;
            else  
                T{1,2}{i,j} = Ts; Phi{1,2}{i,j} = Phis; dP{1,2}{i,j} = dPs; k{1,2}{i,j} = ks; S{1,2}{i,j} = Ss; U{1,2}{i,j} = Uss; Q{1,2}{i,j} = Qs;
    
            end
            metric_var_freq{2,1}(i,j) = calc_metric(Phis,Zs,Ts,1); % calculate strain metric
            metric_var_freq{2,2}(i,j) = calc_metric(Qs,Zs,Ts,1); % calculate flux metric
            % --

        %  ----------------------------------------------------------- PERM -----------------------------------------------------------
        % % Reset stiff damage to zero
        % params.ds = 0; params.ls = 0;
        % % Vary perm damage parameters
        % params.dp = d(j); params.lp = l(j);
        % 
        %     % -- Applied Load
        %     params.Astar = 0.2;
        %     [params,Ts,~,Phis,U_wall,Uss,Ps,Qs,Ss,dUdZ,Fluxes,dPs,ks] = tendon_uniaxial_cyclic_load_stressdiff(params);
        %     T{2,1}{i,j} = Ts; Phi{2,1}{i,j} = Phis; dP{2,1}{i,j} = dPs; k{2,1}{i,j} = ks; S{2,1}{i,j} = Ss; U{2,1}{i,j} = Uss; Q{2,1}{i,j} = Qs; % store solution
        %      metric_var_freq{3,1}(i,j) = calc_metric(Phis,Zs,Ts,1); % calculate strain metric
        %     metric_var_freq{3,2}(i,j) = calc_metric(Qs,Zs,Ts,1); % calculate flux metric
        % 
        %     % -- Applied Displacement
        %     try
        %     params.Astar = 0.1;
        %     [params,Ts,~,Phis,Uss,Ss,Ps,Qs,dUdZ,dPs,ks] = cylic_uniaxial_Lag_disp_2(params);
        %     catch ME
        %     end
        %     if Ts(end) < params.p*2*pi/params.omega % if code breaks assign NaN
        %        T{2,2}{i,j} = NaN; Phi{2,2}{i,j} = NaN; dP{2,2}{i,j} = NaN; k{2,2}{i,j} = NaN; S{2,2}{i,j} = NaN; U{2,2}{i,j} = NaN; Q{2,2}{i,j} = NaN;
        %     else  
        %         T{2,2}{i,j} = Ts; Phi{2,2}{i,j} = Phis; dP{2,2}{i,j} = dPs; k{2,2}{i,j} = ks; S{2,2}{i,j} = Ss; U{2,2}{i,j} = Uss; Q{2,2}{i,j} = Qs; % store solution
        %     end
        %     metric_var_freq{4,1}(i,j) = calc_metric(Phis,Zs,Ts,1); % calculate strain metric
        %     metric_var_freq{4,2}(i,j) = calc_metric(Qs,Zs,Ts,1); % calculate flux metric
        %     % --

    end

end