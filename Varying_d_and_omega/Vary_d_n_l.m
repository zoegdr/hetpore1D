% This code take a 1D heterogeneous poroelastic material, and applies a cyclic load or a
% cyclic displacement at Z = 0, and no flux through Z = 1.
% Heterogeneity is imposed as a local decrease in stiffness of magnitude ds
% and location ls, or in permeability of magnitude dp and location lp.
% The code runs the solution for various damage magnitudes ds/dp, and
% three locations ls/lp

% ------ NOTE BEFORE RUNNING -----
% Stiffness and permeability are both run here. Comment out if only want to
% run one or the other.
% Can also modify params.damage for other forms of heterogeneity.
% Can also modify other parameters.

function metric_var_d = Vary_d_n_l

% Create default parameters
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

% Parameters to vary
d = 0:0.02:0.9;
c = length(d);
l = [0.25 0.5 0.75];

% Cells to store solutions {n,m}:
% n \ m =   1 AL   |   2 AD
% ----------------------------------
% 1     |         Stiff
% ----------------------------------
% 2     |         Perm
T{1,1} = cell(c,3); Phi{1,1} = cell(c,3); k{1,1} = cell(c,3); S{1,1} = cell(c,3); U{1,1} = cell(c,3); Q{1,1} = cell(c,3); % Stiff, AL
T{1,2} = cell(c,3); Phi{1,2} = cell(c,3); k{1,2} = cell(c,3); S{1,2} = cell(c,3); U{1,2} = cell(c,3); Q{1,2} = cell(c,3); % Stiff, AD
T{2,1} = cell(c,3); Phi{2,1} = cell(c,3); k{2,1} = cell(c,3); S{2,1} = cell(c,3); U{2,1} = cell(c,3); Q{2,1} = cell(c,3); % Perm, AL
T{2,2} = cell(c,3); Phi{2,2} = cell(c,3); k{2,2} = cell(c,3); S{2,2} = cell(c,3); U{2,2} = cell(c,3); Q{2,2} = cell(c,3); % Perm AD

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
metric_var_d = cell(4,2);
for p = 1:4
    for q = 1:2
        metric_var_d{p,q} = zeros(c,3); %(d's,locations)
    end
end

% Solve varying d and l 

% Calculate no damage first
params.ds = 0; params.dp = 0;
% AL
    params.Astar = 0.2;
    [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_stress(params);
    T{1,1}{1,1} = Ts; Phi{1,1}{1,1} = Phis; k{1,1}{1,1} = ks; S{1,1}{1,1} = Ss; U{1,1}{1,1} = Uss; Q{1,1}{1,1} = Qs; % store solution, Stiff no damage
    T{2,1}{1,1} = Ts; Phi{2,1}{1,1} = Phis; k{2,1}{1,1} = ks; S{2,1}{1,1} = Ss; U{2,1}{1,1} = Uss; Q{2,1}{1,1} = Qs; % Perm no damage
    metric_var_d{1,1}(1,:) = calc_metric(Phis,Zs,Ts,1); % the 1 here is the choice of metric. 1 corresponds to int(int(abs)dt)dZ - See function.
    metric_var_d{1,2}(1,:) = calc_metric(Qs,Zs,Ts,1);
    metric_var_d{3,1}(1,:) = calc_metric(Phis,Zs,Ts,1); 
    metric_var_d{3,2}(1,:) = calc_metric(Qs,Zs,Ts,1);

% AD
    params.Astar = 0.1;
    [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_disp(params);
    T{1,2}{1,1} = Ts; Phi{1,2}{1,1} = Phis; k{1,2}{1,1} = ks; S{1,2}{1,1} = Ss; U{1,2}{1,1} = Uss; Q{1,2}{1,1} = Qs;
    T{2,2}{1,1} = Ts; Phi{2,2}{1,1} = Phis; k{2,2}{1,1} = ks; S{2,2}{1,1} = Ss; U{2,2}{1,1} = Uss; Q{2,2}{1,1} = Qs;
    metric_var_d{2,1}(1,:) = calc_metric(Phis,Zs,Ts,1); 
    metric_var_d{2,2}(1,:) = calc_metric(Qs,Zs,Ts,1);
    metric_var_d{4,1}(1,:) = calc_metric(Phis,Zs,Ts,1); 
    metric_var_d{4,2}(1,:) = calc_metric(Qs,Zs,Ts,1);

 % Vary d and l       
for i = 2:c
    disp(['i = ' num2str(i)]); % keep track of progress

    %  ----------------------------------------------------------- STIFF -----------------------------------------------------------
    params.ds = d(i);

    for j = 1:3
        params.ls = l(j);

        % AL
        params.Astar = 0.2;
        [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_stress(params);
        T{1,1}{i,j} = Ts; Phi{1,1}{i,j} = Phis; k{1,1}{i,j} = ks; S{1,1}{i,j} = Ss; U{1,1}{i,j} = Uss; Q{1,1}{i,j} = Qs;
        metric_var_d{1,1}(i,j) = calc_metric(Phis,Zs,Ts,1); 
        metric_var_d{1,2}(i,j) = calc_metric(Qs,Zs,Ts,1);
        % AD
        try
            params.Astar = 0.1;
            [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_disp(params);
        catch ME
        end
        if Ts(end) < params.p*2*pi/params.omega % if code breaks assign NaN
            T{1,2}{i,j} = NaN; Phi{1,2}{i,j} = NaN; k{1,2}{i,j} = NaN; S{1,2}{i,j} = NaN; U{1,2}{i,j} = NaN; Q{1,2}{i,j} = NaN;
        else  
            T{1,2}{i,j} = Ts; Phi{1,2}{i,j} = Phis; k{1,2}{i,j} = ks; S{1,2}{i,j} = Ss; U{1,2}{i,j} = Uss; Q{1,2}{i,j} = Qs;

        end
        metric_var_d{2,1}(i,j) = calc_metric(Phis,Zs,Ts,1); 
        metric_var_d{2,2}(i,j) = calc_metric(Qs,Zs,Ts,1);

    end

    %  ----------------------------------------------------------- PERM -----------------------------------------------------------
    params.ds = 0; params.ls = 0; % Reset
    params.dp = d(i);

    for j = 1:3
        params.lp = l(j);

        % AL
        params.Astar = 0.2;
        [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_stress(params);
        T{2,1}{i,j} = Ts; Phi{2,1}{i,j} = Phis; k{2,1}{i,j} = ks; S{2,1}{i,j} = Ss; U{2,1}{i,j} = Uss; Q{2,1}{i,j} = Qs;
        metric_var_d{3,1}(i,j) = calc_metric(Phis,Zs,Ts,1); 
        metric_var_d{3,2}(i,j) = calc_metric(Qs,Zs,Ts,1);
        % AD
        try
            params.Astar = 0.1;
            [params,Ts,Zs,Phis,Uss,Ss,~,Qs,~,ks] = cyclic_uniaxial_disp(params);
        catch ME
        end
        if Ts(end) < params.p*2*pi/params.omega % if code breaks assign NaN
            T{2,2}{i,j} = NaN; Phi{2,2}{i,j} = NaN; k{2,2}{i,j} = NaN; S{2,2}{i,j} = NaN; U{2,2}{i,j} = NaN; Q{2,2}{i,j} = NaN;
        else  
            T{2,2}{i,j} = Ts; Phi{2,2}{i,j} = Phis; k{2,2}{i,j} = ks; S{2,2}{i,j} = Ss; U{2,2}{i,j} = Uss; Q{2,2}{i,j} = Qs;

        end
        metric_var_d{4,1}(i,j) = calc_metric(Phis,Zs,Ts,1); 
        metric_var_d{4,2}(i,j) = calc_metric(Qs,Zs,Ts,1);

    end

end

end