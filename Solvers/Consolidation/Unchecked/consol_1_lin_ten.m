function [t,zs,phis,wss,ps,strain,flux] = consol_1_lin_ten

% SOLVES FOR NORMALISED POROSITY

% consol_1_lin_ten solves the linearised consolidation problem,
% permeability k0 or k(phi), and uniform stiffness M
% the left wall is subject to a tension s_star.
% Grid [a,b] is constant in time

% Equation of the form
% dphi/dt - d/dz(k*dphi/dz) = 0

% Boundary conditions
% s=s_star at z=1
% ds/dz=0 at z=0

% phi is normalised porosity 
% w is displacement 
% s is stress 
% k is permeability

% Using a finite volume method, and implicit method in time
% (ode15s).

% Note on lambertw function:
% y = W(x)/x solves x=-ln(y)/y iff x>=-1/e
% Here this requires s_star <= 1/e

%% Set up

N = 100; % space discretisation
a = 0;
b = 1;

% Initiate grid
dz = (b-a)/N;
zs = a+dz/2:dz:b-dz/2; % space

% Initial condition, uniform porosity
phi_0s = zeros(1,size(zs,2));

% stress at boundary
sig_star = 0.2;
phi_star = sig_star; % uniform stiffness

% Solve in time

tspan = [0 1];
options = odeset('RelTol',1E-5,'AbsTol',1E-5);

[t,phis] = ode15s(@odefun,tspan,phi_0s',options);

% Calculate displacement
% dws/dz = phi with ws(a,t) = 0 
% or ws(z,t) = -int_a^z{phi dz}
wss = zeros(size(phis));
for j = 1:size(phis,1)
    % Integrate to find the integral values at the i+1/2 walls
    ws = cumsum(dz*(phis(j,1:end)));
    ws = [0,ws]; % add the left boundary value
    % Interpolate to find the integral values at the cell centres i
    wss(j,:) = interp1(zs(1)-dz/2:dz:zs(end)+dz/2, ws, zs);
end

% % Calculate pressure
% dp/dz = d/dz(phi) M=1 here
ps = phis - phi_star;

% Calculate strain dws/dz = phis (normalised porosity)
strain = phis;

% Calculate flux profile
flux = zeros(size(phis));
for f = 1:size(phis,1)
    % Calculate values at i+1/2 walls 
    flux_walls = log(1+phis(f,2:N))./(1+phis(f,2:N)) - log(1+phis(f,1:N-1))./(1+phis(f,1:N-1));
    % Add boundary values. Ghost points P_0 = P_1 ; P_N+1 = 2*phi_star - P_N
    flux_walls = [0,flux_walls,log(1+2*phi_star-phis(f,N))./(1+2*phi_star-phis(f,N)) - log(1+phis(f,N))./(1+phis(f,N))];
    % Interpolate to find flux values at cell centres i
    flux(f,:) = interp1(zs(1)-dz/2:dz:zs(end)+dz/2, flux_walls, zs);
end

% ODEs for phi
    function phidots = odefun(t,phis)
           
        % Calculate fluxes
        
        F_left = zeros(1,N+1);
        F_right = zeros(1,N+1);

        for i = 2:N
            F_left(1,i) = -(phis(i)-phis(i-1))/dz;
            F_right(1,i) = -(phis(i)-phis(i-1))/dz;
        end
        
        % Boundary conditions
        % phi(b,t) = phi_star
        % dphi(a,t)/dz = 0
        % Ghost points P_0 = P_1 ; P_N+1 = 2*phi_star - P_N
        F_left(1,1) = 0;
        F_right(1,N+1) = -2*(phi_star-phis(N))/dz;

        phidots = -(F_right(2:N+1) - F_left(1:N))/dz;
        phidots = phidots';

    end

end