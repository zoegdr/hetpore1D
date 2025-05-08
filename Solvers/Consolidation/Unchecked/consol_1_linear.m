function [t,xs,phis,uss,ps] = consol_1_linear(phi_0)

% consol_1_Lag solves the consolidation problem, for linear and non-linear
% elasticity, and permeability k0 or k(phi), and uniform stiffness M
% in a Lagrangian framework
% The left wall is subject to a stress s_star.
% Grid [a,b] is constant in time (Lagrangian)

% Equation of the form
% dP/dt - d/dx (D(phi)*ds/dx) = 0
% with relevant expression for D(phi) depending on
% elasticity and permeability law

% Boundary conditions
% s = s_star at x=0,1

% phi is porosity in Lagrangian
% U is displacement in Lagrangian
% s is stress in Lagrangian

% Using a finite volume method, and implicit method in time
% (ode15s).

% Note on lambertw function:
% y = W(x)/x solves x=-ln(y)/y iff x>=-1/e
% Here this requires s_star <= 1/e

%% Set up

N = 100; % space discretisation
a=0;
b=1;

% Initiate grid
dx = (b-a)/N;
xs = a+dx/2:dx:b-dx/2; % space

% Initial condition, uniform porosity
phi_0s = phi_0*ones(1,size(xs,2));

% stress at boundary
sig_star = -0.5;
phi_star = sig_star*(1-phi_0)+phi_0;

% Function for permeability
k = @(phi) (1-phi_0)^2/phi_0^3 * (phi./(1+phi-phi_0)).^3./(1-phi./(1+phi-phi_0)).^2;

% Solve in time

tspan = [0 1];
options = odeset('RelTol',1E-5,'AbsTol',1E-5);

[t,phis] = ode15s(@odefun,tspan,phi_0s',options);

% Calculate displacement
% dus/dx = phi-phi_0 with us(b,t) = 0 
% or us(x,t) = -int_x^b{phi-phi_0 dx}
uss = zeros(size(phis));
for j = 1:size(phis,1)
    % Integrate to find the integral values at the i-1/2 walls
    us = -(sum(dx*(phis(j,:)-phi_0)) - [0,cumsum(dx*(phis(j,1:end-1)-phi_0))]);
    us = [us,0]; % add the right boundary value
    % Interpolate to find the integral values at the cell centres i
    uss(j,:) = interp1(xs(1)-dx/2:dx:xs(end)+dx/2, us, xs);
end

% % Calculate pressure
% dp/dZ = d/dZ(log(1+phi-phi_0)/(1+phi-phi_0)
ps =  log(1+phis-phi_0)./(1+phis-phi_0) - log(1+phi_star-phi_0)/(1+phi_star-phi_0);

% ODEs for phi
    function phidots = odefun(t,phis)

        % Calculate fluxes
        
        F_left = zeros(1,N+1);
        F_right = zeros(1,N+1);

        for i = 2:N
            F_left(1,i) = -k((phis(i)+phis(i-1))/2)*(phis(i)-phis(i-1))/dx;
            F_right(1,i) = -k((phis(i)+phis(i-1))/2)*(phis(i)-phis(i-1))/dx;
        end
        
        % Boundary conditions
        % phi(a,t) = phi(b,t) = phi_star
        % Ghost points P_0 = 2*phi_star - P_1 ; P_N+1 = 2*phi_star - P_N
        F_left(1,1) = -k(phi_star)*2*(phis(1)-phi_star)/dx;
        F_right(1,N+1) = -k(phi_star)*2*(phi_star-phis(N))/dx;

        phidots = -(F_right(2:N+1) - F_left(1:N))/dx;
        phidots = phidots';

    end

end