function [T,Xs,Phis,Uss] = consol_het_stiff_Lag(a,b,Phi_0)

% set up as in Macminn

% consol_het_stiff_Lag solves the consolidation problem, for linear and 
% non-linear elasticity, and permeability k0 or k(Phi), in a Lagrangian 
% framework, with heterogeneous stiffness M
% The left wall is subject to a stress s_star.
% Grid [a,b] is constant in time (Lagrangian)

% Equation of the form
% dP/dt - d/dX (A(Phi) + D(Phi)*ds/dX) = 0
% with relevant expression for D(Phi) depending on
% elasticity and permeability law

% Phi is porosity in Lagrangian
% U is displacement in Lagrangian
% s is stress in Lagrangian

% Using a finite volume method, and implicit method in time
% (ode15s).

% Note on lambertw function:
% y = W(x)/x solves x=-ln(y)/y iff x>=-1/e
% Here this requires s_star <= 1/e

%% Set up

N = 100; % space discretisation

% Initiate grid
dX = (b-a)/N;
Xs = a+dX/2:dX:b-dX/2; % space

% Initial condition, uniform porosity
Phi_0s = Phi_0*ones(1,size(Xs,2));

% Function for permeability
k = @(Phi) (1-Phi_0)^2/Phi_0^3 * (Phi./(1+Phi-Phi_0)).^3./(1-Phi./(1+Phi-Phi_0)).^2;

% Function for stiffness
s = 0.4;
c = 0.5;
M = @(X) 1-s*exp(-(X-c)^2/(2*(1/16)^2));
dM = @(X) 2*s*(X-c)/2/(1/16)^2*exp(-(X-c)^2/2/(1/16)^2);

% Function for diffusivity
% D = @(Phi) 1/(1+Phi-Phi_0).^3; % Linear elasticity, k0
% D = @(Phi) k(Phi)./(1+Phi-Phi_0).^3; % Linear elasticity, k(Phi)
% D = @(Phi) (1-log(1+Phi-Phi_0))./(1+Phi-Phi_0).^3; % Henky elasticity, k0
D = @(Phi,X) M(X)*k(Phi).*(1-log(1+Phi-Phi_0))./(1+Phi-Phi_0).^3; % Henky elasticity, k(Phi)

% Function for A (extra term from heterogeneous stiffness)
A = @(Phi,X) k(Phi)*dM(X)*log(1+Phi)/(1+Phi)^2;

% Solve in time

tspan = [0 1];
options = odeset('RelTol',1E-5,'AbsTol',1E-5);

[T,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

% Calculate displacement
% dUs/dX = Phi-Phi_0 with Us(b,t) = 0 
% or Us(X,t) = -int_X^b{Phi-Phi_0 dx}

Uss = zeros(size(Phis));

for j = 1:size(Phis,1)
    % Integrate to find the integral values at the i-1/2 walls
    Us = -(sum(dX*(Phis(j,:)-Phi_0)) - [0,cumsum(dX*(Phis(j,1:end-1)-Phi_0))]);
    Us = [Us,0]; % add the right boundary value
    Uss(j,:) = interp1(Xs(1)-dX/2:dX:Xs(end)+dX/2, Us, Xs);
end

% ODEs for Phi
    function Phidots = odefun(T,Phis)

        s_star = -0.5;

        % Linear elasticity
        % Phi_star = 1/(1-s_star)+Phi_0-1; 
        
        % Henky elasticity
        Phi_star = Phi_0-1-lambertw(-s_star/M(Zs(end)+dZ/2))/s_star/M(Zs(end)+dZ/2);

        % Calculate fluxes
        
        F_left = zeros(1,N+1);
        F_right = zeros(1,N+1);

        for i = 2:N
            F_left(1,i) = -A((Phis(i)+Phis(i-1))/2,Xs(i)-dX/2) - D((Phis(i)+Phis(i-1))/2,Xs(i)-dX/2)*(Phis(i)-Phis(i-1))/dX;
            F_right(1,i) = -A((Phis(i)+Phis(i-1))/2,Xs(i)-dX/2) - D((Phis(i)+Phis(i-1))/2,Xs(i)-dX/2)*(Phis(i)-Phis(i-1))/dX;
        end
        
        % Boundary conditions
        % Phi(a,t) = Phi(b,t) = Phi_star
        % Ghost points P_0 = 2*Phi_star - P_1 ; P_N+1 = 2*Phi_star - P_N
        % F_left(1,1) = -A(Phi_star,Xs(1)-dX/2) - D(Phi_star,Xs(1)-dX/2)*2*(Phis(1)-Phi_star)/dX;
        % F_right(1,N+1) = -A(Phi_star,Xs(N)+dX/2) - D(Phi_star,Xs(N)+dX/2)*2*(Phi_star-Phis(N))/dX;

        % Boundary conditions
        % Phi(b,t) = Phi_star
        % dPhi(a,t)/dZ = 0
        % Ghost points P_0 = 2*Phi_star - P_1 ; P_N+1 = P_N
        F_left(1,1) = -A(Phi_star,Xs(1)-dX/2) - D(Phi_star,Xs(1)-dX/2)*2*(Phis(1)-Phi_star)/dX;
        F_right(1,N+1) = -A(Phis(N),Xs(N)+dX/2);

        Phidots = -(F_right(2:N+1) - F_left(1:N))/dX;
        Phidots = Phidots';

    end

end