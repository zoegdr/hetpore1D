function [T,Xs,Phis,Uss,Ps] = consol_uniaxial_Lag(params)

% Scaling:
    %   f is porosity
    %   x is x/L
    %   t is t/T
    %   u is u/L
    %   p is p/M
    %   sig is sig/M
    %   q is mu*q*L/(k0*M)
    %   dp is dp/M = (p(delta)-p(L))/M
    %   a is a/L
    %   adot is (T/L)*adot = (T/L)*da/dt
    %
    % Where
    %   f0 is the initial porosity
    %   L is the initial length of the sponge
    %   T = mu*(L^2)/(k0*M);
    %
    % Geometry:
    %   This code takes the left and right edges of the sponge to be at x=a(t) and x=b=L, respectively,
    %   where the left edge is free and the right edge is fixed. The code assumes 
    %   without loss of generality that a0/L=0 and b0/L=b/L=1.
    %
    % METHOD
    %   This code solves a non-linear diffusion equation of the form
    %   dPhi/dt - d/dx ( D(Phi) ds/dx ) = 0
    %   Using a finite volume method in space and Runge-Kutta (ode15s) in time
    % Where
    %   D(Phi) = k(Phi)/(1+Phi-Phi_0
    %
    % INPUTS
    %   params: Structure with parameters and options.
    %     *.q: Dimensionless flow rate
    %     *.dp: Dimensionless pressure drop
    %       - Exactly one of the two above must be passed in. The other should be absent or NaN.
    %     *.sigstar: Dimensionless effective stress applied at x=a(t)
    %     *.f0: Initial (relaxed) porosity
    %     *.stress_law: 
    %       - linear: Linear elasticity [default]
    %       - log: Hencky elasticity
    %     *.perm_law:
    %       - const: Constant permeability, where k is k/k0 = 1 [default]
    %       - KC: normalized Kozeny-Carman, where k is k/k0 = (((1-f0)^2)/(f0^3))*((f^3)/((1-f)^2))
    %     *.p: time at which to return solution

% Unpack physical parameters
Phi_0 = params.f0;
dp = params.dp;
S_star = params.sigstar;
perm_law = params.perm_law;
stress_law = params.stress_law;
M = 1;
nu = params.nu;
L = M*nu/(1-nu);

% Unpack control parameters
N = params.N;
params.tspan = [0 params.p];
tspan = params.tspan;

a=0;
b=1;

% Initiate grid
dX = (b-a)/N;
Xs = a+dX/2:dX:b-dX/2; % space

% Initial condition, uniform porosity
Phi_0s = Phi_0*ones(1,size(Xs,2));

% Define a function for the permeability law
if strcmp(perm_law,'const')
     k = @(Phi) ones(size(Phi)); % Constant
elseif strcmp(perm_law,'KC')
    k = @(Phi) (1-Phi_0)^2/Phi_0^3 * (Phi./(1+Phi-Phi_0)).^3./(1-Phi./(1+Phi-Phi_0)).^2;
else
    error('Unknown permeability law.')
end

J = @(Phi) 1+Phi-Phi_0;

% Define a function stress(f) for the stress law, and its inverse f_from_sig(sig).
if strcmp(stress_law,'linear')
    stress = @(Phi) Phi; % Linear elasticity (pretty sure this is the expression in Lagrangian for uniaxial?)
    Phi_from_s = @(s) s;
elseif strcmp(stress_law,'log')
    stress = @(Phi) log(J(Phi))./(J(Phi));
    Phi_from_s = @(s) Phi_0-1-lambertw(-s)/s;
elseif strcmp(stress_law,'neo')
    stress = @(Phi) 0.5*M*(J(Phi)-1./(J(Phi)))+0.5*L*(J(Phi)+1./(J(Phi))-2); % Neo-Hookean elasticity
    Phi_from_s = @(s) ( L+s+ sqrt((L+s).^2 + (M+L)*(M-L)) ) / ( M+L ) - 1 + Phi_0;
else
        error('Unknown stress law.')
end

% Function for diffusivity
D = @(Phi) k(Phi)./(1+Phi-Phi_0);

% Inner boundary condition
Phi_star = Phi_from_s(S_star);

% Outer boundary condition
S_b = S_star - dp;
Phi_b = Phi_from_s(S_b);


%% Solve in time

% Options for ODE solver
    options = odeset('RelTol',1E-10, ...
                     'AbsTol',1E-10);

[T,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

base1 = zeros(size(Phis));

[Uss,Ss,Ps,kss] = ode_post(T,Phis);

     % -----------------------------------------------------
     % -----------------------------------------------------
    function Phidots = odefun(T,Phis)

        S = stress(Phis);

        % Calculate fluxes
        
        F_left = zeros(1,N+1);
        F_right = zeros(1,N+1);

        for i = 2:N
            F_left(1,i) = -D((Phis(i)+Phis(i-1))/2)*(S(i)-S(i-1))/dX;
            F_right(1,i) = -D((Phis(i)+Phis(i-1))/2)*(S(i)-S(i-1))/dX;
  
        end
        
        % Boundary conditions
        % Phi(a,t) = Phi(b,t) = Phi_star
        F_left(1,1) = -D(Phi_star)*(S(1)-S_star)/(dX/2);
        F_right(1,N+1) = -D(Phi_star)*(S_star-S(end))/(dX/2);

        Phidots = -(F_right(2:N+1) - F_left(1:N))/dX;
        Phidots = Phidots';

    end

     % -----------------------------------------------------
     % -----------------------------------------------------
    function [Us] = U_from_Phi(Phi)
        % Calculate the displacement field from the porosity
       
        % dUs/dX = Phi-Phi_0 with Us(1,t) = 0 
        % or Us(X,t) = -int_X^1{Phi-Phi_0 dx}
      
        % Integrate to find the integral values at the i-1/2 walls
        Us = -(sum(dX*(Phi-Phi_0)) - [0,cumsum(dX*(Phi(1:end-1)-Phi_0))]);
        Us = [Us,0]; % add the right boundary value
        % Interpolate to find the integral values at the cell centres i
        Us = interp1(Xs(1)-dX/2:dX:Xs(end)+dX/2, Us, Xs);

    end

    % -----------------------------------------------------
    % -----------------------------------------------------
    function [Uss,Ss,Ps,kss] = ode_post(T,Phis)
        
        % Empty structures for saving
        Uss = base1;
        Ss = base1;
        Ps = base1;
        kss = base1;

        
        % Loop over time
        for it=1:size(Phis,1)

            Phi = Phis(it,:);
            
            % a = params.a(T(it));
            
            % Calculate displacement field from porosity field
            [Us] = U_from_Phi(Phi);
            
            % Calculate stress field from porosity field
            S = stress(Phi);
            
            % Calculate permeability field from porosity field
            ks = k(Phi);
            
            % Save
            Uss(it,:) = Us;
            Ss(it,:) = S;
            Ps(it,:) = S - S_star; % this is exact
            kss(it,:) = ks;
            
        end
            
    end


end