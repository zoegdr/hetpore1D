function [T,Xs,Phis,Uss,Ss,Ps,kss] = cylic_uniaxial_Lag_disp_2_COPY(params)

% WORKING 11/07/2024 13:43

% Nomenclature - all Lagrangian quantities
    %   Phi is Lagrangian porosity
    %   X is X/L
    %   T is T/Tpe
    %   U is U/L
    %   k is k/k0
    %   P is P/M
    %   S is S/M
    %   a is a/L
    %   adot is (T/L)*adot = (T/L)*da/dt
    %
    % Where
    %   Phi_0 is the initial porosity
    %   L is the resting length of the sponge
    %   Tpe = mu*(L^2)/(k0*M);
    %
    % Geometry:
    %   This code takes the left and right edges of the sponge to be at X=0 and X=1, respectively,
    %   where the left edge is free and the right edge is fixed. The code assumes 
    %   without loss of generality that a0/L=0 and b0/L=b/L=1.
    %
    % METHOD
    %   This code solves a non-linear diffusion equation of the form
    %   dPhi/dt - d/dx ( D(Phi) ds/dx ) = 0
    %   Using a finite volume method in space and Runge-Kutta (ode15s) in time
    % Where
    %   D(Phi) = k(Phi)/(1+Phi-Phi_0)
    % Boundary conditions
    %   Vs(0,T) = adot;
    %   Vf(1,T) = Vs(1,T) = 0;
    %
    % INPUTS
    %   params: Structure with parameters and options
    %     *.Astar: Dimensionless applied displacement amplitude
    %     *.omega: Dimensionless applied displacement frequency
    %     *.f0: Initial (relaxed) porosity
    %     *.stress_law: 
    %       - linear: Linear elasticity [default]
    %       - log: Hencky elasticity
    %       - neo: Neo-Hookean elasticity
    %     *.perm_law:
    %       - const: Constant permeability, where k is k/k0 = 1 [default]
    %       - KC: normalized Kozeny-Carman, where k is k/k0 = (((1-f0)^2)/(f0^3))*((f^3)/((1-f)^2))
    %     *.p: number of periods to run

% Unpack physical parameters
Phi_0 = params.Phi0;
Astar = params.Astar;
omega = params.omega;
perm_law = params.perm_law;
stress_law = params.stress_law;
M0 = params.M;
nu = params.nu;
L0 = M0*nu/(1-nu);
params.a = @(t) -(Astar/2)*(1-cos(omega*t));
params.adot = @(t) -(Astar*omega/2)*sin(omega*t);

% Unpack control parameters
N = params.N;
params.tspan = [0 params.p*2*pi/omega];
tspan = params.tspan;

a0=0;
b0=1;

% Initial condition: Start in a relaxed state
dX = (b0-a0)/N;
Xs = a0+dX/2:dX:b0-dX/2; % space
Phi_0s = Phi_0*ones(1,size(Xs,2)); % uniform porosity


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
    stress = @(Phi) 0.5*M0*(J(Phi)-1./(J(Phi)))+0.5*L0*(J(Phi)+1./(J(Phi))-2); % Neo-Hookean elasticity
    Phi_from_s = @(s) ( L0+s+ sqrt((L0+s).^2 + (M0+L0)*(M0-L0)) ) / ( M0+L0 ) - 1 + Phi_0;
else
        error('Unknown stress law.')
end

% Define a function for diffusivity
D = @(Phi) k(Phi)./(1+Phi-Phi_0);

%% Solve in time

% Options for ODE solver
    options = odeset('RelTol',1E-10, ...
                     'AbsTol',1E-10);

[T,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

base1 = zeros(size(Phis));

[Uss,Ss,Ps,kss] = ode_post(T,Phis);

%% Calculate pressure
% dp/dZ = d/dZ(log(1+Phi-Phi_0)/(1+Phi-Phi_0)
% Ps =  stress(Phis) - s_star;

    % -----------------------------------------------------
    % -----------------------------------------------------
    function Phidots = odefun(T,Phis)

        adot = params.adot(T);

        S = stress(Phis);

        % Calculate f.v. fluxes
        
        Fs = zeros(1,N+1);

        for i = 2:N
            Fs(1,i) = -D((Phis(i)+Phis(i-1))/2)*(S(i)-S(i-1))/dX;
        end
        
        % Boundary conditions
        FL = -adot;
        FR = 0;

        Fs(1,1) = FL;
        Fs(1,N+1) = FR;

        Phidots = -(Fs(2:N+1) - Fs(1:N))/dX;
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
            
            a = params.a(T(it));
            
            % Calculate displacement field from porosity field
            [Us] = U_from_Phi(Phi);
            
            % Calculate stress field from porosity field
            S = stress(Phi);
            
            % Calculate permeability field from porosity field
            ks = k(Phi);
            
            % Save
            Uss(it,:) = Us;
            Ss(it,:) = S;
            Ps(it,:) = S - S(1); % this is an approximation
            kss(it,:) = ks;
            
        end
            
    end


end