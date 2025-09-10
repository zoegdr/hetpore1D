function [params,Ts,Zs,Phis,Uss,Ss,Ps,Qs,dUdZ,ks] = cyclic_uniaxial_stress(params)

% UNIAXIAL NONDIMENSIONAL LAGRANGIAN FRAMEWORK

% VARIABLES AND OUTPUTS
%   * T is time
%   * Z is space
%   * Phi is NORMALISED porosity
%   * k is permeability
%   * U is solid displacement
%   * P is fluid pressure
%   * Q is flux
%   * S is Piola-Kirchoff stress tensor
%   * dUdZ is strain

% GEOMETRY / BOUNDARY CONDITIONS
% Lagrangian coordinates
% Lower boundary Z=0: free flow and applied cyclic stress S_star
% Upper boundary Z=1: no flow and no displacement (Vf=Vs=Us=Q=0)
%
%
% METHOD
%   * This code solves a non-linear advection-diffusion equation of the form
%     dPhi/dt - d/dZ (k(Phi)/(1+Phi)*dS/dZ) = 0 where 
%       - S(Phi) is calculated from stress law
%       - k(Phi) is calculated from permeability law (Kozeny Carman)
%   * Finite Volume method discretisation in space
%   * Runge-Kutta (ode15s) in time
%
% INPUTS
%   params: Structure with parameters
%       *.N: space discretisation
%       *.p: number of periods to run
%       *.Phi0: initial porosity
%       *.nu: Poisson ratio
%       *.Astar: magnitude of applied load
%       *.omega: frequency of applied load
%       *.stress_law:
%           - linear: Linear elasticity
%           - neo: Neo-Hookean elasticity
%       *.perm_law:
%           - const: Constant permeability k/k0 = 1
%           - KC: normalised Kozeny-Carman k/k0 = (Phi/Phi_0+1).^3./(1+Phi)
%       *.damage:
%           - local-/+: local decrease/increase in stiffness
%           - dec_neg/inc_neg: decrease/increase in stiffness from Z=1 to Z=0
%           - dec_pos/inc_pos: decrease/increase in stiffness from Z=0 to Z=1
%       *.lp: location of permeability damage
%       *.dp: magnitude of permeability damage
%       *.ls: location of stiffness damage
%       *.ds: magnitude of stiffness damage

%% Set up -------------------------------------------------------------

% Input parameters

if strcmp(params,'default') % default params
    clear params
    params = default_parameters('stress');
elseif strcmp(params,'edit') % option to edit parameters
    params = create_parameters;
else
    % use existing params
end

% ---------------------------------------------------------------------

% Unpack physical parameters
Phi_0 = params.Phi0; % initial uniform porosity
nu = params.nu; % Poisson ratio
A_star = params.Astar;
omega = params.omega;
params.S_star = @(t) (A_star/2)*(1-cos(omega*t));

lp = params.lp; % Damage
dp = params.dp;
ls = params.ls;
ds = params.ds;
v = params.v;

% Unpack control parameters
N = params.N;
params.tspan = [0 params.p*2*pi/omega];
tspan = params.tspan;

% Grid
a = 0;
b = 1;
dZ = (b-a)/N;
Zs = a+dZ/2:dZ:b-dZ/2; % space
base1 = zeros(1,N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a function for the permeability law
if strcmp(params.perm_law,'const')
    k = @(Phi,Z) ones(size(Phi)) .* (1-dp*exp(-(Z-lp).^2/(2*(v)^2))); % Constant
elseif strcmp(params.perm_law,'KC')
    if strcmp(params.damage,'local-')  
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (1-dp*exp(-(Z-lp).^2/(2*(v)^2))); % Local dec
    elseif strcmp(params.damage,'dec_neg')
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (-dp*Z+1);
    elseif strcmp(params.damage,'dec_pos')
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (dp*Z+1-dp);
    elseif strcmp(params.damage,'local+')  
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (1+dp*exp(-(Z-lp).^2/(2*(v)^2))); % Normalized Kozeny-Carman, simplified expression for k(Phi/J) where Phi is normalised porosity
    elseif strcmp(params.damage,'inc_neg')
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (-dp*Z+1+dp);
    elseif strcmp(params.damage,'inc_pos')
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (dp*Z+1);
    else
        k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi);
    end
else
    error('Unknown permeability law.')
end
params.k = k;

% Define a function for stiffness
if strcmp(params.damage,'local-')  
    % Local decrease
    params.M = @(Z) 1-ds*exp(-(Z-ls).^2/(2*v^2));
elseif strcmp(params.damage,'local+')
    % Local increase
    params.M = @(Z) 1+ds*exp(-(Z-ls).^2/(2*(1/16)^2));
elseif strcmp(params.damage,'dec_neg')
    % 1-->1-ds (overall decrease) 
    % LINEAR
    params.M = @(Z) -ds*Z+1;
elseif strcmp(params.damage,'dec_pos')
    % 1-ds-->1 (overall decrease)
    % LINEAR
    params.M = @(Z) ds*Z+1-ds;
elseif strcmp(params.damage,'inc_neg')
    % 1+ds-->1 (overall increase)
    params.M = @(Z) 1+ds-ds*cos(Z*pi/2);
elseif strcmp(params.damage,'inc_pos')
    % 1-->1+ds (overall increase)
    params.M = @(Z) 1+ds*cos(Z*pi/2);
else
    % no damage
    params.M = @(Z) 1;
end


% Lame's second parameter
params.L = @(Z) nu/(1-nu) * params.M(Z);

M = params.M; % Oedometric stiffness or P-wave modulus
L = params.L; % Lame's second parameter
E = @(Z) params.M(Z); % Young's modulus
C1 = @(Z) E(Z)/4/(1+nu);
D1 = @(Z) E(Z)/6/(1-2*nu);

% Define a function stress(f) for the stress law, and its inverse f_from_sig(sig).
if strcmp(params.stress_law,'linear')
    stress = @(Phi,Z) M(Z).*Phi; % Linear elasticity (pretty sure this is the expression in Lagrangian for uniaxial?)
    Phi_from_s = @(s,Z) s./M(Z);
elseif strcmp(params.stress_law,'log')
    stress = @(Phi,Z) M(Z).*log(1+Phi)./(1+Phi);
    Phi_from_s = @(s,Z) -1-lambertw(-s./M(Z))/(s/M(Z));
elseif strcmp(params.stress_law,'neo')
    stress = @(Phi,Z) 0.5*M(Z) .*( 1+Phi -1./(1+Phi) ) + 0.5*L(Z) .*( 1 + Phi+1./(1+Phi) -2 ); % Neo-Hookean elasticity
    Phi_from_s = @(s,Z) ( L(Z)+ s + sqrt((L(Z)+s).^2 +  (M(Z)+L(Z)).*(M(Z)-L(Z)) ) ) ./ ( M(Z)+L(Z) ) -1;
elseif strcmp(params.stress_law,'neo2')
    stress = @(Phi,Z) C1(Z) .*(1+Phi -1./(1+Phi)) + D1(Z) .*(1+Phi -1);
    Phi_from_s = @(s,Z) ( s + 2*D1(Z) + sqrt((s + 2*D1(Z)).^2 + 16*C1(Z).*(C1(Z)+D1(Z)) ) ) ./ 4/(C1(Z)+D1(Z)) -1;
else
    error('Unknown stress law.')
end

% Define a function for diffusivity
D = @(Phi,Z) k(Phi,Z)./(1+Phi);

%% Solve in time ------------------------------------------------------

global count;
count = 0;

options = odeset('RelTol',1E-5, 'AbsTol',1E-5);

% options = odeset('OutputFcn',@(Ts,Phis,flag) ode_check_save(Ts,Phis,flag,params),...
                     % 'RelTol',1E-5, ...
                     % 'AbsTol',1E-5);

% Initial condition: normalised initial porosity is just zero
Phi_0s = zeros(1,N);

[Ts,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

% Calculate other quantities
[Uss,Ps,Qs,Ss,dUdZ,ks] = main_ode_post(Ts,Phis);


%% ODEs for Phi -------------------------------------------------------
function Phidots = odefun(Ts,Phis)
    
    % Calculate boundary condition
    S_star = params.S_star(Ts); 
    % Phi_star = Phi_from_s(S_star,Zs(end)+dZ/2);
    
    Phi_star = Phi_from_s(S_star,Zs(1)-dZ/2);
   
    % Calculate f.v. fluxes
    
    Fs = base1;

    for i = 2:N
        Fs(1,i) = - D((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2)*(stress(Phis(i),Zs(i))-stress(Phis(i-1),Zs(i-1)))/dZ;
    end
    
    % Boundary conditions // applied load at Z=0, no flux at Z=1
        % Phi(1,t) = Phi_star
        FU = - D(Phi_star,Zs(1)-dZ/2)*2*(stress(Phis(1),Zs(1))-S_star)/dZ;
        % Vf(0,t) = Vs(0,t) = 0
        FD = 0; %0;
    
    Fs(1,1) = FU; 
    Fs(1,end) = FD;

    Phidots = -(Fs(2:N+1) - Fs(1:N))/dZ;
    Phidots = Phidots';

end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
function [Uss,Ps,Qs,Ss,dUdZ,ks] = main_ode_post(Ts,Phis)

    base2 = zeros(size(Phis));
    % Empty structures for saving
    % U_wall = base3; % displacement at cell walls
    Uss = base2; % displacement at cell centres
    Ps = base2; % pressure at cell centres
    Qs = base2; % flux at cell centres
    Ss = base2; % stress at cell centres
    dPs = base2;
        
    for it = 1:length(Ts)

        Phi = Phis(it,:);
        T = Ts(it);

        S_star = params.S_star(T); 
        Phi_star = Phi_from_s(S_star,Zs(1)-dZ/2);
        
        % Stress
        S = stress(Phi,Zs);
        
        % Displacement
        % dUs/dZ = Phi with Us(0,t) = 0 or Us(Z,t) = -int_0^Z{Phi dx}s
        % Integrate to find the integral values at the i+1/2 walls
        U = [0,cumsum(dZ*(Phi(1:end)))];
        % Interpolate to find the integral values at the cell centres i
        Us = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2,U,Zs);
        Us = flip(Us); % integration was done with flipped axis

        % Pressure
        P = S - S_star;
        
        % Flux
        % Calculate dP/dZ at i+1/2 walls 
        dP_inner_walls = (S(2:N)-S(1:N-1))/dZ;
        % Add boundary values. P(0) = 0 ; dP/dZ(1) = 0
        dP_walls = [(P(1)-0)/dZ*2,dP_inner_walls,0];
        % dP_walls = [(S(1)-S_star)/dZ*2,dP_inner_walls,0];
        % Interpolate to find dP/dZ at cell centres i
        dP = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2, dP_walls, Zs);
        Q = -k(Phi,Zs).*dP./(1+Phi);

        % Save
        Uss(it,:) = Us;
        Ss(it,:) = S;
        Qs(it,:) = Q;
        dPs(it,:) = dP;

    end
    dUdZ = Phis;
    ks = k(Phis,Zs);
end

end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Function for progress
function status = ode_check_save(Ts,Phis,flag,params)

global count;

if strcmp(flag,'done')
    disp([' - Done.'])
else
    t = Ts(1);
    count = count + 1;
    if mod(count,30)==0
        disp([' - Solving at time t = ' num2str(t) ' (' num2str(100*t/max(params.tspan)) '%) ...'])
    end
end

status = 0; % status = 0 --> Continue running

end

