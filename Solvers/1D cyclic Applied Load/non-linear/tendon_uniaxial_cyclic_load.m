function [params,Ts,Zs,Phis,U_wall,Uss,Ps,Qs,Ss,dUdZ,Fluxes,Adv,Diff] = tendon_uniaxial_cyclic_load(params)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UNIAXIAL LAGRANGIAN FRAMEWORK
    %
    % Variables and Outputs
    %   * T is time
    %   * Z is space
    %   * Phi is NORMALISED porosity
    %   * U is solid displacement
    %   * P is fluid pressure
    %   * Q is flux
    %   * S is Piola-Kirchoff stress tensor
    %   * dUdZ is strain
    %   * Fluxes: structure with additional calculations of flux
    %       -.flux_cum: cumulative flux 
    %       -.net_flux_Z: net flux over penultimate period at each Z
    %       -.net_flux_p: total net flux over penultimate period
    %       -.net_flux_tot: total net flux over all computed periods
    %
    % Finite volume method in space
    % Implicit in time (ode15s).

    % GEOMETRY
    %   This code takes the upper and lower edges of the material to be at Lagrangian coordinates 
    %   Z=1 and Z=0 respectively, where the upper edge is free and the lower edge is fixed. 
    %   The upper edge is subject to an applied cyclic load S_star
    %
    % METHOD
    %   * This code solves a non-linear advection-diffusion equation of the form
    %     dPhi/dt - d/dZ (A(Phi)+D(Phi)*dPhi/dZ) = 0
    %     where A(Phi) and D(Phi) are calculated from stress law
    %   * Finite Volume method in space
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
    %           - dec1/inc1: decrease/increase in stiffness from Z=1 to Z=0
    %           - dec2/inc2: decrease/increase in stiffness from Z=0 to Z=1
    %       *.lp: location of permeability damage
    %       *.dp: magnitude of permeability damage
    %       *.ls: location of stiffness damage
    %       *.ds: magnitude of stiffness damage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Set up -------------------------------------------------------------
    
     % Input parameters

    if isstruct(params) == 0
        clear params
        params.N = input("Space discretisation = ");
        params.p = input("Number of periods to run = ");
        params.Phi0 = input("Initial porosity = ");
        params.Astar = input("Magnitude of applied load = ");
        params.omega = input("Frequency of applied load = ");
        params.lp = input("Location of permeability damage = ");
        params.dp = input("Magnitude of permeability damage = ");
        params.ls = input("Location of stiffness damage = ");
        params.ds = input("Magnitude of stiffness damage = ");
        params.damage = input("Type of damage (local, inc1, dec1, inc2, dec2) = ");
        params.perm_law = input("Permeability law = ");
        params.stress_law = input("Stress law = ");
        params.nu = input("Poisson ratio = ");
    else
        question = input("Do you want to change any parameters? ");
        if strcmp(question,'yes')
            params.N = input("Space discretisation = ");
            params.p = input("Number of periods to run = ");
            params.Phi0 = input("Initial porosity = ");
            params.Astar = input("Magnitude of applied load = ");
            params.omega = input("Frequency of applied load = ");
            params.lp = input("Location of permeability damage = ");
            params.dp = input("Magnitude of permeability damage = ");
            params.ls = input("Location of stiffness damage = ");
            params.ds = input("Magnitude of stiffness damage = ");
            params.damage = input("Type of damage (local, inc1, dec1, inc2, dec2) = ");
            params.perm_law = input("Permeability law = ");
            params.stress_law = input("Stress law = ");
            params.nu = input("Poisson ratio = ");
        else
            %
        end
    end

    % ---------------------------------------------------------------------

    % Unpack physical parameters
    Phi_0 = params.Phi0; % initial uniform porosity
    nu = params.nu; % Poisson ratio
    A_star = params.Astar;
    omega = params.omega;
    lp = params.lp; % Damage
    dp = params.dp;
    ls = params.ls;
    ds = params.ds;

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
            k = @(Phi,Z) ones(size(Phi)) .* (1-dp*exp(-(Z-lp).^2/(2*(1/16)^2))); % Constant
        elseif strcmp(params.perm_law,'KC')
            k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) .* (1-dp*exp(-(Z-lp).^2/(2*(1/16)^2))); % Normalized Kozeny-Carman, simplified expression for k(Phi/J) where Phi is normalised porosity
        else
            error('Unknown permeability law.')
        end
        params.k = k;
    
    % Define a function for stiffness
        if strcmp(params.damage,'local-')  
            % Local decrease
            M = @(Z) 1-ds*exp(-(Z-ls).^2/(2*(1/16)^2));
            dM = @(Z) 2*ds*(Z-ls)/2/(1/16)^2.*exp(-(Z-ls).^2/2/(1/16)^2);
        elseif strcmp(params.damage,'local+')
            % Local increase
            M = @(Z) 1+ds*exp(-(Z-ls).^2/(2*(1/16)^2));
            dM = @(Z) -2*ds*(Z-ls)/2/(1/16)^2.*exp(-(Z-ls).^2/2/(1/16)^2);
        elseif strcmp(params.damage,'dec2')
            % ds-->1 (overall decrease) 
            M = @(Z) exp(log(1-ds)*Z);
            dM = @(Z) log(1-ds)*M(Z);
            % M = @(Z) 1-ds*cos(Z*pi/2);
            % dM = @(Z) ds*pi/2*sin(Z*pi/2);
        elseif strcmp(params.damage,'dec1')
            % 1-->ds (overall decrease)
            M = @(Z) exp(log(1-ds)*(1-Z));
            dM = @(Z) -log(1-ds)*M(Z);
            % M = @(Z) (1-ds)+ds*cos(Z*pi/2);
            % dM = @(Z) -ds*pi/2*sin(Z*pi/2);
        elseif strcmp(params.damage,'inc1')
            % 1-->1+ds (overall increase)
            M = @(Z) 1+ds-ds*cos(Z*pi/2);
            dM = @(Z) ds*pi/2*sin(Z*pi/2);
        elseif strcmp(params.damage,'inc2')
            % 1+ds-->1 (overall increase)
            M = @(Z) 1+ds*cos(Z*pi/2);
            dM = @(Z) -ds*pi/2*sin(Z*pi/2);
        else
            % no damage
            M = @(Z) 1;
            dM = @(Z) 0;
        end
 
    % Lame's second parameter
    L = @(Z) nu/(1-nu) * M(Z);
    dL = @(Z) nu/(1-nu) * dM(Z);

    % Define a function stress(f) for the stress law, and its inverse f_from_sig(sig).
    if strcmp(params.stress_law,'linear')
        stress = @(Phi,Z) M(Z)*Phi; % Linear elasticity (pretty sure this is the expression in Lagrangian for uniaxial?)
        Phi_from_s = @(s,Z) s./M(Z);
    elseif strcmp(params.stress_law,'log')
        stress = @(Phi,Z) M(Z)*log(1+Phi)./(1+Phi);
        Phi_from_s = @(s,Z) -1-lambertw(-s./M(Z))/(s/M(Z));
    elseif strcmp(params.stress_law,'neo')
        stress = @(Phi,Z) 0.5*M(Z) .*( 1+ Phi-1./(1+Phi) ) + 0.5*L(Z) .*( 1 + Phi+1./(1+Phi) -2 ); % Neo-Hookean elasticity
        Phi_from_s = @(s,Z) ( L(Z)+ s + sqrt((L(Z)+s).^2 + ( M(Z)+L(Z)).*(M(Z)-L(Z))) ) ./ ( M(Z)+L(Z) ) -1;
    
    else
        error('Unknown stress law.')
    end

    % Define a function for diffusivity
        if strcmp(params.stress_law,'neo')
            D = @(Phi,Z) ( 0.5*k(Phi,Z)./(1+Phi) ) .*( M(Z)+L(Z)+(M(Z)-L(Z))./(1+Phi).^2 ); % Neo-Hookean elasticity
        elseif strcmp(params.stress_law,'linear')
            D = @(Phi,Z) M(Z).*k(Phi,Z)./(1+Phi); % Linear elasticity law in non-linear dynamics
        else
            error(['Unknown stress law: ' stress_law])
        end
    
    % Define a function for advection (extra advection term from heterogenous stiffness)
        if strcmp(params.stress_law,'neo')
            A = @(Phi,Z) 0.5*k(Phi,Z)./(1+Phi).*(dM(Z).*(1+Phi-1./(1+Phi))+dL(Z).*(1+Phi+1./(1+Phi)-2)); % Neo-Hookean elasticity
        elseif strcmp(params.stress_law,'linear')
            A = @(Phi,Z) dM(Z).*k(Phi,Z).*Phi./(1+Phi); % Linear elasticity law in non-linear dynamics
        else
            error(['Unknown stress law: ' stress_law])
        end

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

    base2 = zeros(size(Phis));
    base3 = zeros(length(Ts),N+1);

    % Calculate other quantities
    [U_wall,Uss,Ps,Qs,Ss,dUdZ] = main_ode_post(Ts,Phis);
    [Fluxes] = calc_flux_quantities(Qs,params.p);
    [Adv,Diff] = calc_adv_diff_coeffs(Phis,Zs);


    %% ODEs for Phi -------------------------------------------------------
    function Phidots = odefun(Ts,Phis)
        
        % Calculate boundary condition
        S_star = 0.5*A_star*(1-cos(omega*Ts)); 
        Phi_star = (L(Zs(N)+dZ/2)+S_star+sqrt((L(Zs(N)+dZ/2)+S_star)^2+(M(Zs(N)+dZ/2)+L(Zs(N)+dZ/2))*(M(Zs(N)+dZ/2)-L(Zs(N)+dZ/2))))/(M(Zs(N)+dZ/2)+L(Zs(N)+dZ/2)) - 1;
       
    
        % Calculate fluxes
        
        Fs = base1;
    
        for i = 2:N
            Fs(1,i) = -A((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2) - D((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2)*(Phis(i)-Phis(i-1))/dZ;
        end
        
        % Boundary conditions
        % Phi(1,t) = Phi_star
        FU = -A(Phi_star,Zs(N)+dZ/2) - D(Phi_star,Zs(N)+dZ/2)*2*(Phi_star-Phis(N))/dZ; % = -adot (Vs=adot F(N+1)=-Vs but would need to calculate adot...)
        % Vf(0,t) = Vs(0,t) = 0
        FD = 0;

        Fs(1,1) = FD;
        Fs(1,N+1) = FU;
    
        Phidots = -(Fs(2:N+1) - Fs(1:N))/dZ;
        Phidots = Phidots';
    
    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % Function for displacement, pressure, flux, stress, strain
    function [U_wall,Uss,Ps,Qs,Ss,dUdZ] = main_ode_post(Ts,Phis)


        % Empty structures for saving
        U_wall = base3; % displacement at cell walls
        Uss = base2; % displacement at cell centres
        Ps = base2; % pressure at cell centres
        Qs = base2; % flux at cell centres
        Ss = base2; % stress at cell centres
            
        for it = 1:length(Ts)

            Phi = Phis(it,:);
            T = Ts(it);

            S_star = params.S_star(T); 
            
            % Stress
            S = stress(Phi,Zs);
            
            % Displacement
                % dUs/dZ = Phi with Us(0,t) = 0 or Us(Z,t) = -int_0^Z{Phi dx}s
                % Integrate to find the integral values at the i+1/2 walls
                U = [0,cumsum(dZ*(Phi(1:end)))];
                % Interpolate to find the integral values at the cell centres i
                Us = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2,U,Zs);
            
            % Flux
                % Calculate -dP/dZ at i+1/2 walls 
                dP_inner_walls = -(S(2:N)-S(1:N-1))/dZ;
                % Add boundary values. Ghost points P_0 = P_1 ; P_N+1 = 2*Phi_star - P_N
                dP_walls = [0,dP_inner_walls,-(S_star-S(end))/dZ/2];
                % Interpolate to find dP/dZ at cell centres i
                dP = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2, dP_walls, Zs);
                Q = k(Phi,Zs).*dP./(1+Phi);

            % Pressure
            P = S - S_star;

            % Save
            Uss(it,:) = Us;
            U_wall(it,:) = U;
            Ss(it,:) = S;
            Qs(it,:) = Q;
            Ps(it,:) = P;
            dUdZ = Phis;

        end

    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % Function for flux quantities
    function [fluxes] = calc_flux_quantities(Qs,p)

        Flux_cum = zeros(size(Phis));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cumulative Flux
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:length(Zs)
            Flux_cum(:,m) = cumtrapz(Ts(:),Qs(:,m));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Flux over one period
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,index_2]=min(abs(Ts-2*pi/omega)); % find index T=2*pi/omega (first period)
        [~,index_p1]=min(abs(Ts-2*(p-2)*pi/omega)); % find index T=2*(p-1)*pi/omega
        [~,index_p2]=min(abs(Ts-2*p*pi/omega)); % find index T=2p*pi/omega (penultimate period)
        Net_flux_1 = sum(dZ*trapz(Ts(1:index_2),Qs(1:index_2,:)),2); % net flux first period
        Net_flux_Z = trapz(Ts(index_p1:index_p2),Qs(index_p1:index_p2,:)); % net flux penultimate period for each Z
        Net_flux_p = sum(dZ*Net_flux_Z,2); % net flux penultimate period overall

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Flux over all simulated time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Net_flux_total = trapz(Ts,Qs); % for each z
        Net_flux_tot = sum(dZ*Net_flux_total,2); % overall
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fluxes.flux_cum = Flux_cum;
        fluxes.net_flux_Z = Net_flux_Z;
        fluxes.net_flux_1 = Net_flux_1;
        fluxes.net_flux_p = Net_flux_p;
        fluxes.net_flux_tot = Net_flux_tot;

    end

    % Function for advection and diffusion coefficients
    function [Adv,Diff] = calc_adv_diff_coeffs(Phis,Zs)
                Adv = A(Phis,Zs);
                Diff = D(Phis,Zs);
    end

             

end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Function for progress
function status = ode_check_save(T,Phis,flag,params)

    global count;
    
    if strcmp(flag,'done')
        disp([' - Done.'])
    else
        t = T(1);
        count = count + 1;
        if mod(count,30)==0
            disp([' - Solving at time t = ' num2str(t) ' (' num2str(100*t/max(params.tspan)) '%) ...'])
        end
    end
    
    status = 0; % status = 0 --> Continue running
    
    % N = params.N; 
        
    %     status = 0; % status = 0 --> Continue running
        
    %     % Extract the current solution
        
    %     if strcmp(flag,'init')
            
    %     end
    
    %     % Basic error check
    %     ns = Y(1:N,1);
    %     a = params.a(t);
    %     Xs = params.Xs;
    %     xs = (1-a)*Xs + a;
    %     us = u_from_n(ns,a,params);
    %     u_err = us(end)/(1-a);
    %     if abs(u_err)>1E-5
    %         disp(['*** WARNING: rel err in displacement at right boundary = ' num2str(u_err)])
    %     end
    
    % end

end
