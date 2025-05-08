function [params,T,Zs,Z,Phis,Ws,Wss,P,Q,Ss,dWdZ,Fluxes] = tendon_uniaxial_cyclic_load_stressdiff_COPY(params)

% WORKING 11/07/2024 13:47

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up as in Scott, for tendon.
    % Solves for normalised porosity.
    
    % Poroelastic material subject to a cyclic load s_star = A(1-cos(omega*t))
    % at Z=1, and no flux at Z=0. 
    % Uniaxial deformation.
    % Lagrangian framework
    % Heterogenous stiffness M(Z) and permeability k(Phi,Z)
    % Non-linear constitutive law
    
    % Non-linear diffusion equation of the form
    % dPhi/dt - d/dZ (A(Phi) + D(Phi)*ds/dZ) = 0
    % with relevant expression for A(Phi), D(Phi) depending on
    % elasticity and permeability law
    
    % Boundary conditions
    % s=s_star, P=0 at Z=1
    % ds/dZ=0, W=0 at Z=0
    
    % Phi is normalised porosity in Lagrangian
    % W is displacement in Z-axis in Lagrangian
    % s is Piola-Kirchoff stress tensor
    % dp is magnitude of permeability damage
    % lp is location of permeability damage
    % ds is magnitude of stiffness damages
    % ls is location of stiffness damage
    % A_star = magnitude of stress as boundary
    % omega = frequency of loading
    
    % Finite volume method in space
    % Runge-Kutta in time (ode15s).
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
    params.S_star = @(t) (A_star/2)*(1-cos(omega*t));

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
            k = @(Phi,Z) ones(size(Phi))* (1-dp*exp(-(Z-lp)^2/(2*(1/16)^2))); % Constant
        elseif strcmp(params.perm_law,'KC')
            k = @(Phi,Z) (Phi/Phi_0+1).^3./(1+Phi) * (1-dp*exp(-(Z-lp)^2/(2*(1/16)^2))); % Normalized Kozeny-Carman, simplified expression for k(Phi/J) where Phi is normalised porosity
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
        stress = @(Phi,Z) M(Z).*Phi; % Linear elasticity (pretty sure this is the expression in Lagrangian for uniaxial?)
        Phi_from_s = @(s,Z) s./M(Z);
    elseif strcmp(params.stress_law,'log')
        stress = @(Phi,Z) M(Z).*log(1+Phi)./(1+Phi);
        Phi_from_s = @(s,Z) -1-lambertw(-s./M(Z))/(s/M(Z));
    elseif strcmp(params.stress_law,'neo')
        stress = @(Phi,Z) 0.5*M(Z) .*( 1+ Phi-1./(1+Phi) ) + 0.5*L(Z) .*( 1 + Phi+1./(1+Phi) -2 ); % Neo-Hookean elasticity
        Phi_from_s = @(s,Z) ( L(Z)+ s + sqrt((L(Z)+s).^2 + ( M(Z)+L(Z)).*(M(Z)-L(Z))) ) ./ ( M(Z)+L(Z) ) -1;
    
    else
        error('Unknown stress law.')
    end

    % Define a function for diffusivity
    D = @(Phi,Z) k(Phi,Z)./(1+Phi);

    %% Solve in time ------------------------------------------------------
    
    global count;
    count = 0;
    
   options = odeset('RelTol',1E-5, 'AbsTol',1E-5);

    % options = odeset('OutputFcn',@(T,Phis,flag) ode_check_save(T,Phis,flag,params),...
                         % 'RelTol',1E-5, ...
                         % 'AbsTol',1E-5);

    % Initial condition: normalised initial porosity is just zero
    Phi_0s = zeros(1,N);
    
    [T,Phis] = ode15s(@odefun,tspan,Phi_0s',options);

    base2 = zeros(size(Phis));
    base3 = zeros(length(T),N+1);

    % Calculate other quantities
    [Z,Ws,Wss,P,Q,Ss,dWdZ] = main_ode_post(T,Phis);
    [Fluxes] = calc_flux_quantities(Q,params.p);


    %% ODEs for Phi -------------------------------------------------------
    function Phidots = odefun(T,Phis)
        
        % Calculate boundary condition
        S_star = params.S_star(T); 
        Phi_star = Phi_from_s(S_star,Zs(end)+dZ/2);

        S = stress(Phis,Zs);
       
        % Calculate f.v. fluxes
        
        F_left = base1;
        F_right = base1;
    
        for i = 2:N
            F_left(1,i) = - D((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2)*(S(i)-S(i-1))/dZ;
            F_right(1,i) = - D((Phis(i)+Phis(i-1))/2,Zs(i)-dZ/2)*(S(i)-S(i-1))/dZ;
        end
        
        % Boundary conditions
        % Phi(b,t) = Phi_star
        % dPhi(a,t)/dZ = 0
        % Ghost point P_-1/2 = P_1/2 ;
        % F_right(N) = -D(Phi_star)*(Phi_star-Phi(N-1/2))/dZ/2  % one sided difference
        F_left(1,1) = 0; %(Vf = 0 --> F(1)=0 always) %=-A(Phis(1),Zs(1)-dZ/2);
        F_right(1,N+1) = - D(Phi_star,Zs(N)+dZ/2)*2*(S_star-S(N))/dZ; % = -adot (Vs=adot F(N+1)=-Vs but would need to calculate adot...)
    
        Phidots = -(F_right(2:N+1) - F_left(1:N))/dZ;
        Phidots = Phidots';
    
    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % Function for displacement, pressure, flux, stress, strain
    function [Z,Ws,Wss,P,Q,S,dWdZ] = main_ode_post(T,Phis)

        Z = a:dZ:b; % Z at cell walls

        % Empty structures for saving
        Ws = base3; % displacement at cell walls
        Wss = base2; % displacement at cell centres
        P = base2; % pressure at cell centres
        Q = base2; % flux at cell centres
        S = base2; % stress at cell centres
        dWdZ = base2; % strain at cell centres
            
        for n = 1:length(T)

            S_star(n) = 0.5*A_star*(1-cos(omega*T(n))); 
            Phi_star(n) = (L(Zs(N)+dZ/2)+S_star(n)+sqrt((L(Zs(N)+dZ/2)+S_star(n))^2+(M(Zs(N)+dZ/2)+L(Zs(N)+dZ/2))*(M(Zs(N)+dZ/2)-L(Zs(N)+dZ/2))))/(M(Zs(N)+dZ/2)+L(Zs(N)+dZ/2))-1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DISPLACEMENT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % dWs/dZ = Phi with Ws(a,t) = 0 or Ws(Z,t) = -int_a^Z{Phi dx}s

            % Integrate to find the integral values at the i+1/2 walls
            Ws(n,:) = [0,cumsum(dZ*(Phis(n,1:end)))];
            % Interpolate to find the integral values at the cell centres i
            Wss(n,:) = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2, Ws(n,:), Zs);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FLUX
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate -dP/dZ at i+1/2 walls 
            for j = 1:N-1
                dP_inner_walls(n,j) = - 0.5* (((L(Zs(j+1))+M(Zs(j+1)))*(1+Phis(n,j+1))+(L(Zs(j+1))-M(Zs(j+1)))./(1+Phis(n,j+1))-2*L(Zs(j+1)))...
                -((L(Zs(j))+M(Zs(j)))*(1+Phis(n,j))+(L(Zs(j))-M(Zs(j)))./(1+Phis(n,j))-2*L(Zs(j))))/dZ;
            end
            % Add boundary values. Ghost points P_0 = P_1 ; P_N+1 = 2*Phi_star - P_N
            dP_walls(n,:) = [0,dP_inner_walls(n,:),...
            -0.5*(((L(Zs(N)+dZ/2)+M(Zs(N)+dZ/2))*(1+2*Phi_star(n)-Phis(n,N))+(L(Zs(N)+dZ/2)-M(Zs(N)+dZ/2))./(1+2*Phi_star(n)-Phis(n,N))-2*L(Zs(N)+dZ/2))...
            - ((L(Zs(N))+M(Zs(N)))*(1+Phis(n,N))+(L(Zs(N))-M(Zs(N)))./(1+Phis(n,N))-2*L(Zs(N))))/dZ];
            % Interpolate to find dP/dZ at cell centres i
            dP(n,:) = interp1(Zs(1)-dZ/2:dZ:Zs(end)+dZ/2, dP_walls(n,:), Zs);
        

            for s = 1:N

            % Calculate flux = -k(Phi)/J*dP/dZ
            Q(n,s) = k(Phis(n,s),Zs(s))*dP(n,s)/(1+Phis(n,s));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRESSURE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % dP/dZ = d/dZ(S), P(1)=0 
            % => P(Z) = S(Z)-S_star
            P(n,s) = 0.5*((L(Zs(s))+M(Zs(s)))*(1+Phis(n,s))+(L(Zs(s))-M(Zs(s)))./(1+Phis(n,s))-2*L(Zs(s)))...
                - 0.5*((L(Zs(end)+dZ/2)+M(Zs(end)+dZ/2))*(1+Phi_star(n))+(L(Zs(end)+dZ/2)-M(Zs(end)+dZ/2))./(1+Phi_star(n))-2*L(Zs(end)+dZ/2));

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STRESS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
                S(n,s) = .5*M(Zs(s)).*(1+Phis(n,s)-1./(1+Phis(n,s)))+.5*L(Zs(s)).*(Phis(n,s)+1./(Phis(n,s)+1)-2);
                
            end

        end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STRAIN
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dWdZ = Phis;

    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % Function for flux quantities
    function [fluxes] = calc_flux_quantities(Q,p)

        Flux_cum = zeros(size(Phis));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cumulative Flux
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:length(Zs)
            Flux_cum(:,m) = cumtrapz(T(:),Q(:,m));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Flux over one period
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,index_2]=min(abs(T-2*pi/omega)); % find index T=2*pi/omega (first period)
        [~,index_p1]=min(abs(T-2*(p-2)*pi/omega)); % find index T=2*(p-1)*pi/omega
        [~,index_p2]=min(abs(T-2*p*pi/omega)); % find index T=2p*pi/omega (penultimate period)
        Net_flux_1 = sum(dZ*trapz(T(1:index_2),Q(1:index_2,:)),2); % net flux first period
        Net_flux_Z = trapz(T(index_p1:index_p2),Q(index_p1:index_p2,:)); % net flux penultimate period for each Z
        Net_flux_p = sum(dZ*Net_flux_Z,2); % net flux penultimate period overall

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Flux over all simulated time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Net_flux_total = trapz(T,Q); % for each z
        Net_flux_tot = sum(dZ*Net_flux_total,2); % overall
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fluxes.flux_cum = Flux_cum;
        fluxes.net_flux_Z = Net_flux_Z;
        fluxes.net_flux_1 = Net_flux_1;
        fluxes.net_flux_p = Net_flux_p;
        fluxes.net_flux_tot = Net_flux_tot;

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
