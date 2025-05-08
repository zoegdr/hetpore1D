function [params,Ts,Zs,Phis,U_wall,Uss,Ps,Qs,Ss,dUdZ,Fluxes,dPs] = intermediate_test(params)

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
    % Runge-Kutta in time (ode15s).

    % GEOMETRY
    %   This code takes the upper and lower edges of the material to be at Lagrangian coordinates 
    %   Z=0 and Z=1 respectively, where the upper edge is free and the lower edge is fixed. 
    %   The upper edge Z=0 is subject to an applied cyclic load S_star
    %   The lower edge Z=1 is impermeable Vs=Vf=Us=0
    %
    % METHOD
    %   * This code solves a non-linear advection-diffusion equation of the form
    %     dPhi/dt - d/dZ (k(Phi)/(1+Phi)*dS/dZ) = 0
    %     where S(Phi) is calculated from stress law
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
        params = create_parameters;
    % else
    %     question = input("Do you want to change any parameters? ");
    %     if strcmp(question,'yes')
    %         params = create_parameters;
        % else
            %
        % end
    end

    % ---------------------------------------------------------------------

    % Unpack physical parameters
    Phi_0 = params.Phi0; % initial uniform porosity
    nu = params.nu; % Poisson ratio
    A_star = params.Astar;
    omega = params.omega;
    params.S_star = @(t) (A_star/2)*(1-cos(omega*t));

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
            k = @(Phi) ones(size(Phi)); 
        elseif strcmp(params.perm_law,'KC')
            k = @(Phi) (Phi/Phi_0+1).^3./(1+Phi); % Normalized Kozeny-Carman, simplified expression for k(Phi/J) where Phi is normalised porosity
        else
            error('Unknown permeability law.')
        end
        params.k = k;
    

    params.M = 1;

    
    % Lame's second parameter
    params.L =  nu/(1-nu) * params.M;

    M = params.M;
    L = params.L;


    % Define a function stress(f) for the stress law, and its inverse f_from_sig(sig).
    if strcmp(params.stress_law,'linear')
        stress = @(Phi) M.*Phi; 
        Phi_from_s = @(s) s./M;
    elseif strcmp(params.stress_law,'log')
        stress = @(Phi) M.*log(1+Phi)./(1+Phi);
        Phi_from_s = @(s) -1-lambertw(-s./M)/(s/M);
    elseif strcmp(params.stress_law,'neo')
        stress = @(Phi) 0.5*M .*( 1+ Phi-1./(1+Phi) ) + 0.5*L .*( 1 + Phi+1./(1+Phi) -2 ); % Neo-Hookean elasticity
        Phi_from_s = @(s) ( L+ s + sqrt((L+s).^2 + ( M+L).*(M-L)) ) ./ ( M+L ) -1;
    
    else
        error('Unknown stress law.')
    end

    % Define a function for diffusivity
    D = @(Phi) k(Phi)./(1+Phi);

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
    [U_wall,Uss,Ps,Qs,Ss,dUdZ,dPs] = main_ode_post(Ts,Phis);
    [Fluxes] = calc_flux_quantities(Qs,params.p);


    %% ODEs for Phi -------------------------------------------------------
    function Phidots = odefun(Ts,Phis)
        
        % Calculate boundary condition
        S_star = params.S_star(Ts); 
        % Phi_star = Phi_from_s(S_star,Zs(end)+dZ/2);
        
        Phi_star = Phi_from_s(S_star);


        % S = stress(Phis,Zs);
       
        % Calculate f.v. fluxes
        
        Fs = base1;
    
        for i = 2:N
            Fs(1,i) = - D((Phis(i)+Phis(i-1))/2)*(stress(Phis(i))-stress(Phis(i-1)))/dZ;
        end
        
        % Boundary conditions // applied load at Z=0, no flux at Z=1
            % Phi(1,t) = Phi_star
            FU = - D(Phi_star)*2*(stress(Phis(1))-S_star)/dZ;
            % Vf(0,t) = Vs(0,t) = 0
            FD = 0; %0;
        
        Fs(1,1) = FU; 
        Fs(1,end) = FD;
    
        Phidots = -(Fs(2:N+1) - Fs(1:N))/dZ;
        Phidots = Phidots';
    
    end

    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % Function for displacement, pressure, flux, stress, strain
    function [U_wall,Uss,Ps,Qs,Ss,dUdZ,dPs] = main_ode_post(Ts,Phis)


        % Empty structures for saving
        U_wall = base3; % displacement at cell walls
        Uss = base2; % displacement at cell centres
        Ps = base2; % pressure at cell centres
        Qs = base2; % flux at cell centres
        Ss = base2; % stress at cell centres
        dPs = base2;
            
        for it = 1:length(Ts)

            Phi = Phis(it,:);
            T = Ts(it);

            S_star = params.S_star(T); 
            Phi_star = Phi_from_s(S_star);
            
            % Stress
            S = stress(Phi);
            
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
                Q = -k(Phi).*dP./(1+Phi);


            % Save
            Uss(it,:) = Us;
            U_wall(it,:) = U;
            Ss(it,:) = S;
            Qs(it,:) = Q;
            Ps(it,:) = P;
            dPs(it,:) = dP;

        end
            dUdZ = Phis;
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
        Net_flux_pZ = trapz(Ts(index_p1:index_p2),Qs(index_p1:index_p2,:)); % net flux penultimate period for each Z
        Net_flux_p = sum(dZ*Net_flux_pZ,2); % net flux penultimate period overall

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Flux over all simulated time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Net_flux_total = trapz(Ts,Qs); % for each z
        Net_flux_tot = sum(dZ*Net_flux_total,2); % overall

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Net Inflow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Net_inflow = trapz(Ts,Qs(:,1));
        Net_inflow_p = trapz(Ts(index_p1:index_p2),Qs(index_p1:index_p2,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fluxes.flux_cum = Flux_cum;
        fluxes.net_flux_pZ = Net_flux_pZ;
        fluxes.net_flux_1 = Net_flux_1;
        fluxes.net_flux_p = Net_flux_p;
        fluxes.net_flux_tot = Net_flux_tot;
        fluxes.net_inflow = Net_inflow;
        fluxes.Net_inflow_p = Net_inflow_p;


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
    


end
