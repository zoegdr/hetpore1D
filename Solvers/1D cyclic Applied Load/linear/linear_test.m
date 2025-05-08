function [params,t,zs,z,phis,ws,wss,ps,qs,sigs,dwdz,fluxes] = linear_test(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up as in Scott, for tendon.
    % Solves for normalised porosity.
    
    % Poroelastic material subject to a cyclic load s_star = A(1-cos(omega*t))
    % at z=1, and no flux at z=0. 
    % Uniaxial deformation.
    % Lagrangian/Eulerian framework (same in linear (?))
    % Heterogenous stiffness M(z) and permeability k(z)
    % Linear (Hookean) constitutive law
    
    % Non-linear diffusion equation of the form
    % dphi/dt - d/dz(A(z) + k(z)*M(z)dphi/dz) = 0 
    % where A(z) = k(z)*dM(z)*phi
    
    % Boundary conditions
    % s=s_star(t), p=0 at z=0
    % ds/dz=0, w=0 at z=1
    
    % phi is normalised porosity in Lagrangian
    % w is displacement in z-axis in Lagrangian
    % s is Cauchy stress tensor
    % dp is magnitude of permeability damage
    % lp is location of permeability damage
    % ds is magnitude of stiffness damages
    % ls is location of stiffness damage
    % A = magnitude of stress as boundary
    % omega = frequency of loading
    
    % Finite volume method in space
    % Runge-Kutta in time (ode15s).
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
    phi_0 = params.Phi0; % initial uniform porosity
    nu = params.nu; % Poisson ratio
    A_star = params.Astar;
    omega = params.omega;

    % Unpack control parameters
    N = params.N;
    params.tspan = [0 params.p*2*pi/omega];
    tspan = params.tspan;
    
    % Grid
    a = 0;
    b = 1;
    dz = (b-a)/N;
    zs = a+dz/2:dz:b-dz/2; % space
    base1 = zeros(1,N+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Lame's second parameter
Lambda = @(z) nu/(1-nu);


%% Solve in time ------------------------------------------------------

global count;
count = 0;

options = odeset('RelTol',1E-5, 'AbsTol',1E-5);

% options = odeset('OutputFcn',@(t,phis,flag) ode_check_save(t,phis,flag,params),...
%                      'RelTol',1E-5, ...
%                      'AbsTol',1E-5);

% Initial condition: normalised initial porosity is just zero
phi_0s = zeros(1,N);

[t,phis] = ode15s(@odefun,tspan,phi_0s',options);

base2 = zeros(size(phis));
base3 = zeros(length(t),N+1);

% Calculate other quantities
[z,ws,wss,ps,qs,sigs,dwdz] = main_ode_post(t,phis);
[fluxes] = calc_flux_quantities(qs,params.p);


%% ODEs for Phi -------------------------------------------------------
function phidots = odefun(t,phis)
    
      % Calculate stress
        sig_star = 0.5*A_star*(1-cos(omega*t));
        phi_star = sig_star; %sig_star/stiffness
           
        % Calculate fluxes
        
        F_left = base1;
        F_right = base1;

        for i = 2:N
            F_left(1,i) = -(phis(i)-phis(i-1))/dz;
            F_right(1,i) = -(phis(i)-phis(i-1))/dz;
        end
        
        % Boundary conditions
        % phi(b,t) = phi_star
        % dphi(a,t)/dz = 0
        % Ghost points phi_0 = phi_1 ; phi_N+1 = 2*phi_star - phi_N
        F_left(1,1) = -2*(phis(1)-phi_star)/dz;
        F_right(1,N+1) = 0;

        phidots = -(F_right(2:N+1) - F_left(1:N))/dz;
        phidots = phidots';

end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% Function for displacement, pressure, flux, stress, strain
function [z,ws,wss,p,q,sig,dwdz] = main_ode_post(t,phis)

    z = a:dz:b; % z at cell walls

    % Empty structures for saving
    ws = base3; % displacement at cell walls
    wss = base2; % displacement at cell centres
    p = base2; % pressure at cell centres
    q = base2; % flux at cell centres
    sig = base2; % stress at cell centres
    dwdz = base2; % strain at cell centres
        
    for n = 1:length(t)

         sig_star(n) = 0.5*A_star*(1-cos(omega*t(n)));
         phi_star(n) = sig_star(n);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISPLACEMENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dWs/dz = Phi with ws(a,t) = 0 or ws(z,t) = -int_a^z{Phi dx}s

        % Integrate to find the integral values at the i+1/2 walls
        ws(n,:) = [0,cumsum(dz*(phis(n,1:end)))];
        % Interpolate to find the integral values at the cell centres i
        wss(n,:) = interp1(zs(1)-dz/2:dz:zs(end)+dz/2, ws(n,:), zs);
        wss(n,:) = flip(wss(n,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLUX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calculate -dp/dz at i+1/2 walls 
        for j = 1:N-1
            dp_inner_walls(n,j) = -(phis(n,j+1)-phis(n,j))/dz;
        end
        % Add boundary values. Ghost points P_0 = P_1 ; P_N+1 = 2*phi_star - P_N
        dp_walls(n,:) = [ -2*( phis(n,1)-phi_star(n) )/dz , dp_inner_walls(n,:) , 0 ];
        % Interpolate to find dp/dz at cell centres i
        dp_(n,:) = interp1(zs(1)-dz/2:dz:zs(end)+dz/2, dp_walls(n,:), zs);
    

        for s = 1:N

        % Calculate flux = -k(Phi)/J*dp/dz
        q(n,s) = dp_(n,s)/(1+phis(n,s));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PRESSURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % dp/dz = d/dz(s), p(1)=0 
        % => p(z) = s(z)-S_star
        p(n,s) = phis(n,s)-phi_star(n);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STRESS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
            sig(n,s) = phis(n,s);
            
        end

    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STRAIN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dwdz = phis;

end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

% Function for flux quantities
function [fluxes] = calc_flux_quantities(q,p)

    flux_cum = zeros(size(phis));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cumulative Flux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:length(zs)
        flux_cum(:,m) = cumtrapz(t(:),q(:,m));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Net Flux over one period
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,index_2]=min(abs(t-2*pi/omega)); % find index T=2*pi/omega (first period)
    [~,index_p1]=min(abs(t-2*(p-2)*pi/omega)); % find index T=2*(p-1)*pi/omega
    [~,index_p2]=min(abs(t-2*p*pi/omega)); % find index T=2p*pi/omega (penultimate period)
    net_flux_1 = sum(dz*trapz(t(1:index_2),q(1:index_2,:)),2); % net flux first period
    net_flux_Z = trapz(t(index_p1:index_p2),q(index_p1:index_p2,:)); % net flux penultimate period for each z
    net_flux_p = sum(dz*net_flux_Z,2); % net flux penultimate period  overall

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Net Flux over all simulated time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    net_flux_total = trapz(t,q); % for each z
    net_flux_tot = sum(dz*net_flux_total,2); % overall

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Net Inflow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Net_inflow = trapz(t,q(:,1));
    Net_inflow_p = trapz(t(index_p1:index_p2),q(index_p1:index_p2,1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fluxes.flux_cum = flux_cum;
    fluxes.net_flux_1 = net_flux_1;
    fluxes.net_flux_Z = net_flux_Z;
    fluxes.net_flux_p = net_flux_p;
    fluxes.net_flux_tot = net_flux_tot;
    fluxes.net_inflow = Net_inflow;
    fluxes.Net_inflow_p = Net_inflow_p;


end
         

end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

% Function for progress
function status = ode_check_save(t,phis,flag,params)

global count;

if strcmp(flag,'done')
    disp([' - Done.'])
else
    t = t(1);
    count = count + 1;
    if mod(count,10)==0
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