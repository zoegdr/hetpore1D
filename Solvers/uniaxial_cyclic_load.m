function [params,Ts,Zs,Phis,ks,Ss] = uniaxial_cyclic_load(params)

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
    %           - dec_neg/inc_neg: decrease/increase in stiffness from Z=1 to Z=0
    %           - dec_pos/inc_pos: decrease/increase in stiffness from Z=0 to Z=1
    %       *.lp: location of permeability damage
    %       *.dp: magnitude of permeability damage
    %       *.ls: location of stiffness damage
    %       *.ds: magnitude of stiffness damage
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Set up -------------------------------------------------------------
    
     % Input parameters

    if isstruct(params) == 0
        clear params
        params = default_parameters('stress');
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
                k = @(Phi,Z) (1-dp)*(Phi/Phi_0+1).^3./(1+Phi);
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
    elseif strcmp(params.stress_law,'uncrimp')
        stress = @(Phi,Z) uncrimp_stress(Phi,Z,params);
        Phi_from_s = @(s,Z) uncrimp_inverse(s,Z,params);  
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

    % Compute permeability and stress
    ks = k(Phis,Zs);
    Ss = zeros(size(Phis));
    for j = 1:length(Ts)
        Phi = Phis(j,:);
        Ss(j,:) = stress(Phi,Zs);
    end


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

    % Function for uncrimping stress law and inverse
    function stress = uncrimp_stress(Phi,Z,params)

        M = params.M;
        L = params.L;
        c = params.c;
        theta0 = params.theta0;
        E = params.E; % fibril Young modulus
        J_t = sec(theta0); % threshold J beyond which all fibrils are taut

        J = 1+Phi;
        
        if J <= 1
            stress = (1-c)*( 0.5*M(Z) .*( J -1./(J) ) + 0.5*L(Z) .*( J +1./(J) -2 ) );
        elseif (J>1) && (J<J_t)
            stress = (1-c)*( 0.5*M(Z) .*( J -1./(J) ) + 0.5*L(Z) .*( J +1./(J) -2 ) ) + ...
                c*E/3/(sin(theta0)^2)*( 2*J -3 + 1./J.^2 );
        else
            stress = (1-c)*( 0.5*M(Z) .*( J -1./(J) ) + 0.5*L(Z) .*( J +1./(J) -2 ) ) + ...
                c*E*( 2*J*(1-cos(theta0)^3)/3/sin(theta0)^2 -1 ) ;
        end
    end
    
    function Phi = uncrimp_inverse(s,Z,params)

        M = params.M;
        L = params.L;
        c = params.c;
        theta0 = params.theta0;
        E = params.E; % fibril Young modulus
        s_t = uncrimp_stress(sec(theta0)-1,0,params); % threshold s beyond which all fibrils are taut

        s_eff = s/(1-c);
        beta = 2*(1-cos(theta0)^3)/3/sin(theta0)^2;
        
        if s <= 0
            A = (M(Z)+L(Z))/2;
            B = - L(Z) - s_eff;
            C = (L(Z)-M(Z))/2;
            Phi = (-B+sqrt(B^2-4*A*C))/2/A -1;
        elseif (s>0) && (s<s_t)
            fun = @(J) J.^3*(2*c*E/3/(sin(theta0))^2 + (1-c)*(M(Z)+L(Z))/2) - J.^2*(c*E/(sin(theta0))^2 + L(Z)*(1-c) + s) + J*(1-c)*(L(Z)-M(Z))/2 + c*E/3/sin(theta0)^2;
            sol = fzero(fun,[1,2]);
            Phi = sol-1;
        else
            A = c/(1-c)*E*beta + (M(Z)+L(Z))/2;
            B = -c/(1-c)*E - L(Z) - s_eff;
            C = (L(Z)-M(Z))/2;
            Phi = (-B+sqrt(B^2-4*A*C))/2/A -1;
        end
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