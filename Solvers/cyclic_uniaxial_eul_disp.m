function [xss,ts,fss,uss,pss,sigss,F_lag,X_lag] = cyclic_uniaxial_eul_disp(params)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UNIAXIAL EULERIAN FRAMEWORK
    %
    % Variables and Outputs
    %   * t is time
    %   * z is space
    %   * f is porosity
    %   * u is solid displacement
    %   * p is fluid pressure
    %   * sig is Cauchy stress tensor
    %   * a is applied displacement
    %   * F_lag is Lagrangian porosity
    %   * X_lag is Lagrangian coordinate
    %
    % Finite volume method in space
    % Runge-Kutta in time (ode15s).

    % GEOMETRY
    %   This code takes the left and right edges of the sponge to be at x=a(t) and x=b=L, respectively,
    %   where the left edge is free and the right edge is fixed. The code assumes 
    %   without loss of generality that a0/L=0 and b0/L=b/L=1.
    %   The left edge is subject to an applied cyclic displacement a(t)
    %       vs(a,t) = adot
    %   There is no flux through or displacement at the right edge.
    %       vs(b,t)=vf(b,t)=us(b,t)=0
    %
    % METHOD
    %   * This code solves a non-linear advection-diffusion equation of the form
    %       df/dt - d/dz ((1-f)*k(f)*dsig/dz) = 0
    %       with a prescribed moving boundary a(t)
    %     where 
    %       sig(f) is calculated from stress law
    %       a(t) is updated at each timestep 
    %       kstar = k(fstar), fstar = f(sigstar)
    %   * Finite Volume method in space
    %   * Runge-Kutta (ode15s) in time
    %
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
    %           - KC: normalised Kozeny-Carman k/k0 = (1-f0)^2/f0^3 * f^3/(1-f)^2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
    % Unpack physical parameters
    f0 = params.f0; % [-] initial porosity
    Astar = params.Astar;
    omega = params.omega;
    perm_law = params.perm_law;
    stress_law = params.stress_law;
    M0 = params.M;
    nu = params.nu;
    L0 = M0*nu/(1-nu);
    params.a = @(t) -(Astar/2)*(1-cos(omega*t)); % applied displacement
    params.adot = @(t) -(Astar*omega/2)*sin(omega*t);
        
    % Unpack control parameters
    N = params.N;
    params.tspan = [0 params.p*2*pi/omega];
    tspan = params.tspan;
        
    a0 = 0;
	b0 = 1;
    
	% Initial condition: Start in a relaxed state
    a = a0;
    b = b0;
    dx = (b-a)/N;
    xs = (a+dx/2):dx:(b-dx/2);
	f0s = f0*ones(1,N); % uniform porosity
    
    fs = f0s;
    
    % Define a function for the permeability law
    if strcmp(perm_law,'const')
        k = @(f) ones(size(f)); % Constant
    elseif strcmp(perm_law,'KC')
        k = @(f) (((1-f0)^2)/(f0^3)) * ((f.^3)./((1-f).^2)); % Normalized Kozeny-Carman
    else
        error('Unknown permeability law.')
    end

    J = @(f) (1-f0)./(1-f);
    
    % Define a function sig(f) for the stress law, and its inverse f_from_sig(sig).
    if strcmp(stress_law,'linear')
        sig = @(f) (f-f0)/(1-f0); % Linear elasticity
        f_from_sig = @(s) f0 + (1-f0)*s;
    elseif strcmp(stress_law,'log')
        sig = @(f) log((1-f0)./(1-f))./((1-f0)./(1-f)); % Hencky elasticity
        f_from_sig = @(s) 1 - (1-f0)*((-s)./lambertw(-s)); % y=-lambertw(-x)/x solves x=ln(y)/y
    elseif strcmp(stress_law,'neo')
       sig = @(f) 0.5*M0*(J(f)-1./(J(f)))+0.5*L0*(J(f)+1./(J(f))-2); % Neo-Hookean elasticity
        f_from_sig = @(s) -(1-f0)*(M0+L0)./ ( L0+s + sqrt((L0+s).^2 + (M0+L0)*(M0-L0)) ) + 1;
    else
        error('Unknown stress law.')
    end
    
    
    % Options for ODE solver
    options = odeset('RelTol',1E-10, ...
                     'AbsTol',1E-10);

    
    [ts,fs] = ode15s(@odefun,tspan,f0s,options);

    fss = fs;
    base1 = zeros(size(fss));

    [xss,uss,sigss,pss,kss] = ode_post(ts,fss);

    % Convert to Lagrangian to compare with Lagrangian formulation (when
    % plotting plot against X=x-u)
    F_lag = (1-f0)./(1-fss).*fss;
    X_lag = xss - uss;
        
    % -----------------------------------------------------
    % -----------------------------------------------------
    function fdots = odefun(ts,fs)

        fs = fs';

        a = params.a(ts);
        adot = params.adot(ts);
        
        % Calculate grid
        dx = (b-a)/N;
        xs = (a+dx/2):dx:(b-dx/2);
        
        % Calculate stress field from porosity field
        sigs = sig(fs);
        
        % Calculate permeability field from porosity field
        ks = k(fs);
        
        % Finite Volume Fluxes
        Fs = zeros(1,N+1);
        for i = 2:N
            Fs(1,i) = ( -(b-(xs(i)+dx/2))*adot/(b-a) )*(fs(i)+fs(i-1))/2 ...
                - ( (1-fs(i-1)).*ks(i-1) + (1-fs(i)).*ks(i) )/2 * ( (sigs(i)-sigs(i-1))/dx );
        end

        % Boundary conditions
        FL = -adot;
        FR = 0;

        Fs(1,1) = FL;
        Fs(1,N+1) = FR;

        fdots = (adot/(b-a))*fs - (Fs(1,2:end)-Fs(1,1:end-1))/dx;
        fdots =  fdots' ;
        
    end

    % -----------------------------------------------------
    % -----------------------------------------------------
    function [a_from_u,us] = u_from_f(xs,fs)
        % Calculate the displacement field from the porosity
        
        % (J-1)/J = (f-f0)/(1-f0) = du/dx, where J = (1-f0)/(1-f) is the Jacobian determinant
        % Integrate to get u at the i-1/2 interfaces
        us = -( sum(dx*((fs-f0)/(1-f0))) - [0,cumsum(dx*((fs(1:end-1)-f0)/(1-f0)))] );
        % Add the right boundary value
        us = [us,0];
        
        % Save the left boundary value
        a_from_u = a0 + us(1);
        
        % Interpolate to get the cell-centered displacements
        us = interp1([xs-dx/2,xs(end)+dx/2], us, xs);
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------

    % -----------------------------------------------------
    % -----------------------------------------------------
    function [xss,uss,sigss,pss,kss] = ode_post(ts,fss)
        
        % Empty structures for saving
        xss = base1;
        uss = base1;
        sigss = base1;
        pss = base1;
        kss = base1;
        
        % Loop over time
        for it=1:size(fss,1)
            
            fs = fss(it,:);
            a = params.a(ts(it));
            
            % Calculate grid
            dx = (b-a)/N;
            xs = [(a+dx/2):dx:(b-dx/2)];
            
            % Calculate displacement field from porosity field
            [~,us] = u_from_f(xs,fs);
            
            % Calculate stress field from porosity field
            sigs = sig(fs);
            
            % Calculate permeability field from porosity field
            ks = k(fs);
            
            % Save
            xss(it,:) = xs;
            uss(it,:) = us;
            sigss(it,:) = sigs;
            pss(it,:) = sigs - sigs(1); % this is an approximation
            kss(it,:) = ks;
            
        end
            
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    
end