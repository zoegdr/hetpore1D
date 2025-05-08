function [xss,ts,fss,uss,pss,sigss,qs,dps,as,F_lag,X_lag] = consol_uniaxial_eul(params)

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
    % INPUTS
    %   params: Structure with parameters and options.
    %     *.q: Dimensionless flow rate
    %     *.dp: Dimensionless pressure drop
    %       - Exactly one of the two above must be passed in. The other should be absent or NaN.
    %     *.sigstar: Dimensionless effective stress applied at x=a(t)
    %     *.f0: Initial (relaxed) porosity
    %     *.p: Times at which to save and return the solution.
    %
    % OPTIONAL INPUTS
    %   params
    %     *.stress_law: 
    %       - linear: Linear elasticity [default]
    %       - log: Hencky elasticity
    %     *.perm_law:
    %       - const: Constant permeability, where k is k/k0 = 1 [default]
    %       - KC: normalized Kozeny-Carman, where k is k/k0 = (((1-f0)^2)/(f0^3))*((f^3)/((1-f)^2))
    
    % Unpack physical parameters
    stress_law = params.stress_law;
    perm_law = params.perm_law;
    
    f0 = params.f0; % [-] initial porosity
    if ~isfield(params,'q')
        params.q = NaN;
    end
    q = params.q;
    if ~isfield(params,'dp')
        params.dp = NaN;
    end
    dp = params.dp;
    if ~isfield(params,'sigstar')
        params.sigstar = 0;
    end
    if isnan(q) & ~isnan(dp)
        bc_fluid = 'dp-fixed';
    elseif ~isnan(q) & isnan(dp)
        bc_fluid = 'q-fixed';
    else
        error('Must input exactly ONE of q and dp. The other will be determined by the solution.')
    end
    sigstar = params.sigstar;
    M0 = 1;
    nu = params.nu;
    L0 = M0*nu/(1-nu);
        
    % Unpack control parameters
    N = params.N;
    params.tspan = [0 params.p];
    tspan = params.tspan;
        
    a0 = 0;
	b0 = 1;
    
	% Initial condition: Start in a relaxed state
    a = a0;
    b = b0;
    dx = (b-a)/N;
    xs = [(a+dx/2):dx:(b-dx/2)];
	f0s = f0*ones(1,N);
    
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
    
    % Inner boundary condition: Imposed effective stress
    if sigstar==0
        fstar = f0; % Relaxed. Below should give same result, but this is exact.
    else
        fstar = f_from_sig(sigstar);
    end
    kstar = k(fstar);
    
    
    % Options for ODE solver
    options = odeset('RelTol',1E-10, ...
                     'AbsTol',1E-10);
    
    % Initial condition
    Y0 = [(fs');a];
    
    [T,Y] = ode15s(@odefun,tspan,Y0,options);
    
    ts = T;
    fss = Y(:,1:N);
    as = Y(:,N+1);
    
    [xss,uss,sigss,pss,kss,qs,dps] = ode_post(fss,as);

    % Convert to Lagrangian to compare with Lagrangian formulation (when
    % plotting plot against X=x-u)
    F_lag = (1-f0)./(1-fss).*fss;
    X_lag = xss - uss;
        
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [fdots,adot] = calcFF(fs,a)
        
        % Calculate grid
        dx = (b-a)/N;
        xs = [(a+dx/2):dx:(b-dx/2)];
        
        % Calculate stress field from porosity field
        sigs = sig(fs);
        
        % Calculate permeability field from porosity field
        ks = k(fs);
        kstar = k(fstar);
        
        if strcmp(bc_fluid,'dp-fixed')
            q = calc_q(sigs);
        end
        
        adot = q + kstar*(sigs(1)-sigstar)/(dx/2);
        
        % Finite Volume Fluxes
        Fs = zeros(1,N+1);
        for i = 2:N
            Fs(1,i) = ( -(b-(xs(i)+dx/2))*adot/(b-a) + q )*(fs(i)+fs(i-1))/2 ...
                - ( (1-fs(i-1)).*ks(i-1) + (1-fs(i)).*ks(i) )/2 * ( (sigs(i)-sigs(i-1))/dx );
        end

        % Boundary conditions
        FL = q - adot;
        FR = q;

        Fs(1,1) = FL;
        Fs(1,N+1) = FR;

        fdots = (adot/(b-a))*fs - (Fs(1,2:end)-Fs(1,1:end-1))/dx;
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [q] = calc_q(sigs)

            % Calculate q (needed to evolve f)
            sig_b = sigstar - dp; % p(b)-p(a)=sig(b)-sig(a), dp=p(a)-p(b), sig(a)=sigstar --> sig_b = sigstar-dp
            f_b = f_from_sig(sig_b);
            k_b = k(f_b);
            q = -k_b*(sig_b-sigs(end))/(dx/2);
     
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
    function Ydot = odefun(T,Y)   
        
        fs = Y(1:N)';
        a = Y(N+1);
        [fdots,adot] = calcFF(fs,a);
        Ydot = [(fdots');adot];
        
    end

    % -----------------------------------------------------
    % -----------------------------------------------------
    function [xss,uss,sigss,pss,kss,qs,dps] = ode_post(fss,as)
        
        % Empty structures for saving
        xss = zeros(size(fss));
        uss = zeros(size(fss));
        sigss = zeros(size(fss));
        pss = zeros(size(fss));
        kss = zeros(size(fss));
        qs = zeros(size(as));
        dps = zeros(size(as));
        
        % Loop over time
        for it=1:size(fss,1)
            
            fs = fss(it,:);
            a = as(it);
            
            % Calculate grid
            dx = (b-a)/N;
            xs = [(a+dx/2):dx:(b-dx/2)];
            
            % Calculate displacement field from porosity field
            [a_from_u,us] = u_from_f(xs,fs);
            
            % Calculate stress field from porosity field
            sigs = sig(fs);
            
            % Calculate permeability field from porosity field
            ks = k(fs);
            
            q = calc_q(sigs);
            
            % Save
            xss(it,:) = xs;
            uss(it,:) = us;
            sigss(it,:) = sigs;
            pss(it,:) = (sigs-sigstar) + dp; % dp/dx = dsig/dx, dp = p(a)-p(b)
            kss(it,:) = ks;
            qs(it) = q;
            dps(it) = dp;
            
        end
            
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    
end