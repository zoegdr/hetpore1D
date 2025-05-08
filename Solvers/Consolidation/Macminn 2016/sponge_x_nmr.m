function [xss,ts,nss,uss,pss,sigss,qs,dps,as] = sponge_x_nmr(dir_prefix,spongepar)
    % Scaling:
    %   n is porosity
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
    %   n0 is the initial porosity
    %   L is the initial length of the sponge
    %   T = mu*(L^2)/(k0*M);
    %
    % Geometry:
    %   This code takes the left and right edges of the sponge to be at x=a(t) and x=b=L, respectively,
    %   where the left edge is free and the right edge is fixed. As in the paper, the code assumes 
    %   without loss of generality that a0/L=0 and b0/L=b/L=1.
    %
    % INPUTS
    %   dir_prefix: Directory in which to save the results.
    %   spongepar: Structure with parameters and options.
    %     *.q: Dimensionless flow rate
    %     *.dp: Dimensionless pressure drop
    %       - Exactly one of the two above must be passed in. The other should be absent or NaN.
    %     *.sigstar: Dimensionless effective stress applied at x=a(t)
    %     *.n0: Initial (relaxed) porosity
    %     *.tobs: Times at which to save and return the solution.
    %
    % OPTIONAL INPUTS
    %   spongepar
    %     *.plot_switch: Plot occasionally while solving [1, default], or not [0].
    %     *.use_force: If there is an existing results file, ignore it [-1] or use it [any other value, default]
    %     *.stress_law: 
    %       - linear: Linear elasticity [default]
    %       - log: Hencky elasticity
    %     *.perm_law:
    %       - const: Constant permeability, where k is k/k0 = 1 [default]
    %       - KC: normalized Kozeny-Carman, where k is k/k0 = (((1-n0)^2)/(n0^3))*((n^3)/((1-n)^2))
    
    % Set some default values
    if ~isfield(spongepar,'use_force')
       spongepar.use_force = 0;
    end
    use_force = spongepar.use_force;
    if ~isfield(spongepar,'plot_switch')
        spongepar.plot_switch = 0;
    end
    plot_switch = spongepar.plot_switch;
    if ~isfield(spongepar,'stress_law')
        spongepar.stress_law = 'linear';
    end
    stress_law = spongepar.stress_law;
    if ~isfield(spongepar,'perm_law')
        spongepar.perm_law = 'const';
    end
    perm_law = spongepar.perm_law;
    
    % Unpack physical parameters
    n0 = spongepar.n0; % [-] initial porosity
    if ~isfield(spongepar,'q')
        spongepar.q = NaN;
    end
    q = spongepar.q;
    if ~isfield(spongepar,'dp')
        spongepar.dp = NaN;
    end
    dp = spongepar.dp;
    if isnan(q) & ~isnan(dp)
        bc_fluid = 'dp-fixed';
    elseif ~isnan(q) & isnan(dp)
        bc_fluid = 'q-fixed';
    else
        error('Must input exactly ONE of q and dp. The other will be determined by the solution.')
    end
    if ~isfield(spongepar,'sigstar')
        spongepar.sigstar = 0;
    end
    sigstar = spongepar.sigstar;
    
    use = 'n';
    save_filename = ['sponge_x_nmr' '_sig_' stress_law '_k_' perm_law '_q' num2str(q) '_dp' num2str(dp) '_sigstar' num2str(sigstar) '_n0' num2str(n0) '.mat'];
    if exist(strcat(dir_prefix,save_filename))==2 && use_force~=-1
        % Use saved results automatically, unless use_force==-1, in which case ignore them.
        disp('Using existing sponge_nmr results...')
        load(strcat(dir_prefix,save_filename),'spongepar','xss','ts','nss','uss','pss','sigss','qs','dps','as');
    else
        
        % Unpack control parameters
        tobs = spongepar.tobs;
        
        a0 = 0;
		b0 = 1;
		
        N = 400;
        
		% Initial condition: Start in a relaxed state
        a = a0;
        b = b0;
        dx = (b-a)/N;
        xs = [(a+dx/2):dx:(b-dx/2)];
		n0s = n0*ones(1,N);
        
        ns = n0s;
		Vs0 = sum(dx.*(1-n0s)); % For n = porosity
        
        % Define a function for the permeability law
        if strcmp(perm_law,'const')
            k = @(n) ones(size(n)); % Constant
        elseif strcmp(perm_law,'KC')
            k = @(n) (((1-n0)^2)/(n0^3)) * ((n.^3)./((1-n).^2)); % Normalized Kozeny-Carman
        else
            error('Unknown permeability law.')
        end
        
        % Define a function sig(n) for the stress law, and its inverse n_from_sig(sig).
        if strcmp(stress_law,'linear')
            sig = @(n) (n-n0)/(1-n0); % Linear elasticity
            n_from_sig = @(s) n0 + (1-n0)*s;
        elseif strcmp(stress_law,'log')
            sig = @(n) log((1-n0)./(1-n))./((1-n0)./(1-n)); % Hencky elasticity
            n_from_sig = @(s) 1 - (1-n0)*((-s)./lambertw(-s)); % y=-lambertw(-x)/x solves x=ln(y)/y
        else
            error('Unknown stress law.')
        end
        
        % Inner boundary condition: Imposed effective stress
        if sigstar==0
            nstar = n0; % Relaxed. Below should give same result, but this is exact.
        else
            nstar = n_from_sig(sigstar);
        end
        kstar = k(nstar);
        
        if ~isnan(dp)
            % If the pressure is imposed, n(x=b) is fixed in time
            n_at_b = n_from_sig(sigstar-dp);
        else
            % If the flux is imposed n(x=b) changes in time
            n_at_b = NaN;
        end
        
        % Empty structures for saving
        ts = [];
        as = [];
        qs = [];
        dps = [];
        as_anl = [];
        qs_anl = [];
        dps_anl = [];
        
        if plot_switch==1
            fig1 = figure();
        end
        count = 0;
        
        % Analytical steady-state
        [xs_ss,ns_ss,us_ss,ps_ss,sigs_ss,q_ss,dp_ss,delta_ss] = sponge_x_ss(dir_prefix,spongepar);
        
        % Options for ODE solver
        options = odeset('OutputFcn',@ode_check_plot_save, ...
                         'RelTol',1E-10, ...
                         'AbsTol',1E-10);
        
        % Initial condition
        Y0 = [(ns');a];
        
        [T,Y] = ode15s(@odefun,tobs,Y0,options);
        % [T,Y] = ode45(@odefun,tobs,Y0,options);
        
        ts = T;
        nss = Y(:,1:N);
        as = Y(:,N+1);
        
        [xss,uss,sigss,pss,kss,qs,dps] = ode_post(nss,as);
        
        % Save the results
        if exist(dir_prefix)~=7
            mkdir(dir_prefix)
        end
        save(strcat(dir_prefix,save_filename),'spongepar','xss','ts','nss','uss','pss','sigss','qs','dps','as');
        
        if plot_switch==1 && exist('fig1')==1
            close(fig1)
            drawnow()
        end
        
    end
    
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [ndots,adot] = calcFF(ns,a)
        
        % Calculate grid
        dx = (b-a)/N;
        xs = [(a+dx/2):dx:(b-dx/2)];
        
        % Calculate stress field from porosity field
        sigs = sig(ns);
        
        % Calculate permeability field from porosity field
        ks = k(ns);
        kstar = k(nstar);
        
        q = calc_q(q,sigs,ks);
        
        adot = q + (2*ks(1)*kstar/(ks(1)+kstar))*(sigs(1)-sigstar)/(dx/2);  
        % adot = q + kstar*(sigs(1)-sigstar)/(dx/2);  % dx/2 because sig(1) is
        % at dx/2 and sigstar is at a
        % dont really understand other expression so sticking to these 2nd
        % ones
        
        % Transmissibilities (transmissivities?)
        TTs = 2* (1-ns(1:end-1)).*ks(1:end-1).* (1-ns(2:end)).*ks(2:end) ...
            ./(  (1-ns(1:end-1)).*ks(1:end-1) + (1-ns(2:end)).*ks(2:end) );
        % TTs = ( (1-ns(1:end-1)).*ks(1:end-1) + (1-ns(2:end)).*ks(2:end) )/2;
        % (using this 2nd expression)
        
        % Fluxes
        Fs_adv = zeros(1,N+1);
        Fs_adv(1,2:N) = ( q - (((b-(xs(1:end-1)+dx/2))*adot)/(b-a)) ).*((ns(1:end-1)+ns(2:end))/2); % Advective, but not upwinded. More accurate??  %% used a central difference approx here?
        Fs_dff = zeros(1,N+1);
        Fs_dff(1,2:N) = -TTs.*( (sigs(2:end)-sigs(1:end-1))/dx ); % Diffusive: Centered
        Fs = Fs_adv + Fs_dff; % Total
        FL = q - adot;
        FR = q;
        Fs(1,1) = FL; Fs(1,N+1) = FR; % Boundary conditions
        
        ndots = (adot/(b-a))*ns - (Fs(1,2:end)-Fs(1,1:end-1))/dx;
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [q] = calc_q(q,sigs,ks)
        
        % !! Leaves q untouched if bc_fluid is q-fixed.
        
        if strcmp(bc_fluid,'dp-fixed')
            
            % Calculate q (needed to evolve n)
            sig_b = sigstar - dp; % p(b)-p(a)=sig(b)-sig(a), dp=p(a)-p(b), sig(a)=sigstar --> sig_b = sigstar-dp
            n_b = n_from_sig(sig_b);
            k_b = k(n_b);
            q = -(2*ks(end)*k_b/(ks(end)+k_b))*(sig_b-sigs(end))/(dx/2);
            % q = -k_b*(sig_b-sigs(end))/(dx/2);
            
        end
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [dp] = calc_dp(dp,us,sigs)
        
        % !! Leaves dp untouched if bc_fluid is dp-fixed.
        
        if strcmp(bc_fluid,'q-fixed')
            
            % Calculate dp (not needed -- just for reference)
            n_b = n0 + (1-n0)*((b-b0)-us(end))/(dx/2); % (n-n0)/(1-n0) = dudx
            sig_b = sig(n_b);
            dp = sigstar - sig_b; % p(b)-p(a)=sig(b)-sig(a), dp=p(a)-p(b), sig(a)=sigstar --> sig_b = sigstar-dp
            
        end
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [a_from_u,us] = u_from_n(xs,ns)
        % Calculate the displacement field from the porosity
        
        % (J-1)/J = (n-n0)/(1-n0) = du/dx, where J = (1-n0)/(1-n) is the Jacobian determinant
        % Integrate to get u at the i-1/2 interfaces
        us = -( sum(dx*((ns-n0)/(1-n0))) - [0,cumsum(dx*((ns(1:end-1)-n0)/(1-n0)))] );
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
        
        ns = Y(1:N)';
        a = Y(N+1);
        [ndots,adot] = calcFF(ns,a);
        Ydot = [(ndots');adot];
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function status = ode_check_plot_save(T,Y,flag)
        
        % This is called after every successful timestep, but plots only occasionally
        count = count + 1;
        
        if strcmp(flag,'done')
            status = 0; % status = 0 --> Continue running
        else
            
            status = 0; % status = 0 --> Continue running
            
            % Extract the current solution
            t = T(1);
            ns = Y(1:N,1)';
            a = Y(N+1,1);
            
            % Plot stuff
            if strcmp(flag,'init') || ( plot_switch==1 && mod(count,10)==0 )
                odePlot(t,ns,a,b);
            end
            
            % Stop when the solution is very close to steady state
            sigs = sig(ns);
            ks = k(ns);
            [a_from_u,us] = u_from_n(xs,ns);
            q = calc_q(q,sigs,ks);
            dp = calc_dp(dp,us,sigs);
            % if ( strcmp(bc_fluid,'q-fixed') && abs((dp-dp_ss)/dp_ss)<(5E-4) ) || ( strcmp(bc_fluid,'dp-fixed') && abs((q-q_ss)/q_ss)<(5E-4) )
            if abs((a-delta_ss)/delta_ss)<(1E-6)
                disp('Stopping because the numerical solution appears to have reached steady state.')
                status = 1; % status = 1 --> Stop
            end
            
        end
        
    end
    % -----------------------------------------------------
    % -----------------------------------------------------
    function odePlot(t,ns,a,b)
        
        % Calculate grid
        dx = (b-a)/N;
        xs = [(a+dx/2):dx:(b-dx/2)];
        
        % Calculate stress field from porosity field
        sigs = sig(ns);
        
        % Calculate permeability field from porosity field
        ks = k(ns);
        
        % Calculate displacement field from porosity field
        [a_from_u,us] = u_from_n(xs,ns);
        
        q = calc_q(q,sigs,ks);
        dp = calc_dp(dp,us,sigs);
        
        % Record a, q, and dp vs. t so that we can plot them.
        ts = [ts,t];
        as = [as,a];
        qs = [qs,q];
        dps = [dps,dp];
        
        spongepar.tobs = t;
        [xs_anl,t_anl,ns_anl,us_anl,ps_anl,sigs_anl,q_anl,dp_anl,a_anl] = sponge_x_anl_linear(dir_prefix,spongepar);
        as_anl = [as_anl,a_anl];
        qs_anl = [qs_anl,q_anl];
        dps_anl = [dps_anl,dp_anl];
        
        disp([' - Plotting the solution at time t = ' num2str(t) '...'])
        
        figure(fig1)
        
        Xs = xs - us;
        subplot(321)
        hold off
        plot(xs_ss,ns_ss,'k-')
        hold on
        plot(xs,ns,'r-')
        plot(a,nstar,'r*')
        xlabel('x/L')
        ylabel('n')
        set(gca,'XLim',[0,1])
        title(['t = ' num2str(t)])
        
        % plot(xs+us_anl,ns_anl,'g--')
        plot(xs_anl,ns_anl,'g--')
        
        subplot(323)
        hold off
        plot(xs_ss,us_ss,'k-')
        hold on
        plot(xs,us,'r-')
        plot(a_from_u,(a_from_u-a0),'r*')
        plot(a,(a-a0),'ko')
        xlabel('x/L')
        ylabel('u/L')
        set(gca,'XLim',[0,1])
        
        % plot(xs+us_anl,us_anl,'g--')
        plot(xs_anl,us_anl,'g--')
        
        subplot(325)
        hold off
        plot(xs_ss,sigs_ss,'k-')
        hold on
        plot(xs,sigs,'r-')
        plot(a,sigstar,'r*')
        xlabel('x/L')
        ylabel('\sigma/M')
        set(gca,'XLim',[0,1])
        
        % plot(xs+us_anl,sigs_anl,'g--')
        plot(xs_anl,sigs_anl,'g--')
        
        subplot(322)
        hold off
        plot([min(ts),max(ts)],[delta_ss,delta_ss],'k--')
        hold on
        plot(ts,as-a0,'k-')
        xlabel('t/T')
        ylabel('\delta/L')
        
        plot(ts,as_anl-a0,'g--')
        
        subplot(324)
        hold off
        plot([min(ts),max(ts)],[q_ss,q_ss],'k--')
        hold on
        plot(ts,qs,'k-')
        xlabel('t/T')
        ylabel('\mu{}q/k_0M')
        
        plot(ts,qs_anl,'g--')
        
        subplot(326)
        hold off
        plot([min(ts),max(ts)],[dp_ss,dp_ss],'k--')
        hold on
        plot(ts,dps,'k-')
        xlabel('t/T')
        ylabel('dp/M')
        
        plot(ts,dps_anl,'g--')
        
        % Check the result at the boundaries
        a_err = abs(a-a_from_u)/(b0-a0);
        if a_err>(1E-5)
            disp(['*** WARNING: rel err in displacement at left boundary = ' num2str(a_err)]);
        end
        
        % Check volume conservation
        Vs = sum(dx.*(1-ns));
        Vs_err = (Vs-Vs0)/Vs0;
        if abs(Vs_err)>1E-5
            disp(['*** WARNING: rel err in solid volume = ' num2str(Vs_err)])
        end
        
        drawnow
        
    end
    
    % -----------------------------------------------------
    % -----------------------------------------------------
    function [xss,uss,sigss,pss,kss,qs,dps] = ode_post(nss,as)
        
        % Empty structures for saving
        xss = zeros(size(nss));
        uss = zeros(size(nss));
        sigss = zeros(size(nss));
        pss = zeros(size(nss));
        kss = zeros(size(nss));
        qs = zeros(size(as));
        dps = zeros(size(as));
        
        % Loop over time
        for it=1:size(nss,1)
            
            ns = nss(it,:);
            a = as(it);
            
            % Calculate grid
            dx = (b-a)/N;
            xs = [(a+dx/2):dx:(b-dx/2)];
            
            % Calculate displacement field from porosity field
            [a_from_u,us] = u_from_n(xs,ns);
            
            % Calculate stress field from porosity field
            sigs = sig(ns);
            
            % Calculate permeability field from porosity field
            ks = k(ns);
            
            q = calc_q(q,sigs,ks);
            dp = calc_dp(dp,us,sigs);
            
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