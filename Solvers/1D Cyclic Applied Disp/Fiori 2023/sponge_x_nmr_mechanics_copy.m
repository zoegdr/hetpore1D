function [xss,ts,nss,uss,sigss,vfss,vsss,as] = sponge_x_nmr_mechanics_copy(params)

    % if ~strcmp(dir_prefix(end),'/')
    %     dir_prefix = strcat(dir_prefix,'/');
    % end
    % params.dir_prefix = dir_prefix;
    
    % Set some default values
    if ~isfield(params,'post_switch')
        params.post_switch = 1;
    end
    post_switch = params.post_switch;
    if ~isfield(params,'stress_law')
        params.stress_law = 'linear';
    end
    stress_law = params.stress_law;
    if ~isfield(params,'perm_law')
        params.perm_law = 'const';
    end
    perm_law = params.perm_law;
    if ~isfield(params,'diff')
        params.diff = 'compact';
    end
    if ~isfield(params,'tol')
        params.tol = 1E-8;
    end

    % Unpack physical parameters
    n0 = params.n0; % [-] initial porosity
    A = params.Astar;
    T = 2*pi/params.omega;
    params.a = @(t) (A/2)*(1-cos(2*pi*t/T));
    params.adot = @(t) (pi*A/T)*sin(2*pi*t/T);

    % save_filename = [ 'results_nmr_A_' num2str(A) '_T_' num2str(T/pi) 'pi' '_n0_' num2str(n0) '_'  '.mat'];
    
    % if exist(strcat(dir_prefix,save_filename))==2 
    %     % Use saved results automatically, unless use_force==-1, in which case ignore them.
    %     disp(['Using existing sponge_nmr results in ' save_filename '...'])
    %     load(strcat(dir_prefix,save_filename),'params','ts','nss','as');
    %     if abs(ts(end)-params.tobs(end))>1E-8
    %         warning('This simulation did not finish successfully! Re-saving with NaNs...')
    %         ts = NaN;
    %         nss = NaN;
    %         as = NaN;
    %         save(strcat(dir_prefix,save_filename),'params','ts','nss','as');
    %     end
    % else
        
        % Unpack control parameters
        params.tobs = [0 params.p*T];
        tobs = params.tobs;
        N = params.N;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Chebyshev differentiation matrices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if exist('chebdif')~=2
            error('The solver uses Chebyshev spectral differentiation matrices provided by the function chebdif.m by Weideman & Reddy, which itself is readily available on the internet.')
        end
        [x,D] = chebdif(N,2);
        
        % The chebdif function above produces a grid of x = [1,-1].
        % Rescale to x = [0,1].
        Xs = (1-x)/2;
        DDX = -2*D(:,:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        params.Xs = Xs;
        params.DDX = DDX;
        IIX = zeros(size(DDX)); IIX(2:end,2:end) = inv(DDX(2:end,2:end)); % Such that IIX*u integrates u from X(1) to X
        params.IIX = IIX;       
        
        % Initial condition: Start in a relaxed state
        a0 = params.a(0);
        x0s = (1-a0)*Xs + a0;
        
        n0s = n0*ones(1,N);
        params.a0 = a0;       
        
        % Define a function for the permeability law
        if strcmp(perm_law,'const')
            k = @(n) ones(size(n)); % Constant
        elseif strcmp(perm_law,'KC')
            k = @(n) (((1-n0)^2)/(n0^3)) * ((n.^3)./((1-n).^2)); % Normalized Kozeny-Carman
        else
            error('Unknown permeability law.')
        end
        params.k = k;
        
        global count;
        count = 0;

        % Mass matrix for enforcing BCs as algebraic constraints (see odeFun)
        MM = speye(N,N);
        MM(1,:) = 0;
        MM(N,:) = 0;

        JPattern = sparse(N,N);
        JPattern(1:N,1:N) = (DDX*DDX)~=0;

        % Options for ODE solver
        options = odeset('OutputFcn',@(T,Y,flag) ode_check_save(T,Y,flag,params) , ...
                         'Mass',MM , 'MassSingular', 'yes', ...
                         'JPattern', JPattern, ...
                         'RelTol',params.tol, ...
                         'AbsTol',params.tol);

        % Initial condition
        Y0 = [ n0s(:) ];
        
        [T,Y] = ode15s(@(T,Y) odefun(T,Y,params),tobs,Y0,options);
        % if abs(T(end)-params.tobs(end))>1E-8
        %     warning('Simulation did not finish successfully.')
        %     ts = NaN;
        %     nss = NaN;
        %     as = NaN;
        % else
            %Unpack solution
            ts = T;
            nss = Y(:,1:N);
            as = params.a(ts);
        % end

        % % Save the results
        % if exist(dir_prefix)~=7 && ~strcmp(dir_prefix,'./')
        %     mkdir(dir_prefix)
        % end
        % save(strcat(dir_prefix,save_filename),'params','ts','nss','as');
        
    
    % xss = NaN;
    % uss = NaN;
    % sigss = NaN;
    % vfss = NaN;
    % vsss = NaN;
    % 
    % if isnan(ts(1))
    %     warning('Simulation did not finish successfully.')
    % else

        % Minimal post-processing
        Xs = params.Xs;
        xss = zeros(size(nss));
        for it=1:length(ts)
            a = as(it);
            xs = (1-a)*Xs + a;
            xss(it,:) = xs;
        end

        % Full post-processing
        if post_switch==1
            [uss,sigss,vfss,vsss] = ode_post(ts,nss,as,params);
        end

    end

% end
% -----------------------------------------------------
% -----------------------------------------------------
function Ydot = odefun(T,Y,params)
    
    %Unpack Y
    N = params.N;
    t = T;
    ns = Y(1:N); ns = ns(:);
    a = params.a(t);

    %Unpack parameters
    Xs = params.Xs;
    DDx = (1/(1-a))*params.DDX;
    
    [sigs,ks,vfs,vss] = calc_sigs(ns,a,params);

    %Displacement-driven:
    adot = params.adot(t);
    
    FF = ns.*vfs;
    
    ndots = (1-Xs)*adot.*(DDx*ns) - DDx*FF;

    % Apply BCs
    ndots(1) = vss(1) - adot; % vss(1) = adot
    ndots(end) = FF(end); % vf(end) = 0

    Ydot = [ ndots ];
    
end
% -----------------------------------------------------
% -----------------------------------------------------
function [sigs,ks,vfs,vss] = calc_sigs(ns,a,params)

    %Unpack parameters
    DDx = (1/(1-a))*params.DDX;
    k = params.k;
    
    sigs = sigfun(ns,NaN,params);
    ks = k(ns);
    
    vfs = -((1-ns)./ns).*ks.*(DDx*sigs);
    vss = ks.*(DDx*sigs);

end
% -----------------------------------------------------
% -----------------------------------------------------
function [s,n] = sigfun(n,s,params)

    stress_law = params.stress_law;
    n0 = params.n0;
    
    % Define a function sig(n) for the stress law, and its inverse n_from_sig(sig).
    if strcmp(stress_law,'linear')
        sig = @(n) (n-n0)/(1-n0); % Linear elasticity
        n_from_sig = @(s) n0 + (1-n0)*s;
    elseif strcmp(stress_law,'lin-stiff2')
        % Linear elasticity with a steep divergence as phi --> phi_min
        % See Hewitt PRE 2016 eqs 2 and 15a.
        % nmin = 0;
        nmin = 0.01;
        sig = @(n) ((n-n0)./(1-n0))./(1-((n-n0)./(1-n0))/((nmin-n0)/(1-n0)));
        n_from_sig = @(s) n0 + (1-n0)*s./(1+s/((nmin-n0)/(1-n0))) ;
    elseif strcmp(stress_law,'log')
        sig = @(n) log((1-n0)./(1-n))./((1-n0)./(1-n)); % Hencky elasticity
        n_from_sig = @(s) 1 - (1-n0)*((-s)./lambertw(-s)); % y=-lambertw(-x)/x solves x=ln(y)/y
        % n_from_sig = @(s) fzero( @(n) sig(n)-s , [1E-10,1-1E-10] );
        if any( s > 0 )
            warning('n --> s for Hencky is multivalued in tension. Proceed with extreme caution.')
        end
    else
        error(['Unknown stress law: ' stress_law])
    end
    
    if isnan(s)
        s = sig(n);
    elseif isnan(n)
        n = n_from_sig(s);
    end

end
% -----------------------------------------------------
% -----------------------------------------------------
function [us] = u_from_n(ns,a,params)
    % Calculate the displacement field from the porosity
    
    a0 = params.a0;
    n0 = params.n0;
    IIx = (1-a)*params.IIX;
    us = (a-a0) + IIx*((ns(:)-n0)/(1-n0)); % du/dx = (n-n0)/(1-n0)
    
end
% -----------------------------------------------------
% -----------------------------------------------------
function status = ode_check_save(T,Y,flag,params)
    
    global count;
    
    if strcmp(flag,'done')
        disp([' - Done.'])
    else
        t = T(1);
        count = count + 1;
        if mod(count,10)==0
            disp([' - Solving at time t = ' num2str(t) ' (' num2str(100*t/max(params.tobs)) '%) ...'])
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
% -----------------------------------------------------
% -----------------------------------------------------
function [uss,sigss,vfss,vsss] = ode_post(ts,nss,as,params)    
    
    % Empty structures for saving
    uss = zeros(size(nss));
    vfss = zeros(size(nss));
    vsss = zeros(size(nss));
    sigss = zeros(size(nss));

    adot = params.adot(ts);
    
    N = params.N;
    IIX = params.IIX;
    
    % Loop over time
    for it=1:size(nss,1)
        
        t = ts(it);
        ns = nss(it,:);
        a = as(it);
        IIx = (1-a)*IIX;

        % Calculate displacement field from porosity field
        [us] = u_from_n(ns',a,params);
        
        [sigs,ks,vfs,vss] = calc_sigs(ns',a,params);
        
        % Save
        uss(it,:) = us';
        sigss(it,:) = sigs';
        vfss(it,:) = vfs';
        vsss(it,:) = vss';

    end
    
end
% -----------------------------------------------------
% -----------------------------------------------------