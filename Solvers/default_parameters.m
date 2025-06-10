function params_d = default_parameters(loading)

    params_d.N = 400;
    params_d.p = 20;
    params_d.Phi0 = 0.55;
    params_d.nu = 0.3;
    params_d.omega = 10;
    params_d.perm_law = 'KC';
    params_d.stress_law = 'neo';
    params_d.damage = 'local-';
    params_d.lp = 0;
    params_d.ls = 0;
    params_d.dp = 0;
    params_d.ds = 0;
    params_d.v = 1/10;

    if strcmp(loading,'stress')
        params_d.Astar = 0.2;
    elseif strcmp(loading,'disp')
        params_d.Astar = 0.1;
    else
        error('Unknown loading type.')
    end

end