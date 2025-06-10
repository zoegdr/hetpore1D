function [params,Ts,Zs,Phis,Uss,Ss,Ps,Qs,dUdZ,dPs,ks] = run_simulation(loading,params)

if strcmp(loading,'stress')
    [params,Ts,Zs,Phis,~,Uss,Ps,Qs,Ss,dUdZ,~,dPs,ks] = tendon_uniaxial_cyclic_load_stressdiff(params);
elseif strcmp(loading,'disp')
    [params,Ts,Zs,Phis,Uss,Ss,Ps,Qs,dUdZ,dPs,ks] = cylic_uniaxial_Lag_disp(params);
end

end