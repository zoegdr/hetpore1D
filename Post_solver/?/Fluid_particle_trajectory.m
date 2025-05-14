% Fluid particle trajectory

z_eul = Uss + Zs;
vs = -Qs; % vs = -Qs
vf = -vs.*((1-params.Phi0)./(Phis+params.Phi0)); % vf = vs*(1-J/Phi) (unnormalised Phi)

uf = cumtrapz(Ts,vf,1);
