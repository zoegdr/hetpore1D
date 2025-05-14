function [Q_in,Q_out,T_in,T_out] = calc_inflow_outflow(params,Ts,Qs)

% Calculate total inflow and outflow

% Reduce Q to one period
omega = params.omega;
tmin = 6*pi/omega;
tmax = 8*pi/omega;

[~,index1] = min(abs(Ts-tmin));
[~,index2] = min(abs(Ts-tmax));

Q = Qs(index1:index2,:);
T = Ts(index1:index2);
dZ = 1/params.N;

% Total inflow and outflow
Q_pos = Q > 0; % 1 if positive, 0 if negative
Q_zero = abs(Q) < 0.5; % 1 if smaller than 0.5, 0 otherwise
Q_in = 0;
Q_out = 0;


for n = 1:params.N

    t1 = strfind((Q_pos(:,n))',[0 1])+1; % find index when switches from 0 to 1. +1 because finds the position of the 0 not the 1
    t2 = strfind((Q_pos(:,n))',[1 0]); % find index when switches from 1 to 0

    Q_in = Q_in + dZ* trapz(T(t1:t2),Q(t1:t2,n));
    Q_out = Q_out + dZ* ( trapz(T(1:t1),Q(1:t1,n)) + trapz(T(t2:end),Q(t2:end,n)) );

    % Calculate time of outflow and inflow at the boundary (not at 1
    % because appear to get a longer flow at 2...)
    if n == 2
        % Find times when flow is sufficiently small to approximate as zero
        t3 = strfind((Q_zero(:,2))',[0 1])+1; % find index when switches from 0 to 1. +1 because finds the position of the 0 not the 1
        t4 = strfind((Q_zero(:,2))',[1 0]); % find index when switches from 1 to 0
    
        T_in = T(t2)-T(t1) -( T(t4(1))-T(t3(1)) + T(t4(2))-T(t3(2)) );
        T_out = T(t1)-T(1) + T(end)-T(t2) -( T(t4(1))-T(t3(1)) + T(t4(2))-T(t3(2)) );
    end

end



end