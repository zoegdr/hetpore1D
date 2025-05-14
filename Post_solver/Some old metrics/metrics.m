function [HSI,HSI_max,HSI_Z,penetration,max_pen,max_pen_time] = metrics(x,y,z,params)

% metrics to quantify data from poroelastic framework

% x is time
% y is space
% z is which quantity the metric relates to

%% Unpack physical parameters ---------------------------------------------

A = params.Astar;
p = params.p;
omega = params.omega;
T = 2*pi/omega;

%% Calculate metrics ------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interval of high strain (>A-0.03) at each Z (working title)

% matrix of elements = 0 if below = 1 if above
M1 = z > A-A/3;%A-A/10; %A-0.03; Threshold is chosen depending on scenario which will affect the start and finish array sizes at different Zs
M1 = flip(M1,2); % to avoid transient region

% Run over each Z
for i=1:length(y)
    start = strfind((M1(:,i))',[0 1])+1; % find when switches from 0 to 1. +1 because finds the position of the 0 not the 1
    finish = strfind((M1(:,i))',[1 0]); % find when switches from 1 to 0
    if isempty(start) % if no regions above threshold
        start=1;
        finish=1;
    else
        % select value for period that is far enough past transient, but
        % not right at the end to 1) avoid delayed response not captured
        % and 2) initial response offsets some of the values for certain Zs
        % -> for some we are then calculating the metric in two different
        % periods (eg. local damage away from applied load) but as each
        % period is identical this doesn't matter
        % possibly needs to be modified according to scenario
        start = start(2);
        finish = finish(2); 
    end
    HSI(i) = x(finish)-x(start);
    HSI(i) = HSI(i)/T; % calculate corresponding fraction of period
end

[HSI_max,HSI_Z] = max(HSI); % Maximum High Strain Interval and spatial point (Z) at which it occurs
HSI_Z = y(end)-y(HSI_Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Penetration length
% note if penetrates all the way this metric does not work as does not see
% pattern of 0 to 1. Can just see penetrates all the way anyway.

% matrix of elements = 0 if below = 1 if above
M2 = z > A-A/4; %A-0.03; Threshold is chosen depending on scenario which will affect the start and finish array sizes at different Zs

% Run over each timestep
for i=1:length(x)
    p = strfind(M2(i,:),[0 1])+1; % find when switches from 0 to 1. +1 because finds the position of the 0 not the 1
    if isempty(p) % if no regions above threshold i.e. no penetration
        p = length(y);
    else
        p = p(1); % only care about the first one as this is the furthest is goes (Z=0 is away from applied load) in case in a scenario where strain is heterogeneous
    end
    penetration(i) = y(p);
    pen_percent(i) = (y(end)-y(p))*100; % calculate as a percentage
end

[max_pen,max_pen_time] = min(penetration); % max penetration length and time index at which it occurs (which is a minimum furthest is Z=0)
[max_pen_percent,~] = max(pen_percent);
max_pen_time = x(max_pen_time)/T; % at what point in the period
end