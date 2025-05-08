% Plot damaged diffusivity

% plot diffusivity x 3

params = Perm_D_params{1};

omega = params.omega;

D0 = (Perm_D_Phis{1}/params.Phi0+1).^3./((1+Perm_D_Phis{1}).^2);
D_local = (Perm_D_Phis{3}/params.Phi0+1).^3./((1+Perm_D_Phis{3}).^2) .* (1-Perm_D_params{3}.dp*exp(-(Zs-Perm_D_params{3}.lp).^2/(2*(Perm_D_params{3}.c)^2)));
D_pos = (Perm_D_Phis{5}/params.Phi0+1).^3./((1+Perm_D_Phis{5}).^2) .* (-Perm_D_params{5}.dp*Zs+1);
D_neg = (Perm_D_Phis{11}/params.Phi0+1).^3./((1+Perm_D_Phis{11}).^2) .* (Perm_D_params{11}.dp*Zs+1-Perm_D_params{11}.dp);

%% Limits

tmin = 18*pi/omega;
tmid = 19*pi/omega;
tmax = 20*pi/omega;

ymin = 0.7;
ymax = 2.5;


%% Create figure ----------------------------------------------------------

figure

xlim([tmin,tmax]);
xticks([tmin tmid tmax])
xticklabels({'','', ''})

% ylim([1,1.7]);
% yticks([1 2]);

% Find indexes
% [~,index1] = min(abs(Perm_D_Ts{1}-tmin));
% [~,index3] = min(abs(Perm_D_Ts{1}-tmax));

hold on
% Shading (deloading)
fill([tmid tmax tmax tmid], [1 1 2 2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
% Plot
% plot(Perm_D_Ts{1}(index1:index3),trapz(Zs,D0(index1:index3,:),2),'LineWidth',1,'Color','k');
% plot(Perm_D_Ts{3}(index1:index3),trapz(Zs,D_local(index1:index3,:),2),'LineWidth',1,'Color','r');
% plot(Perm_D_Ts{5}(index1:index3),trapz(Zs,D_pos(index1:index3,:),2),'LineWidth',1,'Color','g');
% plot(Perm_D_Ts{11}(index1:index3),trapz(Zs,D_neg(index1:index3,:),2),'LineWidth',1,'Color','b');

% xlabel('time')
% ylabel('D')

% Unpack
x = Perm_D_Ts{1}; z = D0;

% Find time indexes
[~,index1] = min(abs(x-tmin));
[~,index3] = min(abs(x-tmax));

% Plot 10 equally spaced lines over one cycle between 0<Z<1
plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));

xlabel('Z')
ylabel('D')
