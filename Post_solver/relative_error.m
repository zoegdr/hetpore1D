% Relative Error

%% Set parameters

params.p = 5;
params.Phi0 = 0.55;
params.nu = 0.3;
params.c = 1/16;

params.perm_law = 'KC';
params.stress_law = 'neo';
params.damage = '';
params.lp = 0;
params.dp = 0;
params.ls = 0;
params.ds = 0;

params.Astar = 0.1;
params.omega = 10;

%% Calculate errors
% or load error.mat

% Ns = [50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000];
% error = zeros(1,length(Ns));
% 
% % Initiate
% params.N = 25;
% [params,Ts,Zs,Phis,U_wall,Uss,Ps,Qs,Ss,dUdZ,Fluxes] = tendon_uniaxial_cyclic_load_stressdiff(params);
% Uss_new = Uss(end,1);
% 
% % Loop through Ns
% for i=1:length(Ns)
%     Uss_old = Uss_new; % Previous new value becomes old value
%     params.N = Ns(i); % Run for new N
%     % [params,Ts,Zs,Phis,U_wall,Uss,Ps,Qs,Ss,dUdZ,Fluxes] = tendon_uniaxial_cyclic_load_stressdiff(params);
%     [params,Ts,Zs,Phis,Uss,Ss,Ps,Qs,dUdZ,dPs,ks] = cylic_uniaxial_Lag_disp_2(params);
%     Uss_new = Uss(end,1); % Assign new value
%     error(i) = abs((Uss_new-Uss_old)/Uss_new); % Relative error
% end

%% Make plot

fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
set(fig,'Position',[0 0 8.2 8.2]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[8.2,8.2])

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters') 
set(gca, 'InnerPosition',[1.7 1.3 6 6])% This is the relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]

box('on');
hold("on")

plot(log(Ns),log(error_AL),'Color',"k",'LineWidth',1)
plot(log(Ns),log(error_AD),'Color',"k",'LineWidth',1,'LineStyle',':')

set(gca,'FontName','Times','FontSize',12);
xlabel('$\log(N)$','Interpreter','latex','Rotation',0);
ylabel('$\log(\varepsilon_{rel})$','Interpreter','latex');
ylim([-15 5])
xticks([4 5 6 7]);
% yticks([0 0.005 0.01 0.015])

legend('AS','AD','Location','best','NumColumns',1,'Box','off','IconColumnWidth',15);

print('relative-error','-dpdf','-r0')