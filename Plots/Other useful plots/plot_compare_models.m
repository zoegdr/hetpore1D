% plot compare linear, intermediate and non-linear

%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
% set(gcf,'Position',[0 0 9.5 8])
set(gcf,'Position',[0 0 12.2 8.75]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters') 
% set(gca, 'InnerPosition',[2 1.5 5.5 5.5])
set(gca, 'InnerPosition',[2 1.5 9.5 5.5])% This is the relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]

box('on');
set(gca,'FontName','Times','FontSize',12);

%% Choose time to plot

tc = (params.p*pi-pi/2)/params.omega;
[~,index_k0] = min(abs(Ts-tc));
[~,index_lin] = min(abs(t-tc));
[~,index_KC] = min(abs(Ts_KC-tc));

%% Plot

ymin = 0.1;
ymax = 0.23;

hold on
% Shading (deloading)
fill([0 1 1 0], [ymin ymin ymax ymax], [0.97,0.97,0.97],EdgeColor=[0.97,0.97,0.97])
plot(Zs,dUdZ_KC(index_KC,:),'color',[0.5,0.7,0.2],'LineWidth',1)
plot(Zs,dUdZ(index_k0,:),'color',[0.6,0.1,0.2],'LineWidth',1)
plot(Zs,dwdz(index_lin,:),'color',[0.2,0.5,1.0],'LineWidth',1)

%% Axes/Labels

set(gca, "Layer", "top")
ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
xlabel('$Z$','Interpreter','latex','Rotation',0);

ylim([ymin ymax]) % uniform
yticks([0.1 0.15 0.2])
xlim([0 1])
xticks([0 0.5 1])

legend('','Non-linear, KC', 'Non-linear, $k_0$', 'Linear','Interpreter','Latex',Location='best')

%% Export ----------------------------------------------------------------------

% print('Compare-models','-depsc','-loose')

