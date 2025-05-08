% plot compare linear, intermediate and non-linear TEST

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

t1 = (params.p*pi)/params.omega;
t2 = (params.p*pi-pi/2)/params.omega;
t3 = (params.p*pi-pi)/params.omega;

[~,index1_k0] = min(abs(Ts-t1));
[~,index2_k0] = min(abs(Ts-t2));
[~,index3_k0] = min(abs(Ts-t3));

[~,index1_lin] = min(abs(t-t1));
[~,index2_lin] = min(abs(t-t2));
[~,index3_lin] = min(abs(t-t3));

[~,index1_KC] = min(abs(Ts_KC-t1));
[~,index2_KC] = min(abs(Ts_KC-t2));
[~,index3_KC] = min(abs(Ts_KC-t3));


%% Plot

ymin = 0;
ymax = 0.3;

hold on
% Shading (deloading)
fill([0 1 1 0], [ymin ymin ymax ymax], [0.97,0.97,0.97],EdgeColor=[0.97,0.97,0.97])

x2 = [Zs, fliplr(Zs)];

fill(x2, [min(dUdZ_KC,[],1), max(dUdZ_KC,[],1)], [0.5,0.7,0.2],'LineStyle','none');
alpha(0.05)
plot(Zs,dUdZ_KC(index2_KC,:),'color',[0.5,0.7,0.2],'LineWidth',1)

fill(x2, [min(dUdZ,[],1), max(dUdZ,[],1)], [0.6,0.1,0.2],'LineStyle','none');
alpha(0.1)
plot(Zs,dUdZ(index2_k0,:),'color',[0.6,0.1,0.2],'LineWidth',1)

fill(x2, [min(dwdz,[],1), max(dwdz,[],1)], [0.2,0.5,1],'LineStyle','none');
alpha(0.2)
plot(Zs,dwdz(index2_lin,:),'color',[0.2,0.5,1.0],'LineWidth',1)

% fill(x2, [dUdZ_KC(index1_KC,:), dUdZ_KC(index3_KC,:)], [0.5,0.7,0.2],'LineStyle','none');
% alpha(0.05)
% plot(Zs,dUdZ_KC(index2_KC,:),'color',[0.5,0.7,0.2],'LineWidth',1)
% 
% fill(x2, [dUdZ(index1_k0,:), dUdZ(index3_k0,:)], [0.6,0.1,0.2],'LineStyle','none');
% alpha(0.1)
% plot(Zs,dUdZ(index2_k0,:),'color',[0.6,0.1,0.2],'LineWidth',1)
% 
% fill(x2, [dwdz(index1_lin,:), dwdz(index3_lin,:)], [0.2,0.5,1],'LineStyle','none');
% alpha(0.2)
% plot(Zs,dwdz(index2_lin,:),'color',[0.2,0.5,1.0],'LineWidth',1)

%% Axes/Labels

set(gca, "Layer", "top")
ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
xlabel('$Z$','Interpreter','latex','Rotation',0);

ylim([ymin ymax]) % uniform
yticks([ymin ymin+(ymax-ymin)/2 ymax])
xlim([0 1])
xticks([0 0.5 1])

% legend('','Non-linear, KC', 'Non-linear, $k_0$', 'Linear','Interpreter','Latex',Location='best')

%% Export ----------------------------------------------------------------------

% print('Compare-models','-depsc','-loose')
