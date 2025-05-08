% Plot uniform metrics

%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 6 8]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters') 
set(gca, 'InnerPosition',[1.5 1.5 4.2 5])% This is the relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]

box('on');
set(gca,'FontName','Times','FontSize',12);

%% Choose metric to plot
% porosity or flux

x = 1:1:30;

% var = 'flux';
var = 'flux trans';
% var = 'porosity';

hold on

if strcmp(var,'flux')

    plot(x,max_flow,'Color',[0.0,0.4,0.7],'LineWidth',1);
    plot(x,min_flow,'Color',[0.6,0.9,1.0],'LineWidth',1);
    plot(x,delta_Q,'Color',"k","LineStyle","--","LineWidth",1);

    ylabel('$Q$','Interpreter','latex','Rotation',0);
    lgd = legend('$Q^{max}$','$Q^{min}$','$\Delta Q$','Interpreter','latex','NumColumns',3);
    ylim([-1 1])
    yticks([-1 -0.5 0 0.5 1])

elseif strcmp(var,'flux trans')

    plot(x,max_flow_trans,'Color',[0.0,0.4,0.7],'LineWidth',1);
    plot(x,min_flow_trans,'Color',[0.6,0.9,1.0],'LineWidth',1);
    plot(x,delta_Q_trans,'Color',"k","LineStyle","--","LineWidth",1);

    ylabel('$Q$','Interpreter','latex','Rotation',0);
    lgd = legend('$Q^{max}_{trans}$','$Q^{min}_{trans}$','$\Delta Q$','Interpreter','latex','NumColumns',3);
    ylim([-1 1])
    yticks([-1 -0.5 0 0.5 1])

elseif strcmp(var,'porosity')

    plot(x,max_total_porosity,'Color',[0.0,0.4,0.7],'LineWidth',1);
    plot(x,min_total_porosity,'Color',[0.6,0.9,1.0],'LineWidth',1);
    plot(x,delta_phi,'Color',"k","LineStyle","--","LineWidth",1);
    
    ylabel('$\phi$','Interpreter','latex','Rotation',0);
    lgd = legend('$\phi^{max}_{tot}$','$\phi^{min}_{tot}$','$\Delta \phi+\phi_0$','Location','northoutside','Interpreter','latex','NumColumns',3);
    ylim([0.55 0.63])
    yticks([0.55 0.59 0.63])

end

xlabel('$\omega_1$','Interpreter','latex','Rotation',0);
xlim([1 30]);
xticks([1 10 20 30])
lgd.Units = 'centimeters';
lgd.Position = [1.5 6.5 4.2 1.3];
lgd.IconColumnWidth = 7;
legend('boxoff')

%% Export ----------------------------------------------------------------------

% print('Max-Min-Q-omegas','-depsc','-loose')
print('Max-Min-Q-trans-omegas','-depsc','-loose')
% print('Max-Min-Porosity-omegas','-depsc','-loose')