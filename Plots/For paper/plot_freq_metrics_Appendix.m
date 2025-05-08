% Plot strain and flux metrics for varying frequency

% Load metric_var_freq.mat

%% SET UP

omegas = 0.5:0.5:50;
c = slanCM('tab20c');
cmap = [c(52,:);c(78,:);c(92,:)];

%% PLOT ----------------------------------------------------------

% Create figure
fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 6.5 9])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[6.5,9])


% Create subfigures

% Stiff AL, Q_int
ax11 = axes('Units','centimeters','InnerPosition',[1.7 5.5 4 2.5]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',10);
    set(ax11,'ColorOrder',cmap)

    metric = metric_var_freq{1,2};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([0 1])
    % yticks([0 1 2 3 4])
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')
    
    legend('$l = 0.25$','$l = 0.5$','$l=0.75$','Interpreter','latex','Units','centimeters','Position',[3 6 1 5],'NumColumns',4,'Box','off','IconColumnWidth',15);

% Stiff AL, dU_int
ax12 = axes('Units','centimeters','InnerPosition',[1.7 1.5 4 2.5]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    metric = metric_var_freq{2,1};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([0 0.1])
    ylabel('$\Delta U$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')

%% Export
print('Effect-of-freq&l-Appendix','-dpdf','-r0')
  