% Plot strain and flux metrics for varying frequency

% Load metric_var_d.mat

%% SET UP

dss = 0:0.02:0.9;
omegas = 0.5:0.5:50;
% c = slanCM('tab20');
% cmap = [c(105,:);c(120,:)];
c = slanCM('tab20c');
cmap = [c(52,:);c(78,:);c(92,:)];

%% PLOT ----------------------------------------------------------

% Create figure
fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 13 6])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[13 6])


% Create subfigures

% Q_int AD against d
ax11 = axes('Units','centimeters','InnerPosition',[1.5 1 2.8 4]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',10);
    set(ax11,'ColorOrder',cmap)

    annotation("textbox", 'String','a)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [0 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_d{2,2}; % Stiff, Q_int AD
    
    % AD ---
    plot(dss(:),(metric(:,1)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric(:,2)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric(:,3)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([-0.1 1.5])
    ylabel(ax11,{'$\Delta Q$'},'Interpreter','latex','Rotation',0)
    xlabel('$d$','Interpreter','latex')

% Q_int AD against omega
ax12 = axes('Units','centimeters','InnerPosition',[6.5 1 5 4]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    annotation("textbox", 'String','b)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [5 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_freq{2,2};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.1 10.5])
    % yticks([0 1 2 3 4])
    % ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')
    
    legend('$l = 0.25$','$l = 0.5$','$l=0.75$','Interpreter','latex','Units','centimeters','Position',[6.5 3.2 1 5],'NumColumns',4,'Box','off','IconColumnWidth',15,'FontSize',12);




%% Export
print('Effect-of-params-AD-stiff','-dpdf','-r0')    