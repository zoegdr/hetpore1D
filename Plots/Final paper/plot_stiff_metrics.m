% Plot strain and flux metrics for stiffness damage

% Load metric_var_d.mat and metric_var_freq.mat

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
set(fig,'Position',[0 0 18 6])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[18 6])


% Create subfigures

% Q_int AD against d
ax11 = axes('Units','centimeters','InnerPosition',[1.5 1 2.5 3.75]);
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
ax12 = axes('Units','centimeters','InnerPosition',[5 1 3.75 2.5]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    annotation("textbox", 'String','b)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [4.2 4.5 0.5 0.5],'LineStyle','none')

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

% dU_int AS against d
ax13 = axes('Units','centimeters','InnerPosition',[10.5 1 2.5 3.75]);
    hold(ax13,"on")
    box(ax13,"on")
    set(ax13,'FontName','Times','FontSize',10);
    set(ax13,'ColorOrder',cmap)

    annotation("textbox", 'String','c)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [9 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_d{1,1}; % Stiff, dU_int AS
    
    % AD ---
    plot(dss(:),(metric(:,1)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric(:,2)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric(:,3)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([0 1.5])
    % yticks([0 0.025 0.05])
    ylabel(ax13,{'$\Delta U$'},'Interpreter','latex','Rotation',0)
    xlabel('$d$','Interpreter','latex')

% dU_int AS against omega
ax14 = axes('Units','centimeters','InnerPosition',[14 1 3.75 2.5]);
    hold(ax14,"on")
    box(ax14,"on")
    set(ax14,'FontName','Times','FontSize',10);
    set(ax14,'ColorOrder',cmap)

    annotation("textbox", 'String','d)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [13.2 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_freq{1,1}; % Stiff, dU_int AS

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([0 10])
    % yticks([0 1 2 3 4])
    % ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    % ax12.YAxis.Exponent = -2;
    xlabel('$\omega$','Interpreter','latex')
    
    legend('$l = 0.25$','$l = 0.5$','$l=0.75$','Interpreter','latex','Units','centimeters','Position',[9 3.2 1 5],'NumColumns',4,'Box','off','IconColumnWidth',12,'FontSize',10);




%% Export
print('Effect-of-params-AD-stiff','-dpdf','-r0')    