% Plot strain and flux metrics for varying frequency

function plot_perm_metrics(metric_var_d,metric_var_freq)

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
set(fig,'Position',[0 0 17 6])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[17 6])


% Create subfigures

% dU_int AL against d
ax11 = axes('Units','centimeters','InnerPosition',[1.5 1 2.8 4]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',10);
    set(ax11,'ColorOrder',cmap)

    annotation("textbox", 'String','a)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [0 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_d{3,1}; % Perm, dU_int AL
    
    % AL ---
    plot(dss(:),(metric(:,1)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric(:,2)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric(:,3)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([-0.15 0.15])
    yticks([-0.15 -0.1 -0.05 0 0.05 0.1 0.15])
    ax11.YAxis.Exponent = -1;
    xlabel('$d$','Interpreter','latex')
    ylabel('$\Delta U$','Interpreter','latex','Rotation',0)

% Q_int AL against d
ax12 = axes('Units','centimeters','InnerPosition',[6.5 1 2.8 4]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    annotation("textbox", 'String','b)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [5 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_d{3,2}; % Perm, Q_int AL

    plot(dss(:),(metric(:,1)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric(:,2)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric(:,3)-metric(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    xlim([0 0.9])
    ylim([-1.5 0.1])
    xticks([0 0.5 0.9])
    xlabel('$d$','Interpreter','latex') 
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)

% Q_int AL against d
ax13 = axes('Units','centimeters','InnerPosition',[11.5 1 5 4]);
    hold(ax13,"on")
    box(ax13,"on")
    set(ax13,'FontName','Times','FontSize',10);
    set(ax13,'ColorOrder',cmap)

    annotation("textbox", 'String','c)', 'FontName','Times','FontSize',12, 'Units','centimeters','Position', [10 4.5 0.5 0.5],'LineStyle','none')

    metric = metric_var_freq{3,2}; % Perm, Q_int AL

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    ylim([-0.3 0])
    xticks([0.5 10 20 30 40 50])
    % yticks([-0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0])
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')
    
    legend('$l = 0.25$','$l = 0.5$','$l=0.75$','Interpreter','latex','Units','centimeters','Position',[8.5 3.2 1 5],'NumColumns',4,'Box','off','IconColumnWidth',12,'FontSize',10);

end 