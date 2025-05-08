% Plot strain and flux metrics for varying frequency

% Load metric_var_d.mat

%% SET UP

dss = 0:0.02:0.9;
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
set(fig,'Position',[0 0 18 7.5])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[18 7.5])


% Create subfigures

% Stiff, dU_int
ax11 = axes('Units','centimeters','InnerPosition',[1.5 1 2.8 4]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',10);
    set(ax11,'ColorOrder',cmap)

    metric1 = metric_var_d{1,1}; % Stiff, dU_int AL
    metric2 = metric_var_d{2,1}; % Stiff, dU_int AD
    
    % AL ---
    plot(dss(:),(metric1(:,1)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric1(:,2)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric1(:,3)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % AD ---
    plot(dss(:),(metric2(:,1)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.25
    plot(dss(:),(metric2(:,2)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.5
    plot(dss(:),(metric2(:,3)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([-1.5 1.5])
    ylabel(ax11,{'$\Delta U$'},'Interpreter','latex','Rotation',0)
    xlabel('$d$','Interpreter','latex')

    legend('$l = 0.25$','$l = 0.5$', '$l = 0.75$','Interpreter','latex','Units','centimeters','Position',[8.9 3.8 1 5],'NumColumns',4,'Box','off','IconColumnWidth',15);


% Perm, dU_int
ax12 = axes('Units','centimeters','InnerPosition',[5.5 1 2.8 4]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    metric1 = metric_var_d{3,1}; % Perm, dU_int AL
    metric2 = metric_var_d{4,1}; % Perm, dU_int AD
    
    % AL ---
    plot(dss(:),(metric1(:,1)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric1(:,2)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric1(:,3)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % AD ---
    plot(dss(:),(metric2(:,1)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.25
    plot(dss(:),(metric2(:,2)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.5
    plot(dss(:),(metric2(:,3)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([-0.15 0.15])
    yticks([-0.15 -0.1 -0.05 0 0.05 0.1 0.15])
    ax12.YAxis.Exponent = -1;
    xlabel('$d$','Interpreter','latex')

    legend('AL','', '','AD','Interpreter','latex','Units','centimeters','Position',[8.9 4.5 1 5],'NumColumns',4,'Box','off');


% Stiff, Q_int
ax13 = axes('Units','centimeters','InnerPosition',[10.5 1 2.8 4]);
    hold(ax13,"on")
    box(ax13,"on")
    set(ax13,'FontName','Times','FontSize',10);
    set(ax13,'ColorOrder',cmap)

    metric1 = metric_var_d{1,2}; % Stiff, Q_int AL
    metric2 = metric_var_d{2,2}; % Stiff, Q_int AD
    
    % AL ---
    plot(dss(:),(metric1(:,1)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric1(:,2)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric1(:,3)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % AD ---
    plot(dss(:),(metric2(:,1)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.25
    plot(dss(:),(metric2(:,2)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.5
    plot(dss(:),(metric2(:,3)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    xticks([0 0.5 0.9])
    ylim([-1.5 1.5])
    ylabel(ax13,{'$\Delta Q$'},'Interpreter','latex','Rotation',0)
    xlabel('$d$','Interpreter','latex')

% Perm, Q_int
ax14 = axes('Units','centimeters','InnerPosition',[14.5 1 2.8 4]);
    hold(ax14,"on")
    box(ax14,"on")
    set(ax14,'FontName','Times','FontSize',10);
    set(ax14,'ColorOrder',cmap)
    
    metric1 = metric_var_d{3,2}; % Perm, Q_int AL
    metric2 = metric_var_d{4,2}; % Perm, Q_int AD
    
    % AL ---
    plot(dss(:),(metric1(:,1)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.25
    plot(dss(:),(metric1(:,2)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.5
    plot(dss(:),(metric1(:,3)-metric1(1,1)),'linestyle','-','linewidth',1) % ls = 0.75
    % AD ---
    plot(dss(:),(metric2(:,1)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.25
    plot(dss(:),(metric2(:,2)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.5
    plot(dss(:),(metric2(:,3)-metric2(1,1)),'linestyle','--','linewidth',1) % ls = 0.75
    % ------
    xlim([0 0.9])
    ylim([-1.5 1.5])
    xticks([0 0.5 0.9])
    xlabel('$d$','Interpreter','latex')    

%% Export
print('Effect-of-d&l','-dpdf','-r0')    