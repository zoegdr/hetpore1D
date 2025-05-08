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
set(fig,'Position',[0 0 18 12])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[18,12])


% Create subfigures

% Row 1 (AL) --------------------------------------------------------------------

% Stiff, dU_int
ax11 = axes('Units','centimeters','InnerPosition',[1.7 7 4 2.5]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',10);
    set(ax11,'ColorOrder',cmap)

    metric = metric_var_freq{1,1};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.1 10.5])
    % yticks([0 1 2 3 4])
    ylabel('$\Delta U$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')
    
    legend('$l = 0.25$','$l = 0.5$','$l=0.75$','Interpreter','latex','Units','centimeters','Position',[8.2 9 1 5],'NumColumns',4,'Box','off','IconColumnWidth',15);

    % Inset
    ax110 = axes('Units','centimeters','InnerPosition',[3.5 8 1.3 1]);
        hold(ax110,"on")
        box(ax110,"on")
        set(ax110,'FontName','Times','FontSize',7);
        set(ax110,'ColorOrder',cmap)

        plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
        xlim([20 25])
        xticks([20 25])
        ylim([-0.01 0.13])
        yticks([0 0.1])

% Perm, dU_int
ax12 = axes('Units','centimeters','InnerPosition',[7.5 7 4 2.5]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',10);
    set(ax12,'ColorOrder',cmap)

    metric = metric_var_freq{3,1};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.01 0.04])
    % yticks([-0.01 0 0.02 0.04])
    ax12.YAxis.Exponent = -2;
    ylabel('$\Delta U$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')

% Perm, Q_int    
ax13 = axes('Units','centimeters','InnerPosition',[13.5 7 4 2.5]);
    hold(ax13,"on")
    box(ax13,"on")
    set(ax13,'FontName','Times','FontSize',10);
    set(ax13,'ColorOrder',cmap)

    metric = metric_var_freq{3,2};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    ylim([-0.3 0])
    xticks([0.5 10 20 30 40 50])
    % yticks([-0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0])
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')

% Row 2 (AD) --------------------------------------------------------------------

% Stiff, Q_int
ax21 = axes('Units','centimeters','InnerPosition',[1.7 1.5 4 2.5]);
    hold(ax21,"on")
    box(ax21,"on")
    set(ax21,'FontName','Times','FontSize',10);
    set(ax21,'ColorOrder',cmap)

    metric = metric_var_freq{2,2};

    plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
    plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.1 10.5])
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')
    
    % Inset 1
    ax211 = axes('Units','centimeters','InnerPosition',[4.2 2.5 1.3 1]);
        hold(ax211,"on")
        box(ax211,"on")
        set(ax211,'FontName','Times','FontSize',7);
        set(ax211,'ColorOrder',cmap)

        plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
        xlim([20 25])
        xticks([20 25])
        ylim([-0.01 0.13])
        yticks([0 0.1])

    % Inset 2
    ax212 = axes('Units','centimeters','InnerPosition',[2.75 2.5 0.75 1]);
        hold(ax212,"on")
        box(ax212,"on")
        set(ax212,'FontName','Times','FontSize',7);
        set(ax212,'ColorOrder',cmap)

        plot(omegas,(metric(:,2)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,3)-metric(:,1)),'linestyle','-','linewidth',1)
        plot(omegas,(metric(:,4)-metric(:,1)),'linestyle','-','linewidth',1)
        xlim([0.5 0.55])
        xticks([0.5 0.55])
        % ylim([-0.01 0.13])
        yticks([9 10])

% Perm, dU_int
ax22 = axes('Units','centimeters','InnerPosition',[7.5 1.5 4 2.5]);
    hold(ax22,"on")
    box(ax22,"on")
    set(ax22,'FontName','Times','FontSize',10);
    set(ax22,'ColorOrder',cmap)

    metric = metric_var_freq{4,1};

    plot(omegas(1:65),(metric(1:65,2)-metric(1:65,1)),'linestyle','-','linewidth',1)
    plot(omegas(1:65),(metric(1:65,3)-metric(1:65,1)),'linestyle','-','linewidth',1)
    plot(omegas(1:65),(metric(1:65,4)-metric(1:65,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.01 0.04])
    % yticks([-0.01 0 0.02 0.04])
    ax22.YAxis.Exponent = -2;
    ylabel('$\Delta U$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')

% Perm, Q_int    
ax23 = axes('Units','centimeters','InnerPosition',[13.5 1.5 4 2.5]);
    hold(ax23,"on")
    box(ax23,"on")
    set(ax23,'FontName','Times','FontSize',10);
    set(ax23,'ColorOrder',cmap)

    metric = metric_var_freq{4,2};

    plot(omegas(1:65),(metric(1:65,2)-metric(1:65,1)),'linestyle','-','linewidth',1)
    plot(omegas(1:65),(metric(1:65,3)-metric(1:65,1)),'linestyle','-','linewidth',1)
    plot(omegas(1:65),(metric(1:65,4)-metric(1:65,1)),'linestyle','-','linewidth',1)
    xlim([0.5 50])
    xticks([0.5 10 20 30 40 50])
    ylim([-0.3 0])
    % yticks([-0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0])
    ylabel('$\Delta Q$','Interpreter','latex','Rotation',0)
    xlabel('$\omega$','Interpreter','latex')

%% Export
print('Effect-of-freq&l','-dpdf','-r0')
  