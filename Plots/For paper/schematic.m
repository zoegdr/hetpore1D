% Load/import Stiffness+Perm_Damages_Locations.mat
T1 = Stiff_D_Ts{2,1}; dUdZ1 = Stiff_D_dUdZ{2,1}; Q1 = Stiff_D_Qs{2,1};
T2 = Stiff_D_Ts{2,3}; dUdZ2 = Stiff_D_dUdZ{2,3}; Q2 = Stiff_D_Qs{2,3};

% Set stiffness
v = 0.1; % variance
M = @(Z) 1-0.35*exp(-(Z-0.5).^2/(2*v^2)); % stiffness
L = @(Z) 0.3/0.7 * ( 1-0.35*exp(-(Z-0.5).^2/(2*v^2)) );

% Calculate stress
S1 = zeros(size(dUdZ1)); S2 = zeros(size(dUdZ2));
for i = 1:length(Zs)
    S1(:,i) = 0.5*1 .*( 1+ dUdZ1(:,i)-1./(1+dUdZ1(:,i)) ) + 0.5*1 .*( 1 + dUdZ1(:,i) +1./(1+dUdZ1(:,1)) -2 );
    S2(:,i) = 0.5*M(Zs(i)) .*( 1+ dUdZ2(:,i)-1./(1+dUdZ2(:,i)) ) + 0.5*L(Zs(i)) .*( 1 + dUdZ2(:,i) +1./(1+dUdZ2(:,1)) -2 );
end

%% Set up

% Colour scheme
c = slanCM('tab20c');
cmap = [c(92,:);c(52,:)];

%% PLOT

% Create figure ---------------------------

fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
set(fig,'Position',[0 0 9 3.5]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[9 3.5])

ax1 = axes('Units','centimeters','InnerPosition',[3 0.5 2 2]);
        hold(ax1,"on")
        box(ax1,"off")
        set(ax1,'FontName','Times','FontSize',10);
        set(ax1,'ColorOrder',cmap)

        plot(Zs,trapz(T1,S1,1),'LineWidth',1)
        plot(Zs,trapz(T2,S2,1),'LineWidth',1)

        xlim([0 1]);

        xticks([0 0.4 0.5 0.6 1])
        xticklabels({'','','$l$','',''})
        yticks([])

        % ylabel('Stress','FontSize',10)
        ylabel('$\overline{s''}$','Interpreter','latex','Rotation',0,'FontSize',10);

legend('Undamaged','Damaged','Interpreter','latex','Units','centimeters','Position',[4 0.75 1 5],'NumColumns',2,'Box','off','IconColumnWidth',15);

ax2 = axes('Units','centimeters','InnerPosition',[6 0.5 2 2]);
        hold(ax2,"on")
        box(ax2,"off")
        set(ax2,'FontName','Times','FontSize',10);
        set(ax2,'ColorOrder',cmap)
        
        plot(Zs,trapz(T1,Q1,1),'LineWidth',1)
        plot(Zs,trapz(T2,Q2,1),'LineWidth',1)
        
        xlim([0 1]);

        xticks([0 0.4 0.5 0.6 1])
        xticklabels({'','','$l$','',''})
        yticks([])

        % ylabel('Flux','FontSize',10)
        ylabel('$\overline{Q}$','Interpreter','latex','Rotation',0,'FontSize',10);

print(fig,'Stress-flux-schematic','-dpdf','-r0')


