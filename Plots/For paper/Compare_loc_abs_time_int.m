% Integrate over all time (not just one cycle), compare at different
% locations
% Load/import Stiffness+Perm_Damages_Locations.mat or
% Inc_Stiff_3_locations.m

%% Choose Stiff or Perm or Inc + CHANGE NAME OF FILE AT BOTTOM ACCORDINGLY
% T = Stiff_D_Ts; dUdZ = Stiff_D_dUdZ; Q = Stiff_D_Qs;
% T = Perm_D_Ts; dUdZ = Perm_D_dUdZ; Q = Perm_D_Qs;
T = Inc_T; dUdZ = Inc_dUdZ; Q = Inc_Q;

%% Set up

% Colour scheme
c = slanCM('tab20c');
cmap = [c(52,:);c(78,:);c(92,:);c(205,:)];

%% PLOT

% Create figure ---------------------------

fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
set(fig,'Position',[0 0 9 13]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[9 13])

% Applied Load --------------------------

% Unpack + take absolute value
T0 = T{1,1}; dUdZ0 = abs(dUdZ{1,1}); Q0 = abs(Q{1,1}); % uniform
T1 = T{1,2}; dUdZ1 = abs(dUdZ{1,2}); Q1 = abs(Q{1,2}); % ls = 0.25
T2 = T{1,3}; dUdZ2 = abs(dUdZ{1,3}); Q2 = abs(Q{1,3}); % ls = 0.5
T3 = T{1,4}; dUdZ3 = abs(dUdZ{1,4}); Q3 = abs(Q{1,4}); % ls = 0.75
% No absolute value
% T0 = Stiff_D_Ts{1,1}; dUdZ0 = Stiff_D_dUdZ{1,1}; Q0 = Stiff_D_Qs{1,1}; % uniform
% T1 = Stiff_D_Ts{1,2}; dUdZ1 = Stiff_D_dUdZ{1,2}; Q1 = Stiff_D_Qs{1,2}; % ls = 0.2
% T2 = Stiff_D_Ts{1,3}; dUdZ2 = Stiff_D_dUdZ{1,3}; Q2 = Stiff_D_Qs{1,3}; % ls = 0.5
% T3 = Stiff_D_Ts{1,4}; dUdZ3 = Stiff_D_dUdZ{1,4}; Q3 = Stiff_D_Qs{1,4}; % ls = 0.8

% Strain    
    ax11 = axes('Units','centimeters','InnerPosition',[1.25 7 2.8 4]);
        hold(ax11,"on")
        box(ax11,"on")
        set(ax11,'FontName','Times','FontSize',10);
        set(ax11,'ColorOrder',cmap)
        
        plot(Zs,trapz(T1,dUdZ1,1),'LineWidth',1)
        plot(Zs,trapz(T2,dUdZ2,1),'LineWidth',1)
        plot(Zs,trapz(T3,dUdZ3,1),'LineWidth',1)
        plot(Zs,trapz(T0,dUdZ0,1),'LineWidth',1,'LineStyle',':')
        hold(gca,"off")
        
        title(ax11,{'$\overline{|U_Z|} $',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlim([0 1]);
        xticks([0.25 0.5 0.75])

        % Stiffness Damage
        ylim([0.5 2.2]); 
        yticks([0.5 1 1.5 2])

legend('$l = 0.25$','$l = 0.5$','$l=0.75$','$d = 0$','Interpreter','latex','Units','centimeters','Position',[4 10 1 5],'NumColumns',4,'Box','off','IconColumnWidth',15);



% Flux
    ax12 = axes('Units','centimeters','InnerPosition',[5.5 7 2.8 4]);
        hold(ax12,"on")
        box(ax12,"on")
        set(ax12,'FontName','Times','FontSize',10);
        set(ax12,'ColorOrder',cmap)
        
        plot(Zs,trapz(T1,Q1,1),'LineWidth',1)
        plot(Zs,trapz(T2,Q2,1),'LineWidth',1)
        plot(Zs,trapz(T3,Q3,1),'LineWidth',1)
        plot(Zs,trapz(T0,Q0,1),'LineWidth',1,'LineStyle',':')
        hold(gca,"off")
        
        title(ax12,{'$\overline{|Q|}$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlim([0 1]);
        xticks([0.25 0.5 0.75])

        % Stiffness Damage
        ylim([0 6])
        yticks([0 2 4 6])

% % Calculate area under curve
dU1 = trapz(Zs,trapz(T1,dUdZ1,1));
dU2 = trapz(Zs,trapz(T2,dUdZ2,1));
dU3 = trapz(Zs,trapz(T3,dUdZ3,1));

% Applied Displacement ---------------------------------------------

% Unpack + take absolute value
T0 = T{2,1}; dUdZ0 = abs(dUdZ{2,1}); Q0 = abs(Q{2,1}); % uniform
T1 = T{2,2}; dUdZ1 = abs(dUdZ{2,2}); Q1 = abs(Q{2,2}); % ls = 0.25
T2 = T{2,3}; dUdZ2 = abs(dUdZ{2,3}); Q2 = abs(Q{2,3}); % ls = 0.5
T3 = T{2,4}; dUdZ3 = abs(dUdZ{2,4}); Q3 = abs(Q{2,4}); % ls = 0.75
% No absolute value
% T0 = Stiff_D_Ts{2,1}; dUdZ0 = Stiff_D_dUdZ{2,1}; Q0 = Stiff_D_Qs{2,1}; % uniform
% T1 = Stiff_D_Ts{2,2}; dUdZ1 = Stiff_D_dUdZ{2,2}; Q1 = Stiff_D_Qs{2,2}; % ls = 0.2
% T2 = Stiff_D_Ts{2,3}; dUdZ2 = Stiff_D_dUdZ{2,3}; Q2 = Stiff_D_Qs{2,3}; % ls = 0.5
% T3 = Stiff_D_Ts{2,4}; dUdZ3 = Stiff_D_dUdZ{2,4}; Q3 = Stiff_D_Qs{2,4}; % ls = 0.8

% Strain  
    ax21 = axes('Units','centimeters','InnerPosition',[1.25 1.2 2.8 4]);
        hold(ax21,"on")
        box(ax21,"on")
        set(ax21,'FontName','Times','FontSize',10);
        set(ax21,'ColorOrder',cmap)

        plot(Zs,trapz(T1,dUdZ1,1),'LineWidth',1)
        plot(Zs,trapz(T2,dUdZ2,1),'LineWidth',1)
        plot(Zs,trapz(T3,dUdZ3,1),'LineWidth',1)
        plot(Zs,trapz(T0,dUdZ0,1),'LineWidth',1,'LineStyle',':')
        
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xlim([0 1])
        xticks([0.25 0.5 0.75])

        % Stiffness Damage
        ylim([0.5 2.2]); 
        yticks([0.5 1 1.5 2])
    
        hold(gca,"off")


% Flux
    ax22 = axes('Units','centimeters','InnerPosition',[5.5 1.2 2.8 4]);
        hold(ax22,"on")
        box(ax22,"on")
        set(ax22,'FontName','Times','FontSize',10);
        set(ax22,'ColorOrder',cmap)

        plot(Zs,trapz(T1,Q1,1),'LineWidth',1)
        plot(Zs,trapz(T2,Q2,1),'LineWidth',1)
        plot(Zs,trapz(T3,Q3,1),'LineWidth',1)
        plot(Zs,trapz(T0,Q0,1),'LineWidth',1,'LineStyle',':')
        
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xlim([0 1])
        xticks([0.25 0.5 0.75])

        % Stiffness Damage
        ylim([0 6])
        yticks([0 2 4 6])
    
        hold(gca,"off")
    
% % Calculate area under curve
Q1 = trapz(Zs,trapz(T1,Q1,1));
Q2 = trapz(Zs,trapz(T2,Q2,1));
Q3 = trapz(Zs,trapz(T3,Q3,1));


%% Save file
% print(fig,'Compare-abs-Stiff-damage-Locations','-dpdf','-r0')
% print(fig,'Compare-abs-Perm-damage-Locations','-dpdf','-r0')
print(fig,'Compare-abs-Stiff-INC-damage-Locations','-dpdf','-r0')
