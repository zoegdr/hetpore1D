% Plots to compare types of damage

% In this code: 
% Swap stiff and perm, ls and lp
% Formatting (and other changes) when plotting


%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 18 14.7]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(groot,'defaultAxesTickLabelInterpreter','latex');  

% Create colourmap
maps = custom_colourmap; g = maps.green; b = maps.blue; % call custom colourmap function and unpack


%% Applied Load

% Unpack and restrict to one steady cycle
params = Perm_D_params{3};
T0 = Perm_D_Ts{1}; T1 = Perm_D_Ts{3}; T2 = Perm_D_Ts{5}; T3 = Perm_D_Ts{11};
omega = params.omega;
tmin = 18*pi/omega; tmid = 19*pi/omega; tmax = 20*pi/omega;
[~,t01] = min(abs(T0-tmin)); [~,t11] = min(abs(T1-tmin)); [~,t21] = min(abs(T2-tmin)); [~,t31] = min(abs(T3-tmin));
[~,t02] = min(abs(T0-tmax)); [~,t12] = min(abs(T1-tmax)); [~,t22] = min(abs(T2-tmax)); [~,t32] = min(abs(T3-tmax));

T0 = T0(t01:t02); dUdZ0 = Perm_D_dUdZ{1}(t01:t02,:); Q0 = Perm_D_Qs{1}(t01:t02,:); % uniform
T1 = T1(t11:t12); dUdZ1 = Perm_D_dUdZ{3}(t11:t12,:); Q1 = Perm_D_Qs{3}(t11:t12,:); % local
T2 = T2(t21:t22); dUdZ2 = Perm_D_dUdZ{5}(t21:t22,:); Q2 = Perm_D_Qs{5}(t21:t22,:); % positive
T3 = T3(t31:t32); dUdZ3 = Perm_D_dUdZ{11}(t31:t32,:); Q3 = Perm_D_Qs{11}(t31:t32,:); % negative

%% Plot

% Strain
    % Time of max strain at max damage and corresponding index in no damage
    [~,index1] = max(dUdZ1(:,params.lp*params.N)); [~,index01] = min(abs(T0-T1(index1)));
    [~,index2] = max(dUdZ2(:,1)); [~,index02] = min(abs(T0-T2(index2)));
    [~,index3] = max(dUdZ3(:,end)); [~,index03] = min(abs(T0-T3(index3)));
    
    ax11 = axes('Units','centimeters','InnerPosition',[1.9 8.4 3 4.13]);
        hold(ax11,"on")
        box(ax11,"on")
        set(ax11,'FontName','Times','FontSize',12);
    
        % Against time
        plot(T1,dUdZ1(:,params.lp*params.N)-interp1(T0,dUdZ0(:,params.lp*params.N),T1),'LineWidth',1,'Color',g(5,:))
        plot(T2,dUdZ2(:,1)-interp1(T0,dUdZ0(:,1),T2),'LineWidth',1,'Color',g(120,:),'LineStyle','-.')
        plot(T3,dUdZ3(:,end)-interp1(T0,dUdZ0(:,end),T3),'LineWidth',1,'Color',g(end,:),'LineStyle','--')
        % plot(T1,dUdZ1(:,params.lp*params.N),'LineWidth',1,'Color',g(5,:))
        % plot(T2,dUdZ2(:,1),'LineWidth',1,'Color',g(120,:),'LineStyle','-.')
        % plot(T3,dUdZ3(:,end),'LineWidth',1,'Color',g(end,:),'LineStyle','--')
        hold(gca,"off")
        
        title(ax11,{'$\frac{\partial U_d}{\partial Z}$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        % title(ax11,{'$\frac{\partial U_d}{\partial Z} - \frac{\partial U}{\partial Z}$','',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$t$','Interpreter','latex','Rotation',0);
        xlim([tmin,tmax]);
        xticks([tmin tmid tmax])
        xticklabels({'0','T/2','T'})
        ylabel({'$Z|_{M=1-d_s}$',''},'Interpreter','latex','FontSize',12);
        ylim([-0.3 0.3]); 
        yticks([-0.3 0 0.3])
        ax11.YAxis.Exponent = -1;
    
    
        lgd = legend('local ($t_1$)','+ve ($t_2$)','-ve ($t_3$)','Interpreter','latex','FontSize',12,'NumColumns',3);
        lgd.Units = 'centimeters';
        lgd.Position = [8 13.7 2.9 1];
        lgd.IconColumnWidth = 25;
        legend('boxoff')
    
    
    ax21 = axes('Units','centimeters','InnerPosition',[1.9 1.2 3 4.13]);
        hold(ax21,"on")
        box(ax21,"on")
        set(ax21,'FontName','Times','FontSize',12);
    
        % In space
        plot(Zs,dUdZ1(index1,:)-dUdZ0(index01,:),'LineWidth',1,'Color',g(5,:));
        plot(Zs,dUdZ2(index2,:)-dUdZ0(index02,:),'LineWidth',0.8,'Color',g(120,:),'LineStyle','-.');
        plot(Zs,dUdZ3(index3,:)-dUdZ0(index03,:),'LineWidth',0.8,'Color',g(end,:),'LineStyle','--');
        
        title(ax21,{'$\frac{\partial U_d}{\partial Z} - \frac{\partial U}{\partial Z}$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xticks([0 0.5 1])
        xticklabels({'0','$l_s$','1'})
        ylabel({'$t=t_i$',''},'Interpreter','latex','FontSize',12)
        ylim([-0.01 0.1]); 
        yticks([0 0.05 0.1])
        ax21.YAxis.Exponent = -1;
    
        hold(gca,"off")
    
        % Calculate area under curve
        dU_loc_al = trapz(Zs,dUdZ1(index1,:)-dUdZ0(index01,:));
        dU_pos_al = trapz(Zs,dUdZ2(index2,:)-dUdZ0(index02,:));
        dU_neg_al = trapz(Zs,dUdZ3(index3,:)-dUdZ0(index03,:));

% Flux
    % Time of max flux at max damage and corresponding index in no damage
    [~,index1] = max(Q1(:,params.lp*params.N)); [~,index01] = min(abs(T0-T1(index1)));
    [~,index2] = max(Q2(:,1)); [~,index02] = min(abs(T0-T2(index2)));
    [~,index3] = min(abs(T3-tmax)); [~,index03] = min(abs(T0-tmax)); % just plot at T since flux is zero at Z=1
    
    ax12 = axes('Units','centimeters','InnerPosition',[6.1 8.4 3 4.13]);
        hold(ax12,"on")
        box(ax12,"on")
        set(ax12,'FontName','Times','FontSize',12);
    
        % Against time
        plot(T1,Q1(:,params.lp*params.N)-interp1(T0,Q0(:,params.lp*params.N),T1),'LineWidth',1,'Color',b(5,:))
        plot(T2,Q2(:,1)-interp1(T0,Q0(:,1),T2),'LineWidth',1,'Color',b(120,:),'LineStyle','-.')
        plot(T3,Q3(:,end)-interp1(T0,Q0(:,end),T3),'LineWidth',1,'Color',b(end,:),'LineStyle','--')
        % plot(T1,Q1(:,params.lp*params.N),'LineWidth',1,'Color',b(5,:))
        % plot(T2,Q2(:,1),'LineWidth',1,'Color',b(120,:),'LineStyle','-.')
        % plot(T3,Q3(:,end),'LineWidth',1,'Color',b(end,:),'LineStyle','--')
        hold(gca,"off")
        
        title(ax12,{'$Q_d$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        % title(ax12,{'$Q_d - Q$','',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$t$','Interpreter','latex','Rotation',0);
        xlim([tmin,tmax]);
        xticks([tmin tmid tmax])
        xticklabels({'0','T/2','T'})
        ylim([-0.5 0.5])
        yticks([-0.5 0 0.5])
        ax12.YAxis.Exponent = -1;
    
    ax22 = axes('Units','centimeters','InnerPosition',[6.1 1.2 3 4.13]);
        hold(ax22,"on")
        box(ax22,"on")
        set(ax22,'FontName','Times','FontSize',12);
    
        % In space
        plot(Zs,Q1(index1,:)-Q0(index01,:),'LineWidth',1,'Color',b(5,:));
        plot(Zs,Q2(index2,:)-Q0(index02,:),'LineWidth',0.8,'Color',b(120,:),'LineStyle','-.');
        plot(Zs,Q3(index3,:)-Q0(index03,:),'LineWidth',0.8,'Color',b(end,:),'LineStyle','--');
        
        title(ax22,{'$Q_d - Q$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xticks([0 0.5 1])
        xticklabels({'0','$l_s$','1'})
        ylim([-0.25 0.25])
        yticks([-0.25 0 0.25])
        ax22.YAxis.Exponent = -1;
    
        hold(gca,"off")
    
        % Calculate area under curve
        Q_loc_al = trapz(Zs,Q1(index1,:)-Q0(index01,:));
        Q_pos_al = trapz(Zs,Q2(index2,:)-Q0(index02,:));
        Q_neg_al = trapz(Zs,Q3(index3,:)-Q0(index03,:));

%% Applied Displacement

% Unpack and restrict to one steady cycle
params = Perm_D_params{8};
T0 = Perm_D_Ts{6}; T1 = Perm_D_Ts{8}; T2 = Perm_D_Ts{10}; T3 = Perm_D_Ts{12};
omega = params.omega;
tmin = 18*pi/omega; tmid = 19*pi/omega; tmax = 20*pi/omega;
[~,t01] = min(abs(T0-tmin)); [~,t11] = min(abs(T1-tmin)); [~,t21] = min(abs(T2-tmin)); [~,t31] = min(abs(T3-tmin));
[~,t02] = min(abs(T0-tmax)); [~,t12] = min(abs(T1-tmax)); [~,t22] = min(abs(T2-tmax)); [~,t32] = min(abs(T3-tmax));

T0 = T0(t01:t02); dUdZ0 = Perm_D_dUdZ{6}(t01:t02,:); Q0 = Perm_D_Qs{6}(t01:t02,:); % uniform
T1 = T1(t11:t12); dUdZ1 = Perm_D_dUdZ{8}(t11:t12,:); Q1 = Perm_D_Qs{8}(t11:t12,:); % local
T2 = T2(t21:t22); dUdZ2 = Perm_D_dUdZ{10}(t21:t22,:); Q2 = Perm_D_Qs{10}(t21:t22,:); % positive
T3 = T3(t31:t32); dUdZ3 = Perm_D_dUdZ{12}(t31:t32,:); Q3 = Perm_D_Qs{12}(t31:t32,:); % negative

%% Plot

% Strain
    % Time of max strain at max damage and corresponding index in no damage
    [~,index1] = max(dUdZ1(:,params.lp*params.N)); [~,index01] = min(abs(T0-T1(index1)));
    [~,index2] = max(dUdZ2(:,1)); [~,index02] = min(abs(T0-T2(index2)));
    [~,index3] = max(dUdZ3(:,end)); [~,index03] = min(abs(T0-T3(index3)));
    
    
    ax13 = axes('Units','centimeters','InnerPosition',[10.5 8.4 2.9 4.13]);
        hold(ax13,"on")
        box(ax13,"on")
        set(ax13,'FontName','Times','FontSize',12);
    
        % Against time
        plot(T1,dUdZ1(:,params.lp*params.N)-interp1(T0,dUdZ0(:,params.lp*params.N),T1),'LineWidth',1,'Color',g(5,:))
        plot(T2,dUdZ2(:,1)-interp1(T0,dUdZ0(:,1),T2),'LineWidth',1,'Color',g(120,:),'LineStyle','-.')
        plot(T3,dUdZ3(:,end)-interp1(T0,dUdZ0(:,end),T3),'LineWidth',1,'Color',g(end,:),'LineStyle','--')
        % plot(T1,dUdZ1(:,params.lp*params.N),'LineWidth',1,'Color',g(5,:))
        % plot(T2,dUdZ2(:,1),'LineWidth',1,'Color',g(120,:),'LineStyle','-.')
        % plot(T3,dUdZ3(:,end),'LineWidth',1,'Color',g(end,:),'LineStyle','--')
        hold(gca,"off")
        
        title(ax13,{'$\frac{\partial U_d}{\partial Z}$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        % title(ax13,{'$\frac{\partial U_d}{\partial Z}-\frac{\partial U}{\partial Z}$','',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$t$','Interpreter','latex','Rotation',0);
        xlim([tmin,tmax]);
        xticks([tmin tmid tmax])
        xticklabels({'0','T/2','T'})
        ylim([-0.3 0.3]); 
        yticks([-0.3 0 0.3])
        ax13.YAxis.Exponent = -1;
    
    
    ax23 = axes('Units','centimeters','InnerPosition',[10.5 1.2 2.9 4.13]);
        hold(ax23,"on")
        box(ax23,"on")
        set(ax23,'FontName','Times','FontSize',12);
    
        % In space
        plot(Zs,dUdZ1(index1,:)-dUdZ0(index01,:),'LineWidth',1,'Color',g(5,:));
        plot(Zs,dUdZ2(index2,:)-dUdZ0(index02,:),'LineWidth',0.8,'Color',g(120,:),'LineStyle','-.');
        plot(Zs,dUdZ3(index3,:)-dUdZ0(index03,:),'LineWidth',0.8,'Color',g(end,:),'LineStyle','--');
        
        title(ax23,{'$\frac{\partial U_d}{\partial Z}-\frac{\partial U}{\partial Z}$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xticks([0 0.5 1])
        xticklabels({'0','$l_s$','1'})
        ylim([-0.01 0.1]); 
        yticks([0 0.05 0.1])
        ax23.YAxis.Exponent = -1;
    
        hold(gca,"off")
        
        % Calculate area under curve
        dU_loc_ad = trapz(Zs,dUdZ1(index1,:)-dUdZ0(index01,:));
        dU_pos_ad = trapz(Zs,dUdZ2(index2,:)-dUdZ0(index02,:));
        dU_neg_ad = trapz(Zs,dUdZ3(index3,:)-dUdZ0(index03,:));

% Flux
    % Time of max flux at max damage and corresponding index in no damage
    [~,index1] = max(Q1(:,params.lp*params.N)); [~,index01] = min(abs(T0-T1(index1)));
    [~,index2] = max(Q2(:,1)); [~,index02] = min(abs(T0-T2(index2)));
    [~,index3] = min(abs(T3-tmax)); [~,index03] = min(abs(T0-tmax)); % just plot at T since flux is zero at Z=1
    
    
    ax14 = axes('Units','centimeters','InnerPosition',[14.8 8.4 2.9 4.13]);
        hold(ax14,"on")
        box(ax14,"on")
        set(ax14,'FontName','Times','FontSize',12);
    
        % Against time
        plot(T1,Q1(:,params.lp*params.N)-interp1(T0,Q0(:,params.lp*params.N),T1),'LineWidth',1,'Color',b(5,:))
        plot(T2,Q2(:,1)-interp1(T0,Q0(:,1),T2),'LineWidth',1,'Color',b(120,:),'LineStyle','-.')
        plot(T3,Q3(:,end)-interp1(T0,Q0(:,end),T3),'LineWidth',1,'Color',b(end,:),'LineStyle','--')
        % plot(T1,Q1(:,params.lp*params.N),'LineWidth',1,'Color',b(5,:))
        % plot(T2,Q2(:,1),'LineWidth',1,'Color',b(120,:),'LineStyle','-.')
        % plot(T3,Q3(:,end),'LineWidth',1,'Color',b(end,:),'LineStyle','--')
        hold(gca,"off")
        
        title(ax14,{'$Q_d$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        % title(ax14,{'$Q_d - Q$','',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$t$','Interpreter','latex','Rotation',0);
        xlim([tmin,tmax]);
        xticks([tmin tmid tmax])
        xticklabels({'0','T/2','T'})
        ylim([-0.5 0.5])
        yticks([-0.5 0  0.5])
        ax14.YAxis.Exponent = -1;
    
    
    ax24 = axes('Units','centimeters','InnerPosition',[14.8 1.2 2.9 4.13]);
        hold(ax24,"on")
        box(ax24,"on")
        set(ax24,'FontName','Times','FontSize',12);
    
        % In space
        plot(Zs,Q1(index1,:)-Q0(index01,:),'LineWidth',1,'Color',b(5,:));
        plot(Zs,Q2(index2,:)-Q0(index02,:),'LineWidth',0.8,'Color',b(120,:),'LineStyle','-.');
        plot(Zs,Q3(index3,:)-Q0(index03,:),'LineWidth',0.8,'Color',b(end,:),'LineStyle','--');
        
        title(ax24,{'$Q_d - Q$',''},'FontSize',12,'FontWeight','normal','Interpreter','latex')
        xlabel('$Z$','Interpreter','latex','Rotation',0);
        xticks([0 0.5 1])
        xticklabels({'0','$l_s$','1'})
        ylim([-0.25 0.25])
        yticks([-0.25 0 0.25])
        ax24.YAxis.Exponent = -1;
    
        hold(gca,"off")
        
        % Calculate area under curve
        Q_loc_ad = trapz(Zs,Q1(index1,:)-Q0(index01,:));
        Q_pos_ad = trapz(Zs,Q2(index2,:)-Q0(index02,:));
        Q_neg_ad = trapz(Zs,Q3(index3,:)-Q0(index03,:));

%% Labels



%% Export ----------------------------------------------------------------------

% print('Compare-damage-v1','-dpdf','-loose')




