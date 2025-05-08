% Multiple profiles decreased stiffness VS Z
% Load/import Stiffness+Perm_Damages_Locations.mat

%% Choose Stiff or Perm + CHANGE NAME OF FILE AT BOTTOM ACCORDINGLY
% T = Stiff_D_Ts; dUdZ = Stiff_D_dUdZ; Q = Stiff_D_Qs; params = Stiff_D_params;
T = Perm_D_Ts; dUdZ = Perm_D_dUdZ; Q = Perm_D_Qs; params = Perm_D_params;

%% Set up

omega = params{1,1}.omega; Astar1 = params{1,1}.Astar; Astar2 = params{2,1}.Astar;

% Create empty structure for saving time points
ts = zeros(9,1);

% Plot Limits
ymin1 = -0.3; ymax1 = 0.3; % AL
ymin2 = -1; ymax2 = 1; % AD

% Colour scheme
cmap = slanCM('coolwarm',9);
cmap2 = slanCM('coolwarm');

%% PLOT ----------------------------------------------------------

% Create figure
fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'Position',[0 0 16 20.6])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[16 20.6])


% Create subfigures

% Row 1 --------------------------------------------------------------------

ax11 = axes('Units','centimeters','InnerPosition',[2 15.9 2.8 4]);
    hold(ax11,"on")
    box(ax11,"on")
    set(ax11,'FontName','Times','FontSize',12);
    set(ax11,'ColorOrder',cmap)

    % Unpack
    t = T{1,1}(:,:); y = dUdZ{1,1}(:,:);  
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    ylabel('$U_Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    
    % Inset applied load s*
    ax110 = axes('Units','centimeters','InnerPosition',[2.3 16.3 1 1]);
        hold(ax110,"on")
        box(ax110,"on")
        set(ax110,'FontName','Times','FontSize',8);
        set(ax110,'Colormap',cmap2)
        
        t = linspace(0, 2*pi/omega, 100); 
        y = params{1,1}.S_star(t);
        z = zeros(size(t));
        lineColor = t;  % This is the color, it varies with x in this case.
        surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
        t2 = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi]/omega;
        % Add time points
        plot(t2,params{1,1}.S_star(t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")

        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([0 Astar1]); 

        % Labels
        title('$s^*$','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        xticks([0 pi/omega 2*pi/omega])
        xticklabels({'','', ''})
        yticks([0 Astar1/2 Astar1])
        yticklabels({'','',''})

    % Inset ds*/dt
    ax111 = axes('Units','centimeters','InnerPosition',[3.6 16.3 1 1]);
        hold(ax111,"on")
        box(ax111,"on")
        set(ax111,'FontName','Times','FontSize',8);
        set(ax111,'Colormap',cmap2)
        
        t = linspace(0, 2*pi/omega, 100); 
        y = (Astar1*omega/2)*sin(omega*t);
        z = zeros(size(t));
        lineColor = t;  % This is the color, it varies with x in this case.
        surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
        t2 = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi]/omega;
        % Add time points
        plot(t2,(Astar1*omega/2)*sin(omega*t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")

        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-Astar1*omega/2 Astar1*omega/2]); 

        % Labels
        title('$\dot{s}^*$','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        xticks([0 pi/omega 2*pi/omega])
        xticklabels({'','', ''})
        yticks([-Astar1*omega/2 0 Astar1*omega/2])
        yticklabels({'','',''})
        
        
ax12 = axes('Units','centimeters','InnerPosition',[5.2 15.9 2.8 4]);
    hold(ax12,"on")
    box(ax12,"on")
    set(ax12,'FontName','Times','FontSize',12);
    set(ax12,'ColorOrder',cmap)

   % Unpack
    t = T{1,2}(:,:); y = dUdZ{1,2}(:,:);  
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax13 = axes('Units','centimeters','InnerPosition',[8.4 15.9 2.8 4]);
    hold(ax13,"on")
    box(ax13,"on")
    set(ax13,'FontName','Times','FontSize',12);
    set(ax13,'ColorOrder',cmap)

   % Unpack
    t = T{1,3}(:,:); y = dUdZ{1,3}(:,:);  
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax14 = axes('Units','centimeters','InnerPosition',[11.6 15.9 2.8 4]);
    hold(ax14,"on")
    box(ax14,"on")
    set(ax14,'FontName','Times','FontSize',12);
    set(ax14,'ColorOrder',cmap)

   % Unpack
    t = T{1,4}(:,:); y = dUdZ{1,4}(:,:);  
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})


% Row 2 --------------------------------------------------------------------

ax21 = axes('Units','centimeters','InnerPosition',[2 11.4 2.8 4]);
    hold(ax21,"on")
    box(ax21,"on")
    set(ax21,'FontName','Times','FontSize',12);
    set(ax21,'ColorOrder',cmap)

    % Unpack
    t = T{1,1}(:,:); y = Q{1,1}(:,:);  

    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    ylabel('$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    
    % % Inset applied load
    % ax210 = axes('Units','centimeters','InnerPosition',[3.6 11.8 1 1]);
    %     hold(ax210,"on")
    %     box(ax210,"on")
    %     set(ax210,'FontName','Times','FontSize',8);
    %     set(ax210,'Colormap',cmap2)
    % 
    %     t = linspace(0, 2*pi/omega, 100);
    %     y = (Astar1*omega/2)*sin(omega*t);
    %     z = zeros(size(t));
    %     lineColor = t;  % This is the color, it varies with x in this case.
    %     surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
    %     % Add time points
    %     plot(t2,(Astar1*omega/2)*sin(omega*t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")
    % 
    %     set(gca, "Layer", "top")
    %     hold(gca,'off')
    % 
    %     % Limits
    %     ylim([-Astar1*omega/2 Astar1*omega/2]); 
    % 
    %     % Labels
    %     title('$\dot{s}^*$','Interpreter','latex')
    %     xlabel('$t$','Interpreter','latex')
    %     xticks([0 pi/omega 2*pi/omega])
    %     xticklabels({'','', ''})
    %     yticks([-Astar1*omega/2 0 Astar1*omega/2])
    %     yticklabels({'','',''})

ax22 = axes('Units','centimeters','InnerPosition',[5.2 11.4 2.8 4]);
    hold(ax22,"on")
    box(ax22,"on")
    set(ax22,'FontName','Times','FontSize',12);
    set(ax22,'ColorOrder',cmap)

    % Unpack
    t = T{1,2}(:,:); y = Q{1,2}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax23 = axes('Units','centimeters','InnerPosition',[8.4 11.4 2.8 4]);
    hold(ax23,"on")
    box(ax23,"on")
    set(ax23,'FontName','Times','FontSize',12);
    set(ax23,'ColorOrder',cmap)

    % Unpack
    t = T{1,3}(:,:); y = Q{1,3}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax24 = axes('Units','centimeters','InnerPosition',[11.6 11.4 2.8 4]);
    hold(ax24,"on")
    box(ax24,"on")
    set(ax24,'FontName','Times','FontSize',12);
    set(ax24,'ColorOrder',cmap)

    % Unpack
    t = T{1,4}(:,:); y = Q{1,4}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

% Row 3 --------------------------------------------------------------------

ax31 = axes('Units','centimeters','InnerPosition',[2 5.7 2.8 4]);
    hold(ax31,"on")
    box(ax31,"on")
    set(ax31,'FontName','Times','FontSize',12);
    set(ax31,'ColorOrder',cmap)

    % Unpack
    t = T{2,1}(:,:); y = dUdZ{2,1}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    ylabel('$U_Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
        
    % Inset applied load a
    ax310 = axes('Units','centimeters','InnerPosition',[2.3 6.1 1 1]);
        hold(ax310,"on")
        box(ax310,"on")
        set(ax310,'FontName','Times','FontSize',8);
        set(ax310,'Colormap',cmap2)
        
        t = linspace(0, 2*pi/omega, 100); 
        y = params{2,1}.a(t);
        z = zeros(size(t));
        lineColor = t;  % This is the color, it varies with x in this case.
        surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
        t2 = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi]/omega;
        % Add time points
        plot(t2,params{2,1}.a(t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")

        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-Astar2 0]); 

        % Labels
        title('$a$','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        xticks([0 pi/omega 2*pi/omega])
        xticklabels({'','', ''})
        yticks([-Astar2 -Astar2/2 0])
        yticklabels({'','',''})

    % Inset da/dt
    ax311 = axes('Units','centimeters','InnerPosition',[3.6 6.1 1 1]);
        hold(ax311,"on")
        box(ax311,"on")
        set(ax311,'FontName','Times','FontSize',8);
        set(ax311,'Colormap',cmap2)
        
        t = linspace(0, 2*pi/omega, 100); 
        y = params{2,1}.adot(t);
        z = zeros(size(t));
        lineColor = t;  % This is the color, it varies with x in this case.
        surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
        t2 = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi]/omega;
        % Add time points
        plot(t2,params{2,1}.adot(t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")

        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-Astar2*omega/2 Astar2*omega/2]); 

        % Labels
        % Labels
        title('$\dot{a}$','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        xticks([0 pi/omega 2*pi/omega])
        xticklabels({'','', ''})
        yticks([-Astar2*omega/2 0 Astar2*omega/2])
        yticklabels({'','',''})

ax32 = axes('Units','centimeters','InnerPosition',[5.2 5.7 2.8 4]);
    hold(ax32,"on")
    box(ax32,"on")
    set(ax32,'FontName','Times','FontSize',12);
    set(ax32,'ColorOrder',cmap)

   % Unpack
    t = T{2,2}(:,:); y = dUdZ{2,2}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax33 = axes('Units','centimeters','InnerPosition',[8.4 5.7 2.8 4]);
    hold(ax33,"on")
    box(ax33,"on")
    set(ax33,'FontName','Times','FontSize',12);
    set(ax33,'ColorOrder',cmap)

   % Unpack
    t = T{2,3}(:,:); y = dUdZ{2,3}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax34 = axes('Units','centimeters','InnerPosition',[11.6 5.7 2.8 4]);
    hold(ax34,"on")
    box(ax34,"on")
    set(ax34,'FontName','Times','FontSize',12);
    set(ax34,'ColorOrder',cmap)

   % Unpack
    t = T{2,4}(:,:); y = dUdZ{2,4}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([0 1]);

    % Labels
    xticks([0 0.5 1])
    xticklabels({'','',''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

% Row 4 --------------------------------------------------------------------

ax41 = axes('Units','centimeters','InnerPosition',[2 1.2 2.8 4]);
    hold(ax41,"on")
    box(ax41,"on")
    set(ax41,'FontName','Times','FontSize',12);
    set(ax41,'ColorOrder',cmap)

    % Unpack
    t = T{2,1}(:,:); y = Q{2,1}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    ylabel('$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])

    % % Inset applied disp da/dt
    % ax410 = axes('Units','centimeters','InnerPosition',[3.6 1.6 1 1]);
    %     hold(ax410,"on")
    %     box(ax410,"on")
    %     set(ax410,'FontName','Times','FontSize',8);
    %     set(ax410,'Colormap',cmap2)
    % 
    %     t = linspace(0, 2*pi/omega, 100);
    %     y = params{2,1}.adot(t);
    %     z = zeros(size(t));
    %     lineColor = t;  % This is the color, it varies with x in this case.
    %     surface([t;t], [y;y], [z;z], [lineColor;lineColor],'FaceColor', 'no','EdgeColor', 'interp','LineWidth', 1);
    %     % Add time points
    %     plot(t2,params{2,1}.adot(t2),'LineStyle','None','Marker','.','MarkerSize',5,'Color',"k")
    % 
    %     set(gca, "Layer", "top")
    %     hold(gca,'off')
    % 
    %     % Limits
    %     ylim([-Astar2*omega/2 Astar2*omega/2]); 
    % 
    %     % Labels
    %     % Labels
    %     title('$\dot{a}$','Interpreter','latex')
    %     xlabel('$t$','Interpreter','latex')
    %     xticks([0 pi/omega 2*pi/omega])
    %     xticklabels({'','', ''})
    %     yticks([-Astar2*omega/2 0 Astar2*omega/2])
    %     yticklabels({'','',''})

ax42 = axes('Units','centimeters','InnerPosition',[5.2 1.2 2.8 4]);
    hold(ax42,"on")
    box(ax42,"on")
    set(ax42,'FontName','Times','FontSize',12);
    set(ax42,'ColorOrder',cmap)

    % Unpack
    t = T{2,2}(:,:); y = Q{2,2}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax43 = axes('Units','centimeters','InnerPosition',[8.4 1.2 2.8 4]);
    hold(ax43,"on")
    box(ax43,"on")
    set(ax43,'FontName','Times','FontSize',12);
    set(ax43,'ColorOrder',cmap)

    % Unpack
    t = T{2,3}(:,:); y = Q{2,3}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax44 = axes('Units','centimeters','InnerPosition',[11.6 1.2 2.8 4]);
    hold(ax44,"on")
    box(ax44,"on")
    set(ax44,'FontName','Times','FontSize',12);
    set(ax44,'ColorOrder',cmap)

    % Unpack
    t = T{2,4}(:,:); y = Q{2,4}(:,:);  

    % % Shading (deloading)
    % fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    for i = 0:7
        [~,ts(i+1)] = min(abs(t-(19*2*pi+i*pi/4)/omega));
        plot(Zs,y(ts(i+1),:),'Linewidth',1)
    end
    [~,ts(9)] = min(abs(t-(20*2*pi)/omega));
    plot(Zs,y(ts(9),:),'Linewidth',1.5,'LineStyle',':')

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([0 1]);

    % Labels
    xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([0 0.5 1])
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

%% Export ----------------------------------------------------------------------

% print('Stiff-damage-profiles-VS-Z','-dpdf','-r0')
print('Perm-damage-profiles-VS-Z','-dpdf','-r0')