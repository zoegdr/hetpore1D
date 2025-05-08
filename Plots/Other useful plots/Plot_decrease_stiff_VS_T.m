% Multiple profiles decreased stiffness VS time
% Load/import Stiffness+Perm_Damages_Locations.mat

%% Set limits

% Frequency of applied load
omega = Stiff_D_params{1,1}.omega;
% Time limits
tmin = 18*pi/omega;
tmid = 19*pi/omega;
tmax = 20*pi/omega;

% Plot Limits
ymin1 = -0.3; ymax1 = 0.3; % AL
ymin2 = -1; ymax2 = 1; % AD


%% Load colourmaps  -------------------------------

% Create colourmap
maps = custom_colourmap; % call custom colourmap function
g = maps.green;
b = maps.blue;

g11 = [g(1,:);g(23,:);g(46,:);g(69,:);g(92,:);g(115,:);g(133,:);g(156,:);g(179,:);g(202,:);g(end,:)];
b11 = [b(1,:);b(23,:);b(46,:);b(69,:);b(92,:);b(115,:);b(133,:);b(156,:);b(179,:);b(202,:);b(end,:)];


%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(gcf,'Position',[0 0 18 20.6])

%% Background boxes 
% better to do in post

% % AL (top two rows)
% ax01 = axes('Units','centimeters','InnerPosition',[0.25 10.6 17.5 9.55]);
% xticks([]); yticks([]);
% box(ax01,"on"); ax01.XColor = [.8,.8,.8]; ax01.YColor = [.8,.8,.8]; ax01.LineWidth = 2;
% text(9,9.75,"AL",Units="centimeters",FontSize=12,FontName='Times'); % annotation
% 
% % AL (bottom two rows)
% ax02 = axes('Units','centimeters','InnerPosition',[0.25 0.5 17.5 9.5]);
% xticks([]); yticks([]);
% box(ax02,"on"); ax02.XColor = [.8,.8,.8]; ax02.YColor = [.8,.8,.8]; ax02.LineWidth = 2;
% text(9,9.75,"AD",Units="centimeters",FontSize=12,FontName='Times'); % annotation



%% Create subfigures

% Row 1 --------------------------------------------------------------------

ax11 = axes('Units','centimeters','InnerPosition',[2 15.9 2.8 4]);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,1}(:,:); z = Stiff_D_dUdZ{1,1}(:,:); params = Stiff_D_params{1,1};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    
    % Inset applied load
    ax110 = axes('Units','centimeters','InnerPosition',[2.1 16.1 0.75 0.75]);
        hold(gca,"on")
        box(gca,"on")
        set(gca,'FontName','Times','FontSize',8);

        % Shading (deloading)
        fill([tmid tmax tmax tmid], [0 0 params.Astar params.Astar], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])

        % Plot
        plot(x(index1:index3),params.S_star(x(index1:index3)),'LineWidth',0.75,'Color','k');
        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([0 params.Astar]); 
        xlim([tmin,tmax]);
    
        % Labels
        title('$s^*$','Interpreter','latex')
        xticks([tmin tmid tmax])
        xticklabels({'','', ''})
        yticks([0 params.Astar/2 params.Astar])
        yticklabels({'','',''})
        
        
ax12 = axes('Units','centimeters','InnerPosition',[5.2 15.9 2.8 4]);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,2}(:,:); z = Stiff_D_dUdZ{1,2}(:,:); params = Stiff_D_params{1,2};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);


    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax13 = axes('Units','centimeters','InnerPosition',[8.4 15.9 2.8 4]);
    set(ax13,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,3}; z = Stiff_D_dUdZ{1,3}; params = Stiff_D_params{1,3};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);

    set(gca, "Layer", "top")
    hold(gca,'off')    
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax14 = axes('Units','centimeters','InnerPosition',[11.6 15.9 2.8 4]);
    set(ax14,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,4}(:,:); z = Stiff_D_dUdZ{1,4}(:,:); params = Stiff_D_params{1,4};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);

    set(gca, "Layer", "top")
    hold(gca,'off')    

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

% Gradient, old version
% ax15 = axes('Units','centimeters','InnerPosition',[14.8 15.9 2.8 4]);
%     set(ax14,'DefaultAxesColorOrder',g11);
%     hold(gca,"on")
%     box(gca,"on")
%     set(gca,'FontName','Times','FontSize',12);
% 
%     % Unpack
%     x = Stiff_D_Ts{5}(:,:); z = Stiff_D_dUdZ{5}(:,:); params = Stiff_D_params{5};
% 
%     % Find time indexes
%     [~,index1] = min(abs(x-tmin));
%     [~,index3] = min(abs(x-tmax));
% 
%     % Shading (deloading)
%     fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
% 
%     % Plot 10 equally spaced lines over one cycle between 0<Z<1
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
%     plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
%     plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
%     plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
%     plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
%     plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
%     plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
%     plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
%     plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
%     plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
%     plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
% 
%     % emphasize points around and at damage
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',[0.97 0.29 0.11]);
% 
%     set(gca, "Layer", "top")
%     hold(gca,'off')    
% 
%     % Limits
%     ylim([ymin1 ymax1]); 
%     xlim([tmin,tmax]);
% 
%     % Labels
%     xticks([tmin tmid tmax])
%     xticklabels({'','', ''})
%     yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
%     yticklabels({'','',''})




% Row 2 --------------------------------------------------------------------

ax21 = axes('Units','centimeters','InnerPosition',[2 11.4 2.8 4]);
    set(ax21,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,1}(:,:); z = Stiff_D_Qs{1,1}(:,:); params = Stiff_D_params{1,1};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])

    % Inset applied load
    ax210 = axes('Units','centimeters','InnerPosition',[2.1 11.6 0.75 0.75]);
        hold(gca,"on")
        box(gca,"on")
        set(gca,'FontName','Times','FontSize',8);

        % Shading (deloading)
        fill([tmid tmax tmax tmid], [-params.Astar*params.omega/2 -params.Astar*params.omega/2 params.Astar*params.omega/2 params.Astar*params.omega/2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])

        % Plot
        S_dot = @(t) (params.Astar*params.omega/2)*sin(params.omega*t);
        plot(x(index1:index3),S_dot(x(index1:index3)),'LineWidth',0.75,'Color','k');
        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-params.Astar*params.omega/2 params.Astar*params.omega/2]); 
        xlim([tmin,tmax]);
    
        % Labels
        title('$\dot{s}^*$','Interpreter','latex')
        xticks([tmin tmid tmax])
        xticklabels({'','', ''})
        yticks([-params.Astar*params.omega/2 0 params.Astar*params.omega/2])
        yticklabels({'','',''})

ax22 = axes('Units','centimeters','InnerPosition',[5.2 11.4 2.8 4]);
    set(ax22,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,2}(:,:); z = Stiff_D_Qs{1,2}(:,:); params = Stiff_D_params{1,2};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax23 = axes('Units','centimeters','InnerPosition',[8.4 11.4 2.8 4]);
    set(ax23,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,3}(:,:); z = Stiff_D_Qs{1,3}(:,:); params = Stiff_D_params{1,3};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off') 

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax24 = axes('Units','centimeters','InnerPosition',[11.6 11.4 2.8 4]);
    set(ax24,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{1,4}(:,:); z = Stiff_D_Qs{1,4}(:,:); params = Stiff_D_params{1,4};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

% Gradient - old version
% ax25 = axes('Units','centimeters','InnerPosition',[14.8 11.4 2.8 4]);
%     set(ax24,'DefaultAxesColorOrder',b11);
%     hold(gca,"on")
%     box(gca,"on")
%     set(gca,'FontName','Times','FontSize',12);
% 
%     % Unpack
%     x = Stiff_D_Ts{5}(:,:); z = Stiff_D_Qs{5}(:,:); params = Stiff_D_params{5};
% 
%     % Find time indexes
%     [~,index1] = min(abs(x-tmin));
%     [~,index3] = min(abs(x-tmax));
% 
%     % Shading (deloading)
%     fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
% 
%     % Plot 10 equally spaced lines over one cycle between 0<Z<1
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
%     plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
%     plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
%     plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
%     plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
%     plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
%     plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
%     plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
%     plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
%     plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
%     plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
% 
%     % emphasize points around and at damage
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',[0.97 0.29 0.11]);
% 
%     set(gca, "Layer", "top")
%     hold(gca,'off')
% 
%     % Limits
%     ylim([ymin2 ymax2]); 
%     xlim([tmin,tmax]);
% 
%     % Labels
%     xticks([tmin tmid tmax])
%     xticklabels({'0','T/2', 'T'})
%     yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
%     yticklabels({'','',''})

% Row 3 --------------------------------------------------------------------

ax31 = axes('Units','centimeters','InnerPosition',[2 5.7 2.8 4]);
    set(ax31,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,1}(:,:); z = Stiff_D_dUdZ{2,1}(:,:); params = Stiff_D_params{2,1};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])

    % Inset applied load
    ax310 = axes('Units','centimeters','InnerPosition',[2.1 5.9 0.75 0.75]);
        hold(gca,"on")
        box(gca,"on")
        set(gca,'FontName','Times','FontSize',8);

        % Shading (deloading)
        fill([tmid tmax tmax tmid], [-params.Astar -params.Astar 0 0], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])

        % Plot
        plot(x(index1:index3),params.a(x(index1:index3)),'LineWidth',0.75,'Color','k');
        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-params.Astar 0]); 
        xlim([tmin,tmax]);
    
        % Labels
        title('$a^*$','Interpreter','latex')
        xticks([tmin tmid tmax])
        xticklabels({'','', ''})
        yticks([-params.Astar -params.Astar/2 0])
        yticklabels({'','',''})

ax32 = axes('Units','centimeters','InnerPosition',[5.2 5.7 2.8 4]);
    set(ax32,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,2}(:,:); z = Stiff_D_dUdZ{2,2}(:,:); params = Stiff_D_params{2,2};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax33 = axes('Units','centimeters','InnerPosition',[8.4 5.7 2.8 4]);
    set(ax33,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,3}(:,:); z = Stiff_D_dUdZ{2,3}(:,:); params = Stiff_D_params{2,3};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);

    set(gca, "Layer", "top")
    hold(gca,'off')   

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

ax34 = axes('Units','centimeters','InnerPosition',[11.6 5.7 2.8 4]);
    set(ax34,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,4}(:,:); z = Stiff_D_dUdZ{2,4}(:,:); params = Stiff_D_params{2,4};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);

    set(gca, "Layer", "top")
    hold(gca,'off')    

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
    yticklabels({'','',''})

% Gradient - old version
% ax35 = axes('Units','centimeters','InnerPosition',[14.8 5.7 2.8 4]);
%     set(ax34,'DefaultAxesColorOrder',g11);
%     hold(gca,"on")
%     box(gca,"on")
%     set(gca,'FontName','Times','FontSize',12);
% 
%     % Unpack
%     x = Stiff_D_Ts{10}(:,:); z = Stiff_D_dUdZ{10}(:,:); params = Stiff_D_params{10};
% 
%     % Find time indexes
%     [~,index1] = min(abs(x-tmin));
%     [~,index3] = min(abs(x-tmax));
% 
%     % Shading (deloading)
%     fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
% 
%     % Plot 10 equally spaced lines over one cycle between 0<Z<1
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',g(1,:));
%     plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',g(23,:)); % Z = 0.1
%     plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',g(46,:)); % Z = 0.2
%     plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',g(69,:)); % Z = 0.3
%     plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',g(92,:)); % Z = 0.4
%     plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',g(115,:)); % Z = 0.5
%     plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',g(133,:)); % Z = 0.6
%     plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',g(156,:)); % Z = 0.7
%     plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',g(179,:)); % Z = 0.8
%     plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',g(202,:)); % Z = 0.9
%     plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',g(end,:));
% 
%     % emphasize points around and at damage
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',[0.97 0.29 0.11]);
% 
%     set(gca, "Layer", "top")
%     hold(gca,'off')    
% 
%     % Limits
%     ylim([ymin1 ymax1]); 
%     xlim([tmin,tmax]);
% 
%     % Labels
%     xticks([tmin tmid tmax])
%     xticklabels({'','', ''})
%     yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])
%     yticklabels({'','',''})

% Row 4 --------------------------------------------------------------------

ax41 = axes('Units','centimeters', 'InnerPosition',[2 1.2 2.8 4]);
    set(ax41,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,1}(:,:); z = Stiff_D_Qs{2,1}; params = Stiff_D_params{2,1};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xlabel('$t$','Interpreter','latex','Rotation',0);
    ylabel('$Q$','Interpreter','latex','Rotation',0,'FontSize',12);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])    

    % Inset applied load
    ax410 = axes('Units','centimeters','InnerPosition',[2.1 1.4 0.75 0.75]);
        hold(gca,"on")
        box(gca,"on")
        set(gca,'FontName','Times','FontSize',8);

        % Shading (deloading)
        fill([tmid tmax tmax tmid], [-params.Astar*params.omega/2 -params.Astar*params.omega/2 params.Astar*params.omega/2 params.Astar*params.omega/2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])

        % Plot
        plot(x(index1:index3),params.adot(x(index1:index3)),'LineWidth',0.75,'Color','k');
        set(gca, "Layer", "top")
        hold(gca,'off')

        % Limits
        ylim([-params.Astar*params.omega/2 params.Astar*params.omega/2]); 
        xlim([tmin,tmax]);
    
        % Labels
        title('$\dot{a}^*$','Interpreter','latex')
        xticks([tmin tmid tmax])
        xticklabels({'','', ''})
        yticks([-params.Astar*params.omega/2 0 params.Astar*params.omega/2])
        yticklabels({'','',''})


ax42 = axes('Units','centimeters','InnerPosition',[5.2 1.2 2.8 4]);
    set(ax42,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,2}(:,:); z = Stiff_D_Qs{2,2}; params = Stiff_D_params{2,2};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xlabel('$t$','Interpreter','latex','Rotation',0);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax43 = axes('Units','centimeters','InnerPosition',[8.4 1.2 2.8 4]);
    set(ax43,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,3}(:,:); z = Stiff_D_Qs{2,3}; params = Stiff_D_params{2,3};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off')    

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xlabel('$t$','Interpreter','latex','Rotation',0);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

ax44 = axes('Units','centimeters','InnerPosition',[11.6 1.2 2.8 4]);
    set(ax44,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Unpack
    x = Stiff_D_Ts{2,4}(:,:); z = Stiff_D_Qs{2,4}(:,:); params = Stiff_D_params{2,4};

    % Find time indexes
    [~,index1] = min(abs(x-tmin));
    [~,index3] = min(abs(x-tmax));

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 10 equally spaced lines over one cycle between 0<Z<1
    plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
    plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
    plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
    plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
    plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
    plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
    plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
    plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
    plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
    plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
    plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
    
    % emphasize points around and at damage
    plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);

    set(gca, "Layer", "top")
    hold(gca,'off')    

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xlabel('$t$','Interpreter','latex','Rotation',0);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
    yticklabels({'','',''})

% Gradient - old version
% ax45 = axes('Units','centimeters','InnerPosition',[14.8 1.2 2.8 4]);
%     set(ax44,'DefaultAxesColorOrder',g11);
%     hold(gca,"on")
%     box(gca,"on")
%     set(gca,'FontName','Times','FontSize',12);
% 
%     % Unpack
%     x = Stiff_D_Ts{10}(:,:); z = Stiff_D_Qs{10}(:,:); params = Stiff_D_params{10};
% 
%     % Find time indexes
%     [~,index1] = min(abs(x-tmin));
%     [~,index3] = min(abs(x-tmax));
% 
%     % Shading (deloading)
%     fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
% 
%     % Plot 10 equally spaced lines over one cycle between 0<Z<1
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',b(1,:));
%     plot(x(index1:index3),z(index1:index3,40),'LineWidth',1,'Color',b(23,:)); % Z = 0.1
%     plot(x(index1:index3),z(index1:index3,80),'LineWidth',1,'Color',b(46,:)); % Z = 0.2
%     plot(x(index1:index3),z(index1:index3,120),'LineWidth',1,'Color',b(69,:)); % Z = 0.3
%     plot(x(index1:index3),z(index1:index3,160),'LineWidth',1,'Color',b(92,:)); % Z = 0.4
%     plot(x(index1:index3),z(index1:index3,200),'LineWidth',1,'Color',b(115,:)); % Z = 0.5
%     plot(x(index1:index3),z(index1:index3,240),'LineWidth',1,'Color',b(133,:)); % Z = 0.6
%     plot(x(index1:index3),z(index1:index3,280),'LineWidth',1,'Color',b(156,:)); % Z = 0.7
%     plot(x(index1:index3),z(index1:index3,320),'LineWidth',1,'Color',b(179,:)); % Z = 0.8
%     plot(x(index1:index3),z(index1:index3,360),'LineWidth',1,'Color',b(202,:)); % Z = 0.9
%     plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',b(end,:));
% 
%     % emphasize points around and at damage
%     plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',[0.97 0.29 0.11]);
% 
%     set(gca, "Layer", "top")
%     hold(gca,'off')    
% 
%     % Limits
%     ylim([ymin2 ymax2]); 
%     xlim([tmin,tmax]);
% 
%     % Labels
%     xlabel('$t$','Interpreter','latex','Rotation',0);
%     xticks([tmin tmid tmax])
%     xticklabels({'0','T/2', 'T'})
%     yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])
%     yticklabels({'','',''})

%% Export ----------------------------------------------------------------------

print('Stiffness-damage-profiles','-dpdf','-loose')
