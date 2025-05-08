% 4x4 template

%% Set limits

% Frequency of applied load
omega = Local_params{1}.omega;
% Time limits
tmin = 18*pi/omega;
tmid = 19*pi/omega;
tmax = 20*pi/omega;

% Plot Limits
ymin1 = -0.3; ymax1 = 0.3; % Local decrease
ymin2 = -1; ymax2 = 1; % Local decrease


%% Create colour schemes  -------------------------------

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

% Top two rows
ax01 = axes('Units','centimeters','InnerPosition',[0.25 10.6 17.5 9.55]);
xticks([]); yticks([]);
box(ax01,"on"); ax01.XColor = [.82,.91,1]; ax01.YColor = [.82,.91,1]; ax01.LineWidth = 2;
text(9,9.75,"AL",Units="centimeters",FontSize=12,FontName='Times'); % annotation

% Bottom two rows
ax02 = axes('Units','centimeters','InnerPosition',[0.25 0.5 17.5 9.5]);
xticks([]); yticks([]);
box(ax02,"on"); ax02.XColor = [1,.89,.84]; ax02.YColor = [1,.89,.84]; ax02.LineWidth = 2;
text(9,9.75,"AD",Units="centimeters",FontSize=12,FontName='Times'); % annotation



%% Create subfigures

% Row 1 --------------------------------------------------------------------

ax11 = axes('Units','centimeters','InnerPosition',[2 15.9 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot
    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])

ax12 = axes('Units','centimeters','InnerPosition',[6 15.9 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax13 = axes('Units','centimeters','InnerPosition',[10 15.9 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot
    
    % emphasize points around and at damage
    % Local
    % plot(x(index1:index3),z(index1:index3,params.ls*params.N-40),'LineWidth',1,'Color',[1.0,0.6,0.5]);
    plot(x(index1:index3),z(index1:index3,params.ls*params.N),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    % plot(x(index1:index3),z(index1:index3,params.ls*params.N+40),'LineWidth',1,'Color',[0.7,0.1,0.0]);
    % Gradient
    % plot(x(index1:index3),z(index1:index3,end),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    % plot(x(index1:index3),z(index1:index3,1),'LineWidth',1,'Color',[0.7,0.1,0.0]);

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

ax14 = axes('Units','centimeters','InnerPosition',[14 15.9 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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



% Row 2 --------------------------------------------------------------------

ax21 = axes('Units','centimeters','InnerPosition',[2 11.4 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

    set(gca, "Layer", "top")
    hold(gca,'off')
    
    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])

ax22 = axes('Units','centimeters','InnerPosition',[6 11.4 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax23 = axes('Units','centimeters','InnerPosition',[10 11.4 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax24 = axes('Units','centimeters','InnerPosition',[14 11.4 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

% Row 3 --------------------------------------------------------------------

ax31 = axes('Units','centimeters','InnerPosition',[2 5.7 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

    set(gca, "Layer", "top")
    hold(gca,'off')

    % Limits
    ylim([ymin1 ymax1]); 
    xlim([tmin,tmax]);

    % Labels
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    yticks([ymin1 ymin1+(ymax1-ymin1)/2 ymax1])

ax32 = axes('Units','centimeters','InnerPosition',[6 5.7 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax33 = axes('Units','centimeters','InnerPosition',[10 5.7 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax34 = axes('Units','centimeters','InnerPosition',[14 5.7 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin1 ymin1 ymax1 ymax1], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

% Row 4 --------------------------------------------------------------------

ax41 = axes('Units','centimeters', 'InnerPosition',[2 1.2 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

    % Limits
    ylim([ymin2 ymax2]); 
    xlim([tmin,tmax]);

    % Labels
    xlabel('$t$','Interpreter','latex','Rotation',0);
    ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',14);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    yticks([ymin2 ymin2+(ymax2-ymin2)/2 ymax2])    

ax42 = axes('Units','centimeters','InnerPosition',[6 1.2 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);


    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot 

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

ax43 = axes('Units','centimeters','InnerPosition',[10 1.2 3 4]);
    set(0,'DefaultAxesColorOrder',b11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);

    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    
    % Plot

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

ax44 = axes('Units','centimeters','InnerPosition',[14 1.2 3 4]);
    set(0,'DefaultAxesColorOrder',g11);
    hold(gca,"on")
    box(gca,"on")
    set(gca,'FontName','Times','FontSize',12);
    
    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin2 ymin2 ymax2 ymax2], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])

    % Plot
    
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

