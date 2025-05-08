% Plot hysterisis loop

% In this code: 
% - Choose a variable to plot (Uss,Ss,Uf) assigned to
%   z. Applied load or displacement is assigned to x
% - Plot is restricted to one period, and plots every pi/4/omega (10 times)
% - Figure position and size are set according to units which can be found
%   in the code inventory
% - Colourmaps are loaded from custom_colourmap.m. It is assigned according
%   to the variable being plotted (eg. green for strain)
% - ylim depends on the variable being plotted and the parameters. Adjust
%   accordingly
% - xlim is [0 1] or [-1 0] depending on whether applied displacement or
%   load
% - Exports to current folder. Change name as appropriate
%
% Formatting (and other changes) when plotting
% - xlabel (a(t) or s^*(t)) on or off
% - ylabel (var) on or off
% - position etc
% - Title of export
% - ylim

%% Choose variable to plot -----------------------------------------------
% Uss, Ss, Uf (solid or fluid displacement, stress)

% var = 'Us';
var = 'S';
% var = 'Uf';

y = Ts;


%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(gcf,'Position',[0 0 4.5 5.5]) % 4 figures
% set(gcf,'Position',[0 0 6 7.2]) % 3 figures

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters') 
% Relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]
set(gca, 'InnerPosition',[1.45 1.1 2.9 4.13]) % 4 figures
% set(gca, 'InnerPosition',[1.5 1.5 4.2 5.5]) % 3 figures

box('on');
set(gca,'FontName','Times','FontSize',12);

%% Create colourmap --------------------------------------------------------

maps = custom_colourmap; % call custom colourmap function
g = maps.green;
p = maps.purple;
s = maps.red;
b = maps.blue;

%% Plot

if strcmp(var,'Us')
    x = params.S_star(Ts);
    z = Uss(:,1);
    
    ylabel('$U(0,t)/A^*_1$','Interpreter','latex','Rotation',90);
    ymin = 0; ymax = 0.8;
    ylim([ymin ymax]);
    yticks([0 0.4 0.8])
    xlabel('$s^*(t)/A^*_1$','Interpreter','latex','Rotation',0);
    xlim([0 1]);
    xticks([0 0.5 1]);
    
    hold on
    % Shading (deloading)
    fill([0.5 1 1 0.5], [ymin ymin ymax ymax], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    % Plot
    plot(x/params.Astar,z/params.Astar,'Color',p(150,:),'LineWidth',1)
    % Emphasise last cycle in dotted blue
    t1 = (params.p-1)*2*pi/params.omega;
    [~,index] = min(abs(y-t1));
    plot(x(index:end)/params.Astar,z(index:end)/params.Astar,'Color',[0.74,0.92,1],'LineStyle',':','LineWidth',2);


    set(gca, "Layer", "top")

elseif strcmp(var,'S')
    x = params.a(Ts);
    z = Ss(:,1);

    ylabel('$S(0,t)/A_2^*$','Interpreter','latex','Rotation',90);
    % xlabel('$a^*(t)/A_2^*$','Interpreter','latex','Rotation',0);
    ymin = -6; ymax = 3;
    ylim([ymin ymax])
    yticks([-6 -3 0  3])

    xlim([-1 0]);
    xticks([-1 -0.5 0]);
    
    hold on
    % Shading (deloading)
    fill([-0.5 0 0 -0.5], [ymin ymin ymax ymax], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    % Plot
    plot(x/params.Astar,z/params.Astar,'Color',[1.0,0.5,0.0],'LineWidth',1);
    % Emphasise last cycle in dotted black
    t1 = (params.p-1)*2*pi/params.omega;
    [~,index] = min(abs(y-t1));
    plot(x(index:end)/params.Astar,z(index:end)/params.Astar,'Color',"k",'LineStyle',':','LineWidth',2);


    set(gca, "Layer", "top")

end


%% Export

% print('U-hysteris-KC','-depsc','-loose')
print('S-hysteris-o4','-depsc','-loose')

