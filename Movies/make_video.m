% Open video writer
v = VideoWriter('Perm-damage-0,8-flux','MPEG-4');
open(v);

% Assign variables
t = 0:Perm_D_Ts{1,1}(end)/500:Perm_D_Ts{1,1}(end);

% x = Zs; y0 = interp1(Perm_D_Ts{1,1},Perm_D_dUdZ{1,1},t); y1 = interp1(Perm_D_Ts{1,2},Perm_D_dUdZ{1,2},t); % strain
x = Zs; y0 = interp1(Perm_D_Ts{1,1},Perm_D_Qs{1,1},t); y1 = interp1(Perm_D_Ts{1,2},Perm_D_Qs{1,2},t); % flux


% Colour scheme
c1 = [0.4693    0.7545    0.5077]; % green
c2 = [0.97 0.29 0.11]; % red
c3 = "#1B97E6"; % blue
c4 = "#FFCE3C"; % yellow

% Create figure and set dimensions
h = figure;
% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 20 26])
ax = axes('Units','centimeters','InnerPosition',[3 3 16 22]);

% Run through time
for k = 1:length(t)

    % Plot
    hold(ax,'off')
    plot(x,y0(k,:),'LineWidth',1,'Color', "k");
    hold(ax,'on')
    plot(x,y1(k,:),'LineWidth',1,'Color',c2);
    % plot(x,y2(k,:),'LineWidth',1,'Color',c2);
    % plot(x,y3(k,:),'LineWidth',1,'Color',c3);
    % plot(x,y4(k,:),'LineWidth',1,'Color',c4);

    % Formatting
    box(ax,"on")
    set(ax,'FontName','Times','FontSize',25);

    % Limits
    % ylim([0 0.3]) % strain
    % yticks([0 0.15 0.3])
    ylim([-0.5 0.5]) % flux
    yticks([-0.5 0 0.5])
    xlim([0 1])
    xticks([0 0.2 0.4 0.6 0.8 1])

    % Labels
    % ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',30);
    xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',30)

    % Legend
    % legend('Neo-Hookean','UC SEF $M=10^{-4}$','UC SEF $M=10^{-2}$','UC SEF $M=10^{-1}$','UC SEF $M=1$','Interpreter','latex','Location','NorthEast','NumColumns',1,'Box','off');
    legend('No damage','Damage','Interpreter','latex','Location','NorthEast','NumColumns',2,'Box','off');


    % Save frame
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);