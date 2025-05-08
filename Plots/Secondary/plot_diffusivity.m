% plot diffusivity x 3

omega = params.omega;

D_c = 1./(1+Phis);
D_1 = 1*ones(size(Phis));
D_KC = (Phis_KC/params_KC.Phi0+1).^3./((1+Phis_KC).^2);

%% Limits

tmin = 18*pi/omega;
tmid = 19*pi/omega;
tmax = 20*pi/omega;

ymin = 0.7;
ymax = 2.5;


%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(gcf,'Position',[0 0 5 8.75])

ax1 = subplot(3,1,1);
    % Set size and position of axes plotting area within figure dimensions. To
    % keep vertical axes aligned for multiple figure keep the horizontal
    % position consistent
    set(ax1,'Units','centimeters') 
    % Relative positioning of the axes within the frame
    % [ inset_from_left, inset_from_bottom, axes_width, axes_height]
    set(ax1, 'InnerPosition',[0.8 6 4 2])
    box('on');
    set(gca,'FontName','Times','FontSize',12);

    xlim([tmin,tmax]);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    
    ylim([ymin,ymax]);
    yticks([1 2]);

    % Find indexes
    [~,index1] = min(abs(Ts_KC-tmin));
    [~,index3] = min(abs(Ts_KC-tmax));
    
    hold(ax1,"on")
    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin ymin ymax ymax], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    % Plot
    plot(Ts_KC(index1:index3),D_KC(index1:index3,1),'LineWidth',1,'Color','k');
    set(ax1, "Layer", "top")
    hold(ax1,'off')

    title('$D(Z=0)$','Interpreter','latex','Rotation',0)

ax2 = subplot(3,1,2);
    % Set size and position of axes plotting area within figure dimensions. To
    % keep vertical axes aligned for multiple figure keep the horizontal
    % position consistent
    set(ax2,'Units','centimeters') 
    % Relative positioning of the axes within the frame
    % [ inset_from_left, inset_from_bottom, axes_width, axes_height]
    set(ax2, 'InnerPosition',[0.8 3.75 4 2])
    box('on');
    set(ax2,'FontName','Times','FontSize',12);

    xlim([tmin,tmax]);
    xticks([tmin tmid tmax])
    xticklabels({'','', ''})
    
    ylim([ymin,ymax]);
    yticks([1 2]);

    % Find indexes
    [~,index1] = min(abs(Ts-tmin));
    [~,index3] = min(abs(Ts-tmax));
    
    hold(ax2,'on')
    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin ymin ymax ymax], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    % Plot
    plot(Ts(index1:index3),D_c(index1:index3,1),'LineWidth',1,'Color','k');
    set(ax2, "Layer", "top")
    hold(ax2,'off')

 ax3 = subplot(3,1,3);
    % Set size and position of axes plotting area within figure dimensions. To
    % keep vertical axes aligned for multiple figure keep the horizontal
    % position consistent
    set(ax3,'Units','centimeters') 
    % Relative positioning of the axes within the frame
    % [ inset_from_left, inset_from_bottom, axes_width, axes_height]
    set(ax3, 'InnerPosition',[0.8 1.5 4 2])
    box('on');
    set(ax3,'FontName','Times','FontSize',12);

    xlim([tmin,tmax]);
    xticks([tmin tmid tmax])
    xticklabels({'0','T/2', 'T'})
    
    ylim([ymin,ymax]);
    yticks([1 2]);

    % Find indexes
    [~,index1] = min(abs(t-tmin));
    [~,index3] = min(abs(t-tmax));
    
    hold(ax3,'on')
    % Shading (deloading)
    fill([tmid tmax tmax tmid], [ymin ymin ymax ymax], [0.95,0.95,0.95],EdgeColor=[0.95,0.95,0.95])
    % Plot
    plot(t(index1:index3),D_1(index1:index3,1),'LineWidth',1,'Color','k');
    set(ax3, "Layer", "top")
    hold(ax3,'off')


%% Export ----------------------------------------------------------------------


print('Diffusion','-depsc','-loose')
