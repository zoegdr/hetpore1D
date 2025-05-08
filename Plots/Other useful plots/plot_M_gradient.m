% Plot gradient stiffness profile

%% Default settings
set(0,'DefaultAxesFontSize',12);
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultFigureColor',"w");
set(0,'DefaultAxesFontName','times');

%% Variables

Zs = 0:0.01:1;
ls = 0.7;
ds = 0.3;

M_pos = ds*Zs+1-ds; % Positive
M_neg = -ds*Zs+1; % Negative

%% Create figure

figure

% Set figure total dimension
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 9 3])

%% Axes 

ax1 = subplot(1,2,1);

    % Set size and position of axes plotting area within figure dimensions. To
    % keep vertical axes aligned for multiple figure keep the horizontal
    % position consistent
    % This is the relative positioning of the axes within the frame
    % [ inset_from_left, inset_from_bottom, axes_width, axes_height]
    set(ax1,'Units','centimeters')
    set(ax1, 'InnerPosition',[1 1 3.6 1.8])

    plot(Zs,M_neg,'LineWidth',1,'Color',"k");
    box(ax1,'off');

    % Limits
    ylim([0.4 1]);
    xlim([0 1]);
    
    % Labels
    
    set(gca,'TickLabelInterpreter','latex')
    
    yticks([0.7 1])
    yticklabels({'$1-d_s$' '1'})
    ylabel('$M$','Rotation',0,'Position',[-0.2 0.85])
    
    xticks([0 1])
    xlabel('$Z$')

ax2 = subplot(1,2,2);
    
    % Set size and position of axes plotting area within figure dimensions. To
    % keep vertical axes aligned for multiple figure keep the horizontal
    % position consistent
    % This is the relative positioning of the axes within the frame
    % [ inset_from_left, inset_from_bottom, axes_width, axes_height]
    set(ax2,'Units','centimeters')
    set(ax2, 'InnerPosition',[5 1 3.6 1.8])

    plot(Zs,M_pos,'LineWidth',1,'Color',"k");
    box(ax2,'off');

    % Limits
    ylim([0.4 1]);
    xlim([0 1]);
    
    % Labels
    
    set(gca,'TickLabelInterpreter','latex')
    
    yticks([0.7 1])
    yticklabels({'' ''})
    
    xticks([0 1])
    xlabel('$Z$')

% % Export
print('Gradient-decrease-stiffness-profiles','-depsc','-loose')

