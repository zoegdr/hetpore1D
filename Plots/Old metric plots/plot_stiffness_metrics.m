% Heatmaps for metrics varying damage, location and omega

% In this code: 
% - Choose a metric to plot. For local damage, choose omega to plot
%   (index=1,2,3 for omega=1,10,30)
% - Figure position and size are set according to units which can be found
%   in the code inventory
% - Colourmaps are loaded from custom_colourmap.m. It is assigned according
%   to the variable being plotted (eg. green for strain)
% - clim depends on the variable being plotted and the parameters. Adjust
%   accordingly
% - Exports to current folder. Change name as appropriate
%
% Formatting (and other changes) when plotting
% - xlabel (location or omega) on or off
% - xtick labels
% - ylabel (ds) on or off
% - colorbar on or off
% - inner position depending on if left, middle or right plot
% - Title of export
% - clim

%% Choose a metric and omega(index) to plot --------------------------------------

index = 1; % omegas = [1 10 30]

% var = 'max strain';
var = 'max nominal strain';
% var = 'max flow into damage';
% var = 'min porosity';
% var = 'max strain gradient';

%% Create figure ------------------------------------------------------------------

figure

% Set figure total dimension
% set(gcf, 'PaperPositionMode', 'manual','Renderer','painters');
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 6 4.5]) % Absolute print dimensions of figure. 
% set(gcf,'PaperSize',[6 4.5])
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
% 
% % Set size and position of axes plotting area within figure dimensions. To
% % keep vertical axes aligned for multiple figure keep the horizontal
% % position consistent
% set(gca,'Units','normalized') % relative to parent figure
% set(gca, 'InnerPosition',[0.2 0.2 0.65 0.75])% This is the relative positioning of the axes within the frame
% % [ inset_from_left, inset_from_bottom, axes_width, axes_height]

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters')
set(gca, 'InnerPosition',[1.2 0.9 3.9 3.375])% This is the relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]
% [1.2 .... ] left figure (index=1)
% [1.05 .... ] middle figure (index=10)
% [0.9 .... ] right figure (index=30)


%% Create colourmaps

% Colourmaps
maps = custom_colourmap;
p = maps.purple;
g = maps.green;
r = maps.red;
b = maps.blue;
t = maps.turquoise;

locations = [0.2 0.5 0.8];
magnitudes = [0 0.2 0.5 0.8];
omegas = [1 10 30];

%% Plot ---------------------------------------------------------------------------

if strcmp(var,'max strain')
    h = heatmap(locations,magnitudes,max_strain_at_damage(:,:,index),Colormap=g);
    
    % Max strain
    clim([-1.2 0]);
    hs = struct(h);
    hs.Colorbar.Ticks = [0 0.6 1.2];

elseif strcmp(var,'max nominal strain')
    h = heatmap(locations,magnitudes,max_nominal_strain(:,:,index),Colormap=g);
    
    % Nominal Strain
    clim([0.1 0.3]);
    hs = struct(h);
    hs.Colorbar.Ticks = [0.1 0.2 0.3];

elseif strcmp(var,'max flow into damage')
    h = heatmap(locations,magnitudes,max_flow_into_damage(:,:,index),Colormap=b);
    
    % Max flow into damage
    % clim([0 3]);
    % hs = struct(h);
    % hs.Colorbar.Ticks = [0 1 2 3];

elseif strcmp(var,'min porosity')
    h = heatmap(locations,magnitudes,min_total_porosity(:,:,index),Colormap=t);

    % Min porosity
    clim([0.55 0.6])
    hs = struct(h);
    hs.Colorbar.Ticks = [0.55 0.57 0.6];

elseif strcmp(var,'max strain gradient')
    h = heatmap(omegas,-magnitudes,max_strain_positive_gradient_disp,Colormap=flip(g));

    clim([-1.2 0]);
    hs = struct(h);
    hs.Colorbar.Ticks = [0 0.6 1.2];

end


%% Properties  ---------------------------------------------------------

h.XLabel = '$l_s$';
h.YLabel = '$d_s$';
% h.XLabel = '$\omega$';

h.XDisplayLabels = {'0.2','0.5','0.8'};
% h.XDisplayLabels = {'1','10','30'};
h.YDisplayLabels = {'0','0.2','0.5','0.8'};

colorbar('off');

sorty(h,"0.2","descend");
% sorty(h,"1","ascend");

set(gca,'FontName','Times','FontSize',9,'interpreter','latex');
h.CellLabelFormat = '%.2f';


%% Export ---------------------------------------------------------------
% print('min-porosity-metric-o30','-depsc','-loose')
% print('max-strain-metric-o30','-depsc','-loose')
% print('max-strain-metric-pos-grad-disp','-depsc','-loose')
print('max-nominal-strain-metric-o1','-depsc','-loose')
% print('max-flow-into-dam-metric-o30','-depsc','-loose')
