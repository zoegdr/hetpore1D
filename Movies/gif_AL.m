%% Set parameters and variables

Zs = 0+1/2/400:1/400:1-1/2/400;


% Unpack
x = Zs; y1 = Phi0; y2 = interp1(Ts,Phis,T0);
c = [0.4693    0.7545    0.5077]; % green

% x = Zs; y1 = interp1(Stiff_D_Ts{6},Stiff_D_dUdZ{6}(:,:),Stiff_D_Ts{1}); y2 = interp1(Stiff_D_Ts{8},Stiff_D_dUdZ{8}(:,:),Stiff_D_Ts{1});
% x = Zs; y1 = interp1(T2{5,1},Phi2{5,1},T0); y2 = interp1(T2{5,4},Phi2{5,4}(:,:),T0);
% x = Zs; y1 = interp1(Perm_D_Ts{2,1},Perm_D_dUdZ{2,1},Perm_D_Ts{1,1}); y2 = interp1(Perm_D_Ts{2,3},Perm_D_dUdZ{2,3},Perm_D_Ts{1,1});
% c = [0.4462    0.7229    0.8994]; % blue


h = figure;

% Set figure total dimension
set(h,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(h,'Position',[0 0 10 13])

ax = axes('Units','centimeters','InnerPosition',[2 2 7 10]);
hold(ax,"on")
plot(x,y1(1,:),'LineWidth',1,'Color',[0.4693    0.7545    0.5077]);
plot(x,y2(1,:),'LineWidth',1,'Color',[0.97 0.29 0.11]);
box(ax,"on")
set(ax,'FontName','Times','FontSize',15);

% Limits
ylim([0 0.15]); 
xlim([0 1]);

% Labels
ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',15);
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.05 0.1 0.15])

gifFile = 'testing.gif';
exportgraphics(h, gifFile);

for i = 2:length(y1)
    
    hold(ax,"off")
    ax.NextPlot = 'replaceChildren';
    plot(x,y1(i,:),'LineWidth',1,'Color', c);
    hold(ax,"on")
    plot(x,y2(i,:),'LineWidth',1,'Color',[0.97 0.29 0.11]);
    
    exportgraphics(h,gifFile,Append=true);

end