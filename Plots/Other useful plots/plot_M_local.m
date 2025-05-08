% Plot local stiffness/perm profile

%% Choose parameters
ls = 0.3;
ds = 0.35;

%% Plot

% Create figure
fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
% set(fig,'Position',[0 0 4 1.3])
% set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[4,1.3])
set(fig,'Position',[0 0 9 3.5])
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[9 3.5])

% Create axes
ax = axes('Units','centimeters','InnerPosition',[2 1 5 1.75]);
    hold(ax,"on")
    box(ax,"on")
    set(ax,'FontName','Times','FontSize',10);
    set(ax,"TickLabelInterpreter",'latex')


Zs = 0:0.01:1;
M = 1-ds*exp(-(Zs-ls).^2/(2*(1/16)^2));
plot(Zs,M,'LineWidth',1,'Color',"k");

box('off');

% Limits
ylim([0.3 1.05]);
xlim([0 1]);

% Labels

yticks([0.65 1])
yticklabels({'$1-d$' '1'})
% ylabel('$M$','Interpreter','latex','Rotation',0,'FontSize',8)

xticks([0 ls 1])
% xticklabels({'0' num2str(ls) '1'})
xticklabels({'0' '$l$' '1'})
xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',10)

title({'$f(Z)$'},'Interpreter','latex','FontSize',10)

%% Export
print('Gaussian','-dpdf','-r0')

