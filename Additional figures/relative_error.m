% Relative Error

%% Set parameters

params1 = default_parameters('stress');
params2 = default_parameters('disp');
params1.p = 5;
params2.p = 5;

params1.Astar = 0.2;
params2.Astar = 0.1;

%% Calculate errors

Ns = [50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000];
error_AS = zeros(1,length(Ns));
error_AD = zeros(1,length(Ns));

% Initiate
params.N = 25;
[params1,~,~,~,Uss,~,~,~,~,~] = cyclic_uniaxial_stress(params1);
Uss_new_AS = Uss(end,1);
[params2,~,~,~,Uss,~,~,~,~,~] = cyclic_uniaxial_disp(params2);
Uss_new_AD = Uss(end,1);

% Loop through Ns
for i=1:length(Ns)

    Uss_old_AS = Uss_new_AS; % Previous new value becomes old value
    Uss_old_AD = Uss_new_AD; % Previous new value becomes old value

    params1.N = Ns(i); % Run for new N
    params2.N = Ns(i); % Run for new N

    [params1,~,~,~,Uss,~,~,~,~,~] = cyclic_uniaxial_stress(params1);
    Uss_new_AS = Uss(end,1);
    error_AS(i) = abs((Uss_new_AS-Uss_old_AS)/Uss_new_AS); % Relative error AS

    [params2,~,~,~,Uss,~,~,~,~,~] = cyclic_uniaxial_disp(params2);
    Uss_new_AD = Uss(end,1);
    error_AD(i) = abs((Uss_new_AD-Uss_old_AD)/Uss_new_AD); % Relative error AD
end

%% Make plot

fig = figure;

% Set figure total dimension
set(fig,'Units','centimeters')
set(fig,'Position',[0 0 8.2 8.2]) % Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[8.2,8.2])

% Set size and position of axes plotting area within figure dimensions. To
% keep vertical axes aligned for multiple figure keep the horizontal
% position consistent
set(gca,'Units','centimeters') 
set(gca, 'InnerPosition',[1.7 1.3 6 6])% This is the relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]

box('on');
hold("on")

plot(log(Ns),log(error_AS),'Color',"k",'LineWidth',1)
plot(log(Ns),log(error_AD),'Color',"k",'LineWidth',1,'LineStyle',':')

set(gca,'FontName','Times','FontSize',12);
xlabel('$\log(N)$','Interpreter','latex','Rotation',0);
ylabel('$\log(\varepsilon_{rel})$','Interpreter','latex');
ylim([-15 5])
xticks([4 5 6 7]);
% yticks([0 0.005 0.01 0.015])

legend('AS','AD','Location','best','NumColumns',1,'Box','off','IconColumnWidth',15);

% print('relative-error','-dpdf','-r0')