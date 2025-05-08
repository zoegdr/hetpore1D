% plot diffusivity or advection x 1

Phi_0 = 0.55;
Phi = 0:0.01:2;

ds = 0.3;
ls = 0.5;

% Local
M = @(Z) 1-ds*exp(-(Z-ls).^2/(2*(1/16)^2));
dM = @(Z) 2*ds*(Z-ls)/2/(1/16)^2.*exp(-(Z-ls).^2/2/(1/16)^2);

% % Negative gradient
% M = @(Z) -ds*Z+1;
% dM = @(Z) -ds;
% 
% % Positive gradient
% M = @(Z) ds*Z-1;
% dM = @(Z) ds;

nu = 0.3;
L = @(Z) nu/(1-nu)*M(Z);
dL = @(Z) nu/(1-nu)*dM(Z);

k = (Phi/Phi_0+1).^3./(1+Phi);
D0 = k./(1+Phi);

J = 1 + Phi - Phi_0;
Z1 = 0.5; Z2 = 0.3;
D1 = D0.*(M(Z1)/2.*(1+1./(J.^2)) + L(Z1)/2.*(1-1./(J.^2)));
D2 = D0.*(M(Z2)/2.*(1+1./(J.^2)) + L(Z2)/2.*(1-1./(J.^2)));
A1 = dM(Z1)/2.*(J-1./J) + dL(Z1)/2.*(J+1./J-2);
A2 = dM(Z2)/2.*(J-1./J) + dL(Z2)/2.*(J+1./J-2);

%% Limits

ymin = 0.7;
ymax = 8;


%% Create figure ----------------------------------------------------------

figure

% Set figure total dimension
set(gcf,'Units','centimeters')
% Absolute print dimensions of figure. 
% [pos_from_left, pos_from_bottom, fig_width, fig_height]
set(gcf,'Position',[0 0 4 5])

set(gca,'Units','centimeters') 
% Relative positioning of the axes within the frame
% [ inset_from_left, inset_from_bottom, axes_width, axes_height]
set(gca, 'InnerPosition',[1.2 1.5 2.5 3])
box('on');
set(gca,'FontName','Times','FontSize',12);

xlim([0 2]);
% xticks([tmin tmid tmax])
% xticklabels({'0','T/2', 'T'})

% ylim([ymin,ymax]);
% yticks([1 2]);

hold(gca,'on')

% Plot
plot(Phi,A2,'LineWidth',1,'Color','k');
set(gca, "Layer", "top")
hold(gca,'off')

% Axes

ylabel('$D$','Interpreter','latex','Rotation',0)
xlabel('$\Phi$','Interpreter','latex')

%% Export ----------------------------------------------------------------------


% print('Full-Diffusion-Phi','-dpdf','-loose')
