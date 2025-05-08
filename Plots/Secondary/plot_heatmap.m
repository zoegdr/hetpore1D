% Plot heatmaps of solution

% In this code: 
% - Choose a variable to plot (Phis,dUdZ,Ss,P,Q...) assigned to
%   z and var. Time is assigned to x and space is assigned to y.
% - Choose a time to plot (all time 'all', initial response 'init', or few periods 'period')

function plot_heatmap(x,y,z,params,var,time)

% Frequency of applied load
omega = params.omega;
Astar = params.Astar;
p = params.p;

%% Choose time period to plot ---------------------------------------------

if strcmp(time,'all')
    tmin = 0;
    tmax = p*2*pi/omega; 
elseif strcmp(time,'init')
    tmin = 0;
    tmax = 6*pi/omega;
elseif strcmp(time,'period')
    tmin = 10*pi/omega;
    tmax = 12*pi/omega;
else
    tmin = 6*pi/omega;
    tmax = 12*pi/omega; %p*pi/omega;
end

% Find index
[~,index1] = min(abs(x-tmin));
[~,index2] = min(abs(x-tmax));

%% Create figure ----------------------------------------------------------

figure1 = figure;
figure1.Position = [0 0 700 400];

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
pcolor(x(index1:index2),y,z(index1:index2,:)');

set(gca,'YDir','reverse'); % if applying boundary condition at Z=0 to get
%same configuration of plots but with reverse Z axis

shading flat;
if strcmp(var,'A')
    colormap("cool");
elseif strcmp(var,'D')
    colormap("spring");
else
    colormap("jet");
end

%% Limits & Title ---------------------------------------------------------

%%% k=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clim([0,1]); %D
% clim([0,7.5]); %strain 
% clim([-4,3]); %pressure 
% clim([-4.5,6.5]); %flux 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% k(Phi) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(params.damage,'local-') || strcmp(params.damage,'local+')
    if strcmp(var,'D')
        clim([0.5,1.6]); %D
        title('D');
    elseif strcmp(var,'A')
        clim([-1.5,1.5]); %A
        title('A');
    elseif strcmp(var,'dUdZ')
        clim([0,2*Astar]); % strain appears to be max 2*Astar for ds=0.6
        title('Strain');
    elseif strcmp(var,'U')
        clim([0,0.08]); %disp
        title('U');
    elseif strcmp(var,'P')
        clim([-2,1.5]); %pressure
        title('P');
    elseif strcmp(var,'Q')
        clim([-5,5]); %flux
        title('Q');
    else 
        % no limit
    end

elseif strcmp(params.damage,'dec2') || strcmp(params.damage,'dec1')
    if strcmp(var,'D')
        clim([0.5,1.6]); %D
        title('D');
    elseif strcmp(var,'A')
        clim([-0.15,0.15]); %A
        title('A');
    elseif strcmp(var,'dUdZ')
        clim([0,2*Astar]); %strain
        title('Strain');
    elseif strcmp(var,'U')
        clim([0,0.08]); %disp
        title('U');
    elseif strcmp(var,'P')
        clim([-2,1.5]); %pressure
        title('P');
    elseif strcmp(var,'Q')
        clim([-2,2]); %flux
        title('Q');
    else 
        % no limit
    end

else 
    if strcmp(var,'D')
        clim([0.5,1.6]); %D
        title('D');
    elseif strcmp(var,'A')
        clim([-0.15,0.15]); %A
        title('A');
    elseif strcmp(var,'dUdZ')
        clim([0,2*Astar]); %strain
        title('Strain');
    elseif strcmp(var,'U')
        clim([0,Astar]); %disp
        title('U');
    elseif strcmp(var,'P')
        clim([-0.8,0.8]); %pressure
        title('P');
    elseif strcmp(var,'Q')
        clim([-3,3]); %flux
        title('Q');
    else 
        % no limit
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time limit 
xlim([x(index1),x(index2)]);


% space (Z) limit
ylim([0,1]);

%% Labels -----------------------------------------------------------------

colorbar

% Create xlabel
xlabel('$t$','Interpreter','latex','Rotation',0,'FontSize',40);

% Create ylabel
ylabel('$Z$','Interpreter','latex','FontSize',40,'Rotation',0);
% ylabel('$s''_{ZZ}$','Interpreter','latex','FontSize',36,'Rotation',0);

box(axes1,'on');
hold(axes1,'off');
set(axes1,'FontName','Times','FontSize',40);

end