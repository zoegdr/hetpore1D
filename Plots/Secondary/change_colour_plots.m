% Create figure
figure1 = figure('Units','normalized','Position',[0.5 0.5 0.8 0.65],'Name','Profile');

% Create axes
axes1 = axes;
hold(axes1,'on');

omega=10;

% find half periods
for n=1:9
    [~,index_k(n)] = min(abs(T_k-n*pi/omega));
    [~,index(n)] = min(abs(T-n*pi/omega));
    % [~,index_d(n)] = min(abs(T_d-n*pi/omega));
end

% make an array of ones and twos for each half period
p1_k = ones(size(T_k(1:index_k(1),1)));
p2_k = 2*ones(size(T_k(index_k(1)+1:index_k(2),1)));
p3_k = ones(size(T_k(index_k(2)+1:index_k(3),1)));
p4_k = 2*ones(size(T_k(index_k(3)+1:index_k(4),1)));
p5_k = ones(size(T_k(index_k(4)+1:index_k(5),1)));
p6_k = 2*ones(size(T_k(index_k(5)+1:index_k(6),1)));
p7_k = ones(size(T_k(index_k(6)+1:index_k(7),1)));
p8_k = 2*ones(size(T_k(index_k(7)+1:index_k(8),1)));
p9_k = ones(size(T_k(index_k(8)+1:index_k(9),1)));
p10_k = 2*ones(size(T_k(index_k(9)+1:end,1)));
ck = [p1_k;p2_k;p3_k;p4_k;p5_k;p6_k;p7_k;p8_k;p9_k;p10_k];

p1 = ones(size(T(1:index(1),1)));
p2 = 2*ones(size(T(index(1)+1:index(2),1)));
p3 = ones(size(T(index(2)+1:index(3),1)));
p4 = 2*ones(size(T(index(3)+1:index(4),1)));
p5 = ones(size(T(index(4)+1:index(5),1)));
p6 = 2*ones(size(T(index(5)+1:index(6),1)));
p7 = ones(size(T(index(6)+1:index(7),1)));
p8 = 2*ones(size(T(index(7)+1:index(8),1)));
p9 = ones(size(T(index(8)+1:index(9),1)));
p10 = 2*ones(size(T(index(9)+1:end,1)));
c = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10];

% p1_d = ones(size(T_d(1:index_d(1),1)));
% p2_d = 2*ones(size(T_d(index_d(1)+1:index_d(2),1)));
% p3_d = ones(size(T_d(index_d(2)+1:index_d(3),1)));
% p4_d = 2*ones(size(T_d(index_d(3)+1:index_d(4),1)));
% p5_d = ones(size(T_d(index_d(4)+1:index_d(5),1)));
% p6_d = 2*ones(size(T_d(index_d(5)+1:index_d(6),1)));
% p7_d = ones(size(T_d(index_d(6)+1:index_d(7),1)));
% p8_d = 2*ones(size(T_d(index_d(7)+1:index_d(8),1)));
% p9_d = ones(size(T_d(index_d(8)+1:index_d(9),1)));
% p10_d = 2*ones(size(T_d(index_d(9)+1:end,1)));
% cd = [p1_d;p2_d;p3_d;p4_d;p5_d;p6_d;p7_d;p8_d;p9_d;p10_d];

cc = [c,c];
cck = [ck ck];
% ccd = [cd,cd];

custom_colormap = [[0.9,0.2,0.4]; [1.0,0.8,0.0]];

xx1 = [T(:,1) T(:,1)];
xxk = [T_k(:,1) T_k(:,1)];
% xxd = [T_d(:,1) T_d(:,1)];

% yy1 = [S_star(1,:)' S_star(1,:)'];
% yy1 = [Permeability_k(:,50) Permeability_k(:,50)];

% yy1 = [Strain_k(:,50)];

% yy1 = [Ps(:,95),Ps(:,95)];
% yy2 = [Ps(:,5),Ps(:,5)];

% yy1 = [Ps_k(:,95),Ps_k(:,95)];
% yy2 = [Ps_k(:,5),Ps_k(:,5)];

% yy1 = [dP(:,100),dP(:,100)];

yy1 = [Flux_k(:,90),Flux_k(:,90)];
yy2 = [Flux_k(:,30),Flux_k(:,30)];

% yy1 = [Flux(:,100),Flux(:,100)];
% yyd = [Flux_d(:,50),Flux_d(:,50)];

% yy1 = [Wss(:,100),Wss(:,100)];
% yy2 = [Wss(:,30),Wss(:,30)];
yyk = [Wss_k(:,100),Wss_k(:,100)];

zz1 = zeros(size(xx1));
zzk = zeros(size(xxk));
% zzd = zeros(size(xxd));

hold on
hs1 = surf(xxk,yy1,zzk,cck,'EdgeColor','interp','FaceColor','none','LineWidth',3);
hs2 = surf(xxk,yy2,zzk,cck,'EdgeColor','interp','FaceColor','none','LineWidth',1,'LineStyle',':');
% hs1 = surf(xx1,yy1,zz1,cc,'EdgeColor','interp','FaceColor','none','LineWidth',1,'LineStyle',':');
% hsd = surf(xxd,yyd,zzd,ccd,'EdgeColor','interp','FaceColor','none','LineWidth',3);
% hsk = surf(xxk,yyk,zzk,cck,'EdgeColor','interp','FaceColor','none','LineWidth',3);

colormap(custom_colormap) ;     %// assign the colormap
shading flat                    %// so each line segment has a plain color
view(2) %// view(0,90)          %// set view in X-Y plane

% create x and y labels and render using latex
hy=ylabel('$Q$','interpreter','latex','fontsize',40,'Rotation',0);
hx=xlabel('$t$','interpreter','latex','fontsize',40,'fontname','Times New Roman');
% set the axis labels and numbers in the fontsize and style we want (gca means get current axis)
set(gca,'fontsize',40); set(gca,'fontname','Times New Roman');

% k(\frac{\hat{\Phi}_f}{J})
% \frac{\partial W}{\partial Z}=\hat{\Phi}_f