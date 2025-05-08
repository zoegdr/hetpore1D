custom_colormap = [[0.9,0.2,0.4]; [1.0,0.8,0.0]];

omega=10;

% find half periods
for n=1:9
    [~,index_k(n)] = min(abs(T_k-n*pi/omega));
    [~,index(n)] = min(abs(T-n*pi/omega));
    % [~,index_d(n)] = min(abs(T_d-n*pi/omega));
end

% make an array of ones and twos for each half period
p1 = ones(size(T_k(1:index_k(1),1)));
p2 = 2*ones(size(T_k(index_k(1)+1:index_k(2),1)));
p3 = ones(size(T_k(index_k(2)+1:index_k(3),1)));
p4 = 2*ones(size(T_k(index_k(3)+1:index_k(4),1)));
p5 = ones(size(T_k(index_k(4)+1:index_k(5),1)));
p6 = 2*ones(size(T_k(index_k(5)+1:index_k(6),1)));
p7 = ones(size(T_k(index_k(6)+1:index_k(7),1)));
p8 = 2*ones(size(T_k(index_k(7)+1:index_k(8),1)));
p9 = ones(size(T_k(index_k(8)+1:index_k(9),1)));
p10 = 2*ones(size(T_k(index_k(9)+1:end,1)));
ck = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10];

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

cc = repmat(c,1,100);
cck = repmat(ck,1,100);
% ccd = repmat(cd,1,100);

% hs = waterfall(T,Zs,Ps',cc');
% hsk = waterfall(T_k,Zs_k,Ps_k',cck');
% hs = waterfall(T,Zs,Stress',cc');
hsk = waterfall(T_k,Zs_k,Stress_k',cck');
% hsk = waterfall(T,Zs,Phis',cc');
% hsk = waterfall(Zs,T,Phis,cc);
% hs1 = waterfall(T_d,Zs_d,Flux_d',ccd');
% hs1 = waterfall(T_d,Zs_d,Flux_d',ccd');

% hs = waterfall(T,Zs,Ps',cc');

ylabel('$Z$','interpreter','latex','fontsize',25,'Rotation',0);
xlabel('$t$','interpreter','latex','fontsize',25,'fontname','T_kimes New Roman');
zlabel('$Q$','interpreter','latex','fontsize',25,'fontname','T_kimes New Roman','Rotation',0);
set(gca,'fontsize',25); set(gca,'fontname','T_kimes New Roman','xdir','reverse');

colormap(custom_colormap) ;     %// assign the colormap
shading flat                    %// so each line segment has a plain color