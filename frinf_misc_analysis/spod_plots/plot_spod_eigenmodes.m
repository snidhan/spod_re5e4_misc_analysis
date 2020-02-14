        %% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures_2.0/';

%% Generating the figure object

close all;
x0=0;
y0=0;
width=15;
height=15;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(3,2,[.02, 0.05],[.1,.05],[.1 .05]);
%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Lk

% axes(ha(1));

filename = './files/eigenmodes_similarity_diff_loc.mat';
load(filename);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 30$';
Legend{3} = '$x/D = 40$';
Legend{4} = '$x/D = 50$';
Legend{5} = '$x/D = 60$';
Legend{6} = '$x/D = 70$';
Legend{7} = '$x/D = 80$';
Legend{8} = '$x/D = 90$';
Legend{9} = '$x/D = 100$';

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ylim([0 1.2]);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
% set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
% hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


hLegend = legend(Legend, 'Numcolumns', 9);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% Inset plot

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(u_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

% ax2 = axes('Position',[ha(1).Position(1)+0.14 ha(1).Position(2)+0.1 .12 .12]); % For paper

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ux mode as a function of downstream distance using Ld

% clearvars -except dirout;
% %close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 6]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 utheta mode as a function of downstream distance using Lk

% axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20;
xlim([0 6]);
xticks([0 2 4 6]);
% set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','i`nterpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
% hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


% Inset plot

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(v_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

% ax2 = axes('Position',[ha(3).Position(1)+0.2 ha(3).Position(2)+0.075 .18 .18]);

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);


[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');
set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_v_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Ld

% clearvars -except dirout;
% % close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% 
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 6]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');


%% Similarity of m = 2, St = 0 ux mode as a function of downstream distance using Lk

% axes(ha(5));
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
xticklabels({'$0$', '$2$', '$4$', '$6$'});
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

% hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',20);
hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [-0.1, 0.45, 0]);
% hTitle  = title('(e)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(w_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

% ax2 = axes('Position',[ha(5).Position(1)+0.2 ha(5).Position(2)+0.075 .18 .18]);

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.27 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Ld

% clearvars -except dirout;
% % close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 6]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Lk

% clearvars -except dirout;
% close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 6]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ur mode as a function of downstream distance using Ld

% axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,6,2,i)), 'Color',  lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
% set(gca, 'Xticklabel', []);
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
% set(gca, 'Yticklabel', []);

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);
% hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(u_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

% ax2 = axes('Position',[ha(2).Position(1)+0.2 ha(2).Position(2)+0.075 .18 .18]);

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;


% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ur mode as a function of downstream distance using Lk

% clearvars -except dirout;
% % close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 6]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 utheta mode as a function of downstream distance using Ld

% axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';



hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);
% hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(v_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

% ax2 = axes('Position',[ha(4).Position(1)+0.2 ha(4).Position(2)+0.075 .18 .18]);

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;



% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'similarity_v_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Lk

% clearvars -except dirout;
% % close all;
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% % C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% % color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% 
% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% 
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% ax = gca;
% ax.FontSize = 20;
% 
% 
% 
% % hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% % hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Ld

% axes(ha(6));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end


ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
% xticklabels({'$0$', '$2$', '$4$', '$6$'});
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
% set(gca, 'Yticklabel', []);
box on;

% hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',20);
hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hTitle  = title('(f)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(w_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

% ax2 = axes('Position',[ha(4).Position(1)+0.2 ha(4).Position(2)+0.075 .18 .18]);

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
yticklabels({'$-10$', '$0$', '$10$'});


xlim([-13 13]);
xticks([-10 0 10]);
xticklabels({'$-10$', '$0$', '$10$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_m1st0136_m2st0_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_m1st0136_m2st0_x_D', '.eps'),'-depsc2','-r600');
%% Three dimensional slices of real part of ur m=1 St=0.135

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));
% 
% 
% % Fix the sign of mode based on the x/D = 10 mode
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ux_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         disp(i);
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, real(ur_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta)); %#ok<*SAGROW>
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ur_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1]);
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ur_real_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ur_real_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ux m=1 St=0.135

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, real(ux_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ux_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1]);
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.eps'),'-depsc2','-r600');
% 
%% Three dimensional slices of imag part of ur m=1 St=0.135

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));
% 
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, imag(ur_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ur_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1]);
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ur_imag_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ur_imag_m1st0135_x_D_slices.eps'),'-depsc2','-r600');
% 
%% Three dimensional slices of imag part of ux m=1 St=0.135

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, imag(ux_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ux_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1]);
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_imag_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_imag_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ur m=2 St=0

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));
% 
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, real(ur_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta)); %#ok<*SAGROW>
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ur_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1]);
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ur_real_m2st0_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ur_real_m2st0_x_D_slices.eps'),'-depsc2','-r600');
% 
%% Three dimensional slices of real part of ux m=2 St=0

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, real(ux_mode(:,i)), '-', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ux_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1])
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_real_m2st0_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_real_m2st0_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ur m=2 St=0

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));
% 
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, imag(ur_mode(:,i)), '--', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ur_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1])
% 
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ur_imag_m2st0_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ur_imag_m2st0_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ux m=2 St=0

% clearvars -except dirout;
% close all;
% 
% filename = './files/eigenmodes_similarity_diff_loc.mat';
% 
% load(filename);
% 
% theta = linspace(0,2*pi,256);
% 
% ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
% ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));
% 
% % Fix the sign of mode based on the x/D = 10 mode
% 
% real_comp_mode_x10 = real(ux_mode(:,2));
% 
% for i = 1:size(ur_mode,2)
%     real_comp_mode = real(ux_mode(:,i));
%     if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
%         ux_mode(:,i) = -1*ux_mode(:,i);
%         ur_mode(:,i) = -1*ur_mode(:,i);
%     end
% end
% 
% figure; hold on;
% for i = 2:3:20
%     plot(rc, imag(ux_mode(:,i)), '-', 'Linewidth', 2);
% end
% close;
% 
% for i = 1:size(ur_mode,2)
%     real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
%     imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
%     imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));
% 
% end
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure; 
% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
% hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
% view(3);
% 
% cameratoolbar('SetCoordSys','x');
% 
% box on;
% hold on;
% 
% ax = gca;
% ax.FontSize = 16;
% 
% C = cell(7,1);
% h = cell(7,1);
% 
% count = 1;
% for ct = 2:3:20
%     [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ux_mode_contour(:,:,ct)),10);
%     axis equal;
%     ax = gca;
%     ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z
% 
% end
% 
% colormap('hot');
% 
% zlim([0 120]);
% ylim([-15 15]);
% xlim([-15 15]);
% 
% zticks([0 20 40 60 80 100]);
% yticks([-10 0 10]);
% xticks([-10 0 10]);
% camorbit(180,0,'camera');
% 
% daspect([1 1 1])
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_imag_m2st0_x_D_slices.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_imag_m2st0_x_D_slices.eps'),'-depsc2','-r600');
% 


%% Plotting the snapshots representative of m=1, 2 mode at x/D = 10;

clearvars -except dirout;

% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;

% Read w velocity field 

nr = 356; ntheta = 258;
var3 = 'wp';
x = 10;
dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
w = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            w(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
w_centered = w(2:nr-1,2:ntheta-1,:);
w_mean = squeeze(mean(w_centered,3));
w_fluc = w_centered - w_mean;
w_fluc_sampled = w_fluc(:,:,1:end);
clear w_fluc; clear w_centered; clear w;

% Reading the grid files

ntheta = 256;
theta = linspace(0,2*pi,ntheta)';
numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Generating the plot object
close all;
x0=0;
y0=0;
width=15;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,3,[.05, 0.05],[.1,.15],[.1 .05]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Plotting the m=2 mode at x/D = 10

axes(ha(1));

load('./files/coefficients_projection_x_D_10.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
%h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
%colormap hot;
%colorbar;

t = 4663;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>

set(h,'Linecolor','none');
colormap hot;
caxis([-0.3 0.3]);


hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(1).Position(1) ha(1).Position(2)+ 0.36 ha(1).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.3:0.15:0.3);
set(hBar, 'TickLabelInterpreter', 'latex');

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);

ax = gca;
ax.FontSize = 20;

xlim([-3 3]);
xticks([-3 -2 -1 0 1 2 3]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-3 3]);
yticks([-3 -2 -1 0 1 2 3]);
yticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});

%print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_10.png'),'-dpng2','-r600');  
%print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_10.eps'),'-depsc2','-r600');
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

box on;

%% Plotting the m=1 mode at x/D = 10

axes(ha(4));

load('./files/coefficients_projection_x_D_10.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;z
% colorbar;


t = 2641;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);

set(h,'Linecolor','none');
colormap hot;
caxis([-0.3 0.3]);
%colorbar;



ax = gca;
ax.FontSize = 20;


xlim([-3 3]);
xticks([-3 -2 -1 0 1 2 3]);
xticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});
ylim([-3 3]);
yticks([-3 -2 -1 0 1 2 3]);
yticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_10.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_10.eps'),'-depsc2','-r600');

%% Plotting the snapshots representative of m=1, 2 mode at x/D = 40;

clearvars w_mean w_fluc_sampled;
% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;

% Read w velocity field 

nr = 356; ntheta = 258;
var3 = 'wp';
x = 40;
dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
w = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            w(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
w_centered = w(2:nr-1,2:ntheta-1,:);
w_mean = squeeze(mean(w_centered,3));
w_fluc = w_centered - w_mean;
w_fluc_sampled = w_fluc(:,:,1:end);
clear w_fluc; clear w_centered; clear w;

% Reading the grid files

ntheta = 256;
theta = linspace(0,2*pi,ntheta)';
numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Plotting the m=2 mode at x/D = 40

% axes(ha(2));


load('./files/coefficients_projection_x_D_40.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;



t = 4556;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>


set(h,'Linecolor','none');
colormap hot;
caxis([-0.10 0.10]);

hBar = colorbar;
% set(hBar, 'Location', 'northoutside');
% set(hBar,'Position',[ha(2).Position(1) ha(2).Position(2)+ 0.36 ha(2).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.1:0.05:0.1);
set(hBar, 'TickLabelInterpreter', 'latex');



ax = gca;
ax.FontSize = 20;


xlim([-6 6]);
xticks([-6 -4 -2 0 2 4 6]);
% set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-6 6]);
yticks([-6 -4 -2 0 2 4 6]);
yticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
% hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_40.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_40.eps'),'-depsc2','-r600');


%% Plotting the m=1 mode at x/D = 40

% axes(ha(5));

load('./files/coefficients_projection_x_D_40.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;


t = 761;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);

set(h,'Linecolor','none');
colormap hot;
caxis([-0.10 0.10]);
%colorbar;

ax = gca;
ax.FontSize = 20;

hBar = colorbar;
% set(hBar, 'Location', 'northoutside');
% set(hBar,'Position',[ha(2).Position(1) ha(2).Position(2)+ 0.36 ha(2).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.1:0.05:0.1);
set(hBar, 'TickLabelInterpreter', 'latex');


xlim([-6 6]);
xticks([-6 -4 -2 0 2 4 6]);
xticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});
ylim([-6 6]);
yticks([-6 -4 -2 0 2 4 6]);
yticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
% hTitle  = title('(e)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);



print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_40.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_40.eps'),'-depsc2','-r600');



%% Plotting the snapshots representative of m=1, 2 mode at x/D = 80;

clearvars w_mean w_fluc_sampled;
% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;

% Read w velocity field 

nr = 356; ntheta = 258;
var3 = 'wp';
x = 80;
dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
w = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            w(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
w_centered = w(2:nr-1,2:ntheta-1,:);
w_mean = squeeze(mean(w_centered,3));
w_fluc = w_centered - w_mean;
w_fluc_sampled = w_fluc(:,:,1:end);
clear w_fluc; clear w_centered; clear w;

% Reading the grid files

ntheta = 256;
theta = linspace(0,2*pi,ntheta)';
numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Plotting the m=2 mode at x/D = 80

axes(ha(3));

load('./files/coefficients_projection_x_D_80.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;

t = 5133;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>

set(h,'Linecolor','none');
colormap hot;
caxis([-0.05 0.05]);

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(3).Position(1) ha(3).Position(2)+ 0.36 ha(3).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.05:0.025:0.05);
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 20;

xlim([-9 9]);
xticks([-9 -6 -3 0 3 6 9]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-9 9]);
yticks([-9 -6 -3 0 3 6 9]);
yticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});

hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);
%% Plotting the m=1 mode at x/D = 80

axes(ha(6));

load('./files/coefficients_projection_x_D_80.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);

% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;


t = 1325;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);


set(h,'Linecolor','none');
colormap hot;
caxis([-0.05 0.05]);

ax = gca;
ax.FontSize = 20;

xlim([-9 9]);
xticks([-9 -6 -3 0 3 6 9]);
xticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});
ylim([-9 9]);
yticks([-9 -6 -3 0 3 6 9]);
yticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(f)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);
% latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);


print(gcf,strcat(dirout, 'ux_m12_', '_instant_realizations.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_m12_', '_instant_realizations.eps'),'-depsc2','-r600');

