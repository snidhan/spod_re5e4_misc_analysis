%% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

clear;
dirout = './';
addpath('./aux_plots/');
%% Generating the figure object

close all;
x0=0;
y0=0;
width=15;
height=15;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(3,2,[.02, 0.05],[.1,.05],[.1 .05]);
%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Lk

axes(ha(1));

filename = './aux_plots/files/eigenmodes_similarity_diff_loc.mat';
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
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(u_eigenmode_allm(:,1,1,3,i)), 'loess', 2), ...
       'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ylim([0 1.2]);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.

ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.12, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.10, 0.80,0]);


%% Legends

clear hLegend;
hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Orientation = 'horizontal';
hLegend.Position = [0.025 0.48 1 1];

%% Inset plot

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(u_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

ax2 = axes('Position',[ha(1).Position(1)+0.22 ha(1).Position(2)+0.04 .22 .22]); % For paper

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT 
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');



%% Similarity of m = 2, St = 0 utheta mode as a function of downstream distance using Lk

axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(v_eigenmode_allm(:,1,1,3,i)), 'loess', 2), ...
        'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20;
xlim([0 6]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

%hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','i`nterpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.12, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.10, 0.80, 0]);


% Inset plot

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(v_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

ax2 = axes('Position',[ha(3).Position(1)+0.22 ha(3).Position(2)+0.04 .22 .22]);

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);


[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');
set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_v_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ux mode as a function of downstream distance using Lk

axes(ha(5));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(w_eigenmode_allm(:,1,1,3,i)), 'loess', 2), ...
       'Color', lineStyles(count,:), 'Linewidth',2);
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
hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20, 'Position', [3, -0.15, 0]);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.12, 0.45, 0]);
hTitle  = title('(e)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.10, 0.80, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(w_eigenmode_allm(:,1,1,3,8)).*exp(sqrt(-1)*2.*theta);

ax2 = axes('Position',[ha(5).Position(1)+0.22 ha(5).Position(2)+0.04 .22 .22]);

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.27 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ur mode as a function of downstream distance using Lk

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(u_eigenmode_allm(:,1,6,2,i)), 'loess', 2), ...
        'Color',  lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
set(gca, 'Xticklabel', []);
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca, 'Yticklabel', []);

box on;

% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.80, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(u_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

ax2 = axes('Position',[ha(2).Position(1)+0.22 ha(2).Position(2)+0.04 .22 .22]);

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;


% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 utheta mode as a function of downstream distance using Ld

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(v_eigenmode_allm(:,1,6,2,i)), 'loess',2), ...
       'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 20; 

xlim([0 6]);
xticks([0 2 4 6]);
set(gca, 'Xticklabel', []);
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca, 'Yticklabel', []);

box on;

% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.80, 0]);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(v_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

ax2 = axes('Position',[ha(4).Position(1)+0.22 ha(4).Position(2)+0.04 .22 .22]);

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');set(h,'Linecolor','none');
set(h,'Linecolor','none');

ax = gca;
ax.FontSize = 15;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_v_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Ld

axes(ha(6));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle


lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  smoothdata(abs(w_eigenmode_allm(:,1,6,2,i)), 'loess',2), ...
       'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end


ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
xticklabels({'$0$', '$2$', '$4$', '$6$'});
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca, 'Yticklabel', []);

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20, 'Position', [3, -0.15, 0]);
hTitle  = title('(f)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.80, 0]);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

theta = linspace(0,2*pi,256);
two_dimensional_mode_shape_x_D_50 = squeeze(w_eigenmode_allm(:,1,6,2,8)).*exp(sqrt(-1)*1.*theta);

ax2 = axes('Position',[ha(6).Position(1)+0.22 ha(6).Position(2)+0.04 .22 .22]);

% AxesHandle=findobj(gcf,'Type','axes'); % For PPT
% pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
% ax2 = axes('Position',[pt1{1}(1)+0.3 pt1{1}(2)+0.3 .45 .45]);

[C,h] = polarcont(rc,theta',squeeze(real(two_dimensional_mode_shape_x_D_50)),10);
colormap('hot');
set(h,'Linecolor','none');


ax = gca;
ax.FontSize = 20;

ylim([-13 13]);
yticks([-10 0 10]);
% yticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Yticklabel', []);


xlim([-13 13]);
xticks([-10 0 10]);
% xticklabels({'$-10$', '$0$', '$10$'});
set(gca, 'Xticklabel', []);

axis equal;

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

box on;

%% Saving plot
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_m1st0136_m2st0_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_m1st0136_m2st0_x_D', '.eps'),'-depsc2','-r600');

%% Saving annotated plots

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_annotated_m1st0136_m2st0_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uvw_eigenmode_annotated_m1st0136_m2st0_x_D', '.eps'),'-depsc2','-r600');
