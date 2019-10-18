%% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';

%% Similarity of m = 2, St = 0 ux mode as a function of downstream distance using Lk

clearvars -except dirout;
close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';


figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ux mode as a function of downstream distance using Ld

clearvars -except dirout;
%close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Lk

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Ld

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');


%% Similarity of m = 2, St = 0 utheta mode as a function of downstream distance using Lk

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 ur mode as a function of downstream distance using Ld

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,1,3,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Lk

clearvars -except dirout;
close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ux mode as a function of downstream distance using Ld

clearvars -except dirout;
%close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(w_eigenmode_allm(:,1,6,2,i)), 'Color',  lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ur mode as a function of downstream distance using Lk

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 ur mode as a function of downstream distance using Ld

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';


figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(u_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 utheta mode as a function of downstream distance using Lk

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.136 utheta mode as a function of downstream distance using Ld

clearvars -except dirout;
% close all;
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
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';

figure;
hold on;
for i = 4:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2),  abs(v_eigenmode_allm(:,1,6,2,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

xlim([0 6]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');


%% Three dimensional slices of real part of ur m=1 St=0.135

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));


% Fix the sign of mode based on the x/D = 10 mode
real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ux_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        disp(i);
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, real(ur_mode(:,i)), '--', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta)); %#ok<*SAGROW>
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ur_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ur_real_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ur_real_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ux m=1 St=0.135

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));

% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, real(ux_mode(:,i)), '--', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ux_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_real_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ur m=1 St=0.135

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));


% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, imag(ur_mode(:,i)), '--', 'Linewidth', 2);
end
close;


for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ur_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ur_imag_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ur_imag_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ux m=1 St=0.135

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,6,2,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,6,2,:));

% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, imag(ux_mode(:,i)), '--', 'Linewidth', 2);
end
close;


for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*1*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*1*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);


count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ux_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ux_imag_m1st0135_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_imag_m1st0135_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ur m=2 St=0

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));


% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, real(ur_mode(:,i)), '--', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta)); %#ok<*SAGROW>
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ur_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ur_real_m2st0_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ur_real_m2st0_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of real part of ux m=2 St=0

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));

% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, real(ux_mode(:,i)), '-', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(real_ux_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1])

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ux_real_m2st0_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_real_m2st0_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ur m=2 St=0

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));


% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, imag(ur_mode(:,i)), '--', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ur_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1])


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ur_imag_m2st0_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ur_imag_m2st0_x_D_slices.eps'),'-depsc2','-r600');

%% Three dimensional slices of imag part of ux m=2 St=0

clearvars -except dirout;
close all;

filename = './files/eigenmodes_similarity_diff_loc.mat';

load(filename);

theta = linspace(0,2*pi,256);

ur_mode = squeeze(u_eigenmode_allm(:,1,1,3,:));
ux_mode = squeeze(w_eigenmode_allm(:,1,1,3,:));

% Fix the sign of mode based on the x/D = 10 mode

real_comp_mode_x10 = real(ux_mode(:,2));

for i = 1:size(ur_mode,2)
    real_comp_mode = real(ux_mode(:,i));
    if real_comp_mode(80,1)/real_comp_mode_x10(80,1) < 0 
        ux_mode(:,i) = -1*ux_mode(:,i);
        ur_mode(:,i) = -1*ur_mode(:,i);
    end
end

figure; hold on;
for i = 2:3:20
    plot(rc, imag(ux_mode(:,i)), '-', 'Linewidth', 2);
end
close;

for i = 1:size(ur_mode,2)
    real_ur_mode_contour(:,:,i) = real(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    real_ux_mode_contour(:,:,i) = real(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

    imag_ur_mode_contour(:,:,i) = imag(ur_mode(:,i)*exp(sqrt(-1)*2*theta));
    imag_ux_mode_contour(:,:,i) = imag(ux_mode(:,i)*exp(sqrt(-1)*2*theta));

end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; 
hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);
hZLabel = zlabel('$x/D$','interpreter','latex','fontsize',15);
view(3);

cameratoolbar('SetCoordSys','x');

box on;
hold on;

ax = gca;
ax.FontSize = 16;

C = cell(7,1);
h = cell(7,1);

count = 1;
for ct = 2:3:20
    [C{count},h{count}] = polarcont(rc,theta',squeeze(imag_ux_mode_contour(:,:,ct)),10);
    axis equal;
    ax = gca;
    ax.Children(1).ContourZLevel = LK_mean_loc_planes(ct,1); %put at correct z

end

colormap('hot');

zlim([0 120]);
ylim([-15 15]);
xlim([-15 15]);

zticks([0 20 40 60 80 100]);
yticks([-10 0 10]);
xticks([-10 0 10]);
camorbit(180,0,'camera');

daspect([1 1 1])

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ux_imag_m2st0_x_D_slices.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_imag_m2st0_x_D_slices.eps'),'-depsc2','-r600');