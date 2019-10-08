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
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

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
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

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
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_u_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');


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

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

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
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_v_eigenmode_m2_st0_x_D_ld', '.eps'),'-depsc2','-r600');

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
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_w_eigenmode_m1_st0136_x_D_ld', '.eps'),'-depsc2','-r600');

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
