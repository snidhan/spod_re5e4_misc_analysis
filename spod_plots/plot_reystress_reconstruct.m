%% Name - Sheel Nidhan
%  Date - 10th September, 2019
%  Plots for the reynolds stress reconstruction for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';

%% Plot reconstructed <uxur> at x/D = 20

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 4;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxur_x_D_20', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxur_x_D_20', '.eps'),'-depsc2','-r600');


%% Plot reconstructed <uxur> at x/D = 40

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 8;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 80

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 16;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.eps'),'-depsc2','-r600');



%% Plot reconstructed <uxur> at x/D = 100

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 20;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxur_x_D_100', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxur_x_D_100', '.eps'),'-depsc2','-r600');



%% ------------------------------------------------------------------<uxux>------------------------------------------------------------------------------------------------

%% Plot reconstructed <uxur> at x/D = 20

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 4;

h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*reystress_ww_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), reystress_ww_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_ww_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxux_x_D_20', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxux_x_D_20', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 40

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 8;

h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*reystress_ww_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), reystress_ww_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_ww_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxux_x_D_40', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxux_x_D_40', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 80

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 16;

h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*reystress_ww_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), reystress_ww_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_ww_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxux_x_D_80', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxux_x_D_80', '.eps'),'-depsc2','-r600');



%% Plot reconstructed <uxur> at x/D = 100

clearvars -except dirout;
close all;
filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;

i = 20;

h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*reystress_ww_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,2,i), 'm-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,3,i), 'b-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,4,i), 'm--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*reystress_ww_combined(:,5,i), 'b--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), reystress_ww_1d(:,i), 'k-', 'Linewidth',2);

xlim([0 5]);
ylim([0 1.2*max(abs(reystress_ww_1d(:,i)))]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
                    '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxux_x_D_100', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxux_x_D_100', '.eps'),'-depsc2','-r600');
