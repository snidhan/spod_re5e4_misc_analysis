%% Name - Sheel Nidhan
%  Date - 10th September, 2019
%  Plots for the reynolds stress reconstruction for prf_spod_re5e4_frinf paper

% dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures_2.0/';
%% --------------------------------------------------------Plot reconstructed <-uxur> at x/D = 20, 40, 80, 100-----------------------------------------------------------------------------------------------------------

clearvars -except dirout;
close all;

filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;
x0=5;
y0=5;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(2,2,[0.1, 0.08],[.2,.05],[.1 .02]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


axes(ha(1));

i = 4;

total_sum = -2*reystress_uw_combined(:,1,i) - 4*reystress_uw_combined(:,2,i) - 4*reystress_uw_combined(:,3,i) ...
            -4*reystress_uw_combined(:,4,i) - 4*reystress_uw_combined(:,5,i);
hold on;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);
yticks([0 0.25*max(abs(reystress_uw_1d(:,i))) 0.5*max(abs(reystress_uw_1d(:,i))) 0.75*max(abs(reystress_uw_1d(:,i))) 1*max(abs(reystress_uw_1d(:,i)))])
yticklabels({'0','0.24','0.48', '0.72', '0.96'});

xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [0.1 0.95 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize', 20, 'Position', [-0.8, 0.0006, 0]); %#ok<*NASGU>
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

hLegend = legend([h1, h2, h3, h4, h5, h6, h7], '(1) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = 0$', ...
                                           '(2) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 1$', ...
                                           '(3)$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 2$', ...
                                           '(4) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 3$', ...
                                           '(5)  $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 4$', ...
                                           'Ensemble averaged', ...
                                           'Summation from (1) to (5) ');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

%% Plot reconstructed <uxur> at x/D = 40

axes(ha(2));

i = 8;

total_sum = -2*reystress_uw_combined(:,1,i) - 4*reystress_uw_combined(:,2,i) - 4*reystress_uw_combined(:,3,i) ...
            -4*reystress_uw_combined(:,4,i) - 4*reystress_uw_combined(:,5,i);

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);
yticks([0 0.25*max(abs(reystress_uw_1d(:,i))) 0.5*max(abs(reystress_uw_1d(:,i))) 0.75*max(abs(reystress_uw_1d(:,i))) 1*max(abs(reystress_uw_1d(:,i)))])
yticklabels({'0','0.07','0.15', '0.23', '0.31'});

xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [0.58 0.95 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 80


axes(ha(3));


i = 16;

total_sum = -2*reystress_uw_combined(:,1,i) - 4*reystress_uw_combined(:,2,i) - 4*reystress_uw_combined(:,3,i) ...
            -4*reystress_uw_combined(:,4,i) - 4*reystress_uw_combined(:,5,i);

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);
yticks([0 0.25*max(abs(reystress_uw_1d(:,i))) 0.5*max(abs(reystress_uw_1d(:,i))) 0.75*max(abs(reystress_uw_1d(:,i))) 1*max(abs(reystress_uw_1d(:,i)))])
yticklabels({'0','0.28','0.56', '0.84', '1.1'});

xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4', '5'});

str = '$\times 10^{-4}$';
annotation('textbox', [0.1 0.525 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize',20);
hTitle  = title('(c)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 100

axes(ha(4));

i = 20;

total_sum = -2*reystress_uw_combined(:,1,i) - 4*reystress_uw_combined(:,2,i) - 4*reystress_uw_combined(:,3,i) ...
            -4*reystress_uw_combined(:,4,i) - 4*reystress_uw_combined(:,5,i);

hold on;

h1 = plot(rc/LK_TKE_loc_planes(i,2), -2*reystress_uw_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), -4*reystress_uw_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(reystress_uw_1d(:,i)))]);
yticks([0 0.25*max(abs(reystress_uw_1d(:,i))) 0.5*max(abs(reystress_uw_1d(:,i))) 0.75*max(abs(reystress_uw_1d(:,i))) 1*max(abs(reystress_uw_1d(:,i)))])
yticklabels({'0','0.20','0.40', '0.60', '0.80'});

xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4', '5'});

str = '$\times 10^{-4}$';
annotation('textbox', [0.58 0.525 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
%hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_uxur_x_D_20_40_80_100', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_uxur_x_D_20_40_80_100', '.eps'),'-depsc2','-r600');
%% --------------------------------------------------------Plot reconstructed TKE at x/D = 20, 40, 80, 100-----------------------------------------------------------------------------------------------------------

%% Plot of K at x/D = 20

clearvars -except dirout;
close all;

filename = './files/reystress_in_ranseq_construct_similarity_diff_loc.mat';
load(filename);

figure;
x0=5;
y0=5;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(2,2,[0.1, 0.08],[.2,.05],[.1 .02]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


axes(ha(1));

i = 4;

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
yticks([0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))])
yticklabels({'0','1.89','3.78','5.67', '7.56', '10.05'});

xlim([0.03 5]);
xticks([0 1 2 3 4 5]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [0.1 0.95 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$K$','interpreter','latex','fontsize', 20, 'Position', [-0.8, 0.005, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = 0$', ...
                                           '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 1$', ...
                                           '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 2$', ...
                                           '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 3$', ...
                                           '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 4$', ...
                                           'From simulation');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

%% Plot or K at x/D = 40

axes(ha(2));

i = 8;

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
yticks([0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))])
yticklabels({'0','0.65','1.30', '1.95', '2.60'});

xlim([0.03 5]);
xticks([0 1 2 3 4 5]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [0.58 0.95 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.eps'),'-depsc2','-r600');


%% Plot of K at x/D = 80


axes(ha(3));

hold on;

i = 16;

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);


ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
yticks([0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))])
[0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))];
yticklabels({'0','2.51','5.03', '7.54', '10.05'});

xlim([0.03 5]);
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});

str = '$\times 10^{-4}$';
annotation('textbox', [0.1 0.525 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K$','interpreter','latex','fontsize',20, 'Position', [-0.8, 0.0007, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.eps'),'-depsc2','-r600');

%% Plot of K at x/D = 100

axes(ha(4));

i = 20;

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
yticks([0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))])
[0 0.25*max(abs(total_averaged_tke_1d(:,i))) 0.5*max(abs(total_averaged_tke_1d(:,i))) 0.75*max(abs(total_averaged_tke_1d(:,i))) 1*max(abs(total_averaged_tke_1d(:,i)))]
yticklabels({'0','1.91','3.81', '5.72', '7.63'});

xlim([0.03 5]);
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});

str = '$\times 10^{-4}$';
annotation('textbox', [0.58 0.525 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
%hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'recon_k_x_D_20_40_80_100', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'recon_k_x_D_20_40_80_100', '.eps'),'-depsc2','-r600');