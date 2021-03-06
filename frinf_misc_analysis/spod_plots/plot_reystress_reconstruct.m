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

%%
figure;
x0=5;
y0=5;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(2,2,[0.05, 0.08],[.08,.05],[.12 .02]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%
% axes(ha(1));

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
tick_val = [0; 0.25*max(abs(reystress_uw_1d(:,i))); 0.5*max(abs(reystress_uw_1d(:,i))); 0.75*max(abs(reystress_uw_1d(:,i))); 1*max(abs(reystress_uw_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 1000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(3,1)),2), num2str(sprintf('%.2f', 1000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(5,1)),2)});

xlim([0 3]);
xticks([0 1 2 3]);
% set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
%annotation('textbox', [ha(1).Position(1) ha(1).Position(2)+ha(1).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
annotation('textbox', [pt1{1}(1) pt1{1}(2)+pt1{1}(4)+0.05 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',15, 'Edgecolor', 'none');

% ax.FontSize = 20; 
ax.FontSize = 15; 


box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize', 15); %#ok<*NASGU>

% hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize', 15, 'Position', [-0.6, 0.0006, 0]); %#ok<*NASGU>
% hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend([h1, h2, h3, h4, h5, h6, h7], '(1) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = 0$', ...
%                                            '(2) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 1$', ...
%                                            '(3)$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 2$', ...
%                                            '(4) $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 3$', ...
%                                            '(5)  $\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 4$', ...
%                                            'Measured', ...
%                                            'Summation from (1) to (5) ');

hLegend = legend([h1, h2, h3, h4, h5, h6, h7], '(1) $m = 0$', ...
                                           '(2) $m = \pm 1$', ...
                                           '(3) $m = \pm 2$', ...
                                           '(4) $m = \pm 3$', ...
                                           '(5) $m = \pm 4$', ...
                                           'Measured', ...
                                           'Summed over (1) to (5) ');

                                       
                                       
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'recon_uxur_x_D_20', '.png'),'-dpng2','-r600');  


%% Plot reconstructed <uxur> at x/D = 40

% axes(ha(2));

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
tick_val = [0; 0.25*max(abs(reystress_uw_1d(:,i))); 0.5*max(abs(reystress_uw_1d(:,i))); 0.75*max(abs(reystress_uw_1d(:,i))); 1*max(abs(reystress_uw_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 1000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(3,1)),2), num2str(sprintf('%.2f', 1000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(5,1)),2)});

xlim([0 3]);
xticks([0 1 2 3]);
% set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
% annotation('textbox', [ha(2).Position(1) ha(2).Position(2)+ha(2).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
annotation('textbox', [pt1{1}(1) pt1{1}(2)+pt1{1}(4)+0.05 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',15, 'Edgecolor', 'none');

% ax.FontSize = 20; 
ax.FontSize = 15; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$$\langle -u^{''}_{x}u^{''}_{r} \rangle$$','interpreter','latex','fontsize',15);
% hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'recon_uxur_x_D_40', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_40', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 80


% axes(ha(3));


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
tick_val = [0; 0.25*max(abs(reystress_uw_1d(:,i))); 0.5*max(abs(reystress_uw_1d(:,i))); 0.75*max(abs(reystress_uw_1d(:,i))); 1*max(abs(reystress_uw_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(10000*tick_val(1,1)), num2str(sprintf('%.2f', 10000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(3,1)),2), num2str(sprintf('%.2f', 10000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(5,1)),2)});

xlim([0 3]);
xticks([0 1 2 3]);
xticklabels({'0','1','2', '3'});

str = '$\times 10^{-4}$';
% annotation('textbox',  [ha(3).Position(1) ha(3).Position(2)+ha(3).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
annotation('textbox', [pt1{1}(1) pt1{1}(2)+pt1{1}(4)+0.05 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',15, 'Edgecolor', 'none');

ax.FontSize = 15; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize',20, 'Position', [-0.6, 0.00007, 0]);
hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize',15);

% hTitle  = title('(c)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'recon_uxur_x_D_80', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_80', '.eps'),'-depsc2','-r600');

%% Plot reconstructed <uxur> at x/D = 100

% axes(ha(4));

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
tick_val = [0; 0.25*max(abs(reystress_uw_1d(:,i))); 0.5*max(abs(reystress_uw_1d(:,i))); 0.75*max(abs(reystress_uw_1d(:,i))); 1*max(abs(reystress_uw_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(10000*tick_val(1,1)), num2str(sprintf('%.2f', 10000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(3,1)),2), num2str(sprintf('%.2f', 10000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(5,1)),2)});

xlim([0 3]);
xticks([0 1 2 3]);
xticklabels({'0','1','2', '3'});

str = '$\times 10^{-4}$';
% annotation('textbox',  [ha(4).Position(1) ha(4).Position(2)+ha(4).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

AxesHandle=findobj(gcf,'Type','axes'); % For PPT
pt1 = get(AxesHandle,{'Position','tightinset','PlotBoxAspectRatio'});
annotation('textbox', [pt1{1}(1) pt1{1}(2)+pt1{1}(4)+0.05 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',15, 'Edgecolor', 'none');


ax.FontSize = 15; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
%hYLabel = ylabel('-$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\langle -u^{''}_{x}u^{''}_{r} \rangle$','interpreter','latex','fontsize',15);

% hTitle  = title('(d)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{f = 0 \rightarrow 0.5 } m = 0$', '$\sum_{f = 0 \rightarrow 0.5 } m =1$','$\sum_{f = 0 \rightarrow 0.5 } m = 2$', ...
%                     '$\sum_{f = 0 \rightarrow 0.5 } m = 3$', '$\sum_{f = 0 \rightarrow 0.5 } m = 4$', 'Averaged over 528 time units');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'northeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'recon_uxur_x_D_100', '.png'),'-dpng2','-r600');  

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'recon_uxur_x_D_20_40_80_100', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'recon_uxur_x_D_20_40_80_100', '.eps'),'-depsc2','-r600');
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

[ha, pos] = tight_subplot(2,2,[0.05, 0.08],[.08,.05],[.12 .02]);


set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


axes(ha(1));

i = 4;


total_sum = 2*tke_combined(:,1,i)  + 4*tke_combined(:,2,i) + 4*tke_combined(:,3,i) ...
            + 4*tke_combined(:,4,i) + 4*tke_combined(:,5,i);

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
tick_val = [0; 0.25*max(abs(total_averaged_tke_1d(:,i))); ...
            0.5*max(abs(total_averaged_tke_1d(:,i))); 0.75*max(abs(total_averaged_tke_1d(:,i))); 1*max(abs(total_averaged_tke_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 1000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(3,1)),2), num2str(sprintf('%.2f', 1000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(5,1)),2)});

xlim([0.03 3]);
xticks([0 1 2 3]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [ha(1).Position(1) ha(1).Position(2)+ha(1).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$K$','interpreter','latex','fontsize', 20, 'Position', [-0.6, 0.005, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend([h1, h2, h3, h4, h5, h6], '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = 0$', ...
%                                            '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 1$', ...
%                                            '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 2$', ...
%                                            '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 3$', ...
%                                            '$\sum_{ \mbox{\textit{St}} = -0.5 \rightarrow 0.5 } m = \pm 4$', ...
%                                            'From simulation');

hLegend = legend([h1, h2, h3, h4, h5, h6, h7], '(1) $m = 0$', ...
                                           '(2) $m = \pm 1$', ...
                                           '(3) $m = \pm 2$', ...
                                           '(4) $m = \pm 3$', ...
                                           '(5) $m = \pm 4$', ...
                                           'Measured', ...
                                           'Summed over (1) to (5) ');
                                       
                                       
                                       
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

%% Plot or K at x/D = 40

axes(ha(2));

i = 8;

total_sum = 2*tke_combined(:,1,i)  + 4*tke_combined(:,2,i) + 4*tke_combined(:,3,i) ...
            + 4*tke_combined(:,4,i) + 4*tke_combined(:,5,i);

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
tick_val = [0; 0.25*max(abs(total_averaged_tke_1d(:,i))); ...
            0.5*max(abs(total_averaged_tke_1d(:,i))); 0.75*max(abs(total_averaged_tke_1d(:,i))); 1*max(abs(total_averaged_tke_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 1000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(3,1)),2), num2str(sprintf('%.2f', 1000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 1000*tick_val(5,1)),2)});

xlim([0.03 3]);
xticks([0 1 2 3]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.

str = '$\times 10^{-3}$';
annotation('textbox', [ha(2).Position(1) ha(2).Position(2)+ha(2).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

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

total_sum = 2*tke_combined(:,1,i)  + 4*tke_combined(:,2,i) + 4*tke_combined(:,3,i) ...
            + 4*tke_combined(:,4,i) + 4*tke_combined(:,5,i);

h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);


ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
tick_val = [0; 0.25*max(abs(total_averaged_tke_1d(:,i))); ...
            0.5*max(abs(total_averaged_tke_1d(:,i))); 0.75*max(abs(total_averaged_tke_1d(:,i))); 1*max(abs(total_averaged_tke_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 10000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(3,1)),2), num2str(sprintf('%.2f', 10000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(5,1)),2)});

xlim([0.03 3]);
xticks([0 1 2 3]);
xticklabels({'0','1','2', '3'});

str = '$\times 10^{-4}$';
annotation('textbox', [ha(3).Position(1) ha(3).Position(2)+ha(3).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

ax.FontSize = 20; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K$','interpreter','latex','fontsize',20, 'Position', [-0.6, 0.00065, 0]);
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

total_sum = 2*tke_combined(:,1,i)  + 4*tke_combined(:,2,i) + 4*tke_combined(:,3,i) ...
            + 4*tke_combined(:,4,i) + 4*tke_combined(:,5,i);

hold on;
h1 = plot(rc/LK_TKE_loc_planes(i,2), 2*tke_combined(:,1,i), 'c-', 'Linewidth',2);
h2 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,2,i), 'b-', 'Linewidth',2);
h3 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,3,i), 'r-', 'Linewidth',2);
h4 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,4,i), 'b--', 'Linewidth',2);
h5 = plot(rc/LK_TKE_loc_planes(i,2), 4*tke_combined(:,5,i), 'r--', 'Linewidth',2);
h6 = plot(rc/LK_TKE_loc_planes(i,2), total_averaged_tke_1d(:,i), 'k-', 'Linewidth',2);
h7 = plot(rc/LK_TKE_loc_planes(i,2),  total_sum, 'k--', 'Linewidth',2);

ax = gca;

ylim([0 1.2*max(abs(total_averaged_tke_1d(:,i)))]);
tick_val = [0; 0.25*max(abs(total_averaged_tke_1d(:,i))); ...
            0.5*max(abs(total_averaged_tke_1d(:,i))); 0.75*max(abs(total_averaged_tke_1d(:,i))); 1*max(abs(total_averaged_tke_1d(:,i)))];
yticks(tick_val);
tick_char = num2str(tick_val);
yticklabels({num2str(1000*tick_val(1,1)), num2str(sprintf('%.2f', 10000*tick_val(2,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(3,1)),2), num2str(sprintf('%.2f', 10000*tick_val(4,1)),2), ...
             num2str(sprintf('%.2f', 10000*tick_val(5,1)),2)});
         
xlim([0.03 3]);
xticks([0 1 2 3]);
xticklabels({'0','1','2', '3'});

str = '$\times 10^{-4}$';
annotation('textbox', [ha(4).Position(1) ha(4).Position(2)+ha(4).Position(4) 0.03 0.03], 'interpreter','latex', 'String', str, 'FontSize',20, 'Edgecolor', 'none');

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