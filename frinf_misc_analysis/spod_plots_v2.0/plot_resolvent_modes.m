%% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

clear;
dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures_2.0/';
%% Generating the figure object

close all;
x0=0;
y0=0;
width=15;
height=15;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(3,2,[.02, 0.05],[.1,.05],[.1 .05]);


%% m = 2, St = 0 ur mode as a function of downstream distance using Ld

axes(ha(1));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m2_kx0_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m2.mat';
load(filename_spod, 'u_eigenmode', 'rc');

% mode_plot = smoothdata(abs(u_eigenmode(:,1,1)),'gaussian'); 
mode_plot = abs(u_eigenmode(:,1,1));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(N+1:2*N,1,1,1))/max(max(abs(SU(N+1:2*N,1,1,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

ylim([0 1.2])
xlim([0 6]);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

box on;

% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0,0)$', '$1^{st}$ SPOD mode; $St=0$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 2, St = 0 utheta mode as a function of downstream distance using Lk

axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m2_kx0_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m2.mat';
load(filename_spod, 'v_eigenmode', 'rc');

% mode_plot = smoothdata(abs(v_eigenmode(:,1,1)),'lowess'); 
mode_plot = abs(v_eigenmode(:,1,1));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(2*N+1:3*N,1,1,1))/max(max(abs(SU(2*N+1:3*N,1,1,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0,0)$', '$1^{st}$ SPOD mode; $St=0$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','i`nterpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

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

filename = './files/m2_kx0_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m2.mat';
load(filename_spod, 'w_eigenmode', 'rc');

% mode_plot = smoothdata(abs(v_eigenmode(:,1,1)),'lowess'); 
mode_plot = abs(w_eigenmode(:,1,1));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(1:N,1,1,1))/max(max(abs(SU(1:N,1,1,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

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
hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.45, 0]);
hTitle  = title('(e)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0,0)$', '$1^{st}$ SPOD mode; $St=0$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% m=1, St=0.135 plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% m = 1, St = 0.135 ur mode as a function of downstream distance using Ld

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m1_kx0.8482_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m1.mat';
load(filename_spod, 'u_eigenmode', 'rc');

mode_plot = smoothdata(abs(u_eigenmode(:,1,6)),'lowess'); 
% mode_plot = abs(u_eigenmode(:,1,1));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(N+1:2*N,1,6,1))/max(max(abs(SU(N+1:2*N,1,6,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

ylim([0 1.2])
xlim([0 6]);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.

box on;

% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0.8482,0.135)$', '$1^{st}$ SPOD mode; $St=0.135$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_w_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.135 utheta mode as a function of downstream distance using Lk

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m1_kx0.8482_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m1.mat';
load(filename_spod, 'v_eigenmode', 'rc');

mode_plot = smoothdata(abs(v_eigenmode(:,1,6)),'lowess'); 
% mode_plot = abs(v_eigenmode(:,1,1));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(2*N+1:3*N,1,6,1))/max(max(abs(SU(2*N+1:3*N,1,6,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

ylim([0 1.2])
xlim([0 6]);

ax = gca;
ax.FontSize = 20;

xlim([0 19]);
xticks([0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.

hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0.8482,0.135)$', '$1^{st}$ SPOD mode; $St=0.135$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{r}|/|\Phi_{r}|_{\infty}$','i`nterpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{\theta}|/|\Phi_{\theta}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'similarity_v_eigenmode_m2_st0_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'similarity_u_eigenmode_m2_st0_x_D_lk', '.eps'),'-depsc2','-r600');

%% Similarity of m = 1, St = 0.135 ux mode as a function of downstream distance using Lk

axes(ha(6));
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m1_kx0.8482_x_D_100.mat';

load(filename);

filename_spod = './files/eigenmodes_x_D_100_m1.mat';
load(filename_spod, 'w_eigenmode', 'rc');

% mode_plot = smoothdata(abs(v_eigenmode(:,1,1)),'lowess'); 
mode_plot = abs(w_eigenmode(:,1,6));
% mode_plot(1,1) = 0;

hold on;
h1 = plot(r/Ld, squeeze(abs(SU(1:N,1,6,1))/max(max(abs(SU(1:N,1,6,1))))), 'k-', 'Linewidth',2);
h2 = plot(rc/Ld, squeeze(mode_plot)/max(max(mode_plot)), 'r-', 'Linewidth',2);

ax = gca;
ax.FontSize = 20;

xlim([0 6]);
xticks([0 2 4 6]);
xticklabels({'$0$', '$2$', '$4$', '$6$'});
ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.

box on;

% hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',20);
hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$|\Phi_{x}|/|\Phi_{x}|_{\infty}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.45, 0]);
hTitle  = title('(f)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

hLegend = legend([h1,h2], '$\sigma_{1}$ mode;$(k_{x},St)=(0.8482,0.135)$', '$1^{st}$ SPOD mode; $St=0.135$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 12;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';

box on;

%% Saving the plots

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'resolvent_modes_m12_x_D_100.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'resolvent_modes_m12_x_D_100.eps'),'-depsc2','-r600');
