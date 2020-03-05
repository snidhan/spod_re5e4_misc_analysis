%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting the decay of eigenvalues 

clear;
addpath('./aux_plots/');
dirout = './';

%% Figure of contour plot of energy distribution

filename = './aux_plots/files/eigvalue_decay_diff_loc.mat';
load(filename);

%% Setting up the figure environment

figure;

close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% Plot m=1, St=0.135 
h1 = loglog(x_sample_closest(1:end-1), squeeze(eigenspectra_allm(6,1,2,1:end-1)), 'ko', 'MarkerSize', 7);
hold on;

loc = x_sample_closest(2:end-1);                         
log_mode1 = log(squeeze(eigenspectra_allm(6,1,2,2:end-1)));
log_loc =  log(loc);
[coeffs_m1_st0135, S] = polyfit(log_loc, log_mode1, 1);
y_fitted = polyval(coeffs_m1_st0135, log_loc);

y1 = 0.035*x_sample_closest(2:end-1).^(coeffs_m1_st0135(1));
h2 = loglog(x_sample_closest(2:end-1,1),y1, 'r--', 'Linewidth', 2);
%% Plot m=2, St=0
h3 = loglog(x_sample_closest(1:end-1), squeeze(eigenspectra_allm(1,1,3,1:end-1)), 'ks','MarkerSize', 7);

loc = x_sample_closest(2:end-1);                         
log_mode1 = log(squeeze(eigenspectra_allm(1,1,3,2:end-1)));
log_loc =  log(loc);
[coeffs_m2_st0, S] = polyfit(log_loc, log_mode1, 1);
y_fitted = polyval(coeffs_m2_st0, log_loc);

y3 = 0.0054*x_sample_closest(2:end-1,1).^(coeffs_m2_st0(1));
h4 = loglog(x_sample_closest(2:end-1,1),y3, 'b--', 'Linewidth', 2);

%h6 = text(12, 0.003,'$x^{-1.13}$','interpreter','latex','FontSize', 12);
%h7 = text(55, 0.00015,'$x^{-1.42}$','interpreter','latex','FontSize', 12);
%h8 = text(18, 0.00055,'$x^{-0.60}$','interpreter','latex','FontSize', 12);

%% Decay of MKE 

filename = './aux_plots/files/MKE_areaI_FINF.dat';

mke_integrated = importdata(filename);

loc = mke_integrated(11:38,1);                         
log_mke = log(mke_integrated(11:38,2));
log_loc =  log(loc);
[coeffs_mke, S] = polyfit(log_loc, log_mke, 1)
y_fitted = polyval(coeffs_mke, log_loc);

figure;
loglog(mke_integrated(:,1), mke_integrated(:,2),'ko')
hold on;
%% Decay of TKE 

filename = './aux_plots/files/TKE_areaI_FINF.dat';

tke_integrated = importdata(filename);

loc = tke_integrated(11:38,1);                         
log_tke = log(tke_integrated(11:38,2));
log_loc =  log(loc);
[coeffs_tke, S] = polyfit(log_loc, log_tke, 1)
y_fitted = polyval(coeffs_tke, log_loc);

loglog(tke_integrated(:,1), tke_integrated(:,2),'ro')


%% Setting up more figure environments
xlim([1 125]);
xticks([1 10 100]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4], '$m=1, \mbox{\textit{St}}=0.135$ (VS)', '$x^{-1.14}$', '$m=2, \mbox{\textit{St}}=0$ (DH)', '$x^{-0.60}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southwest';

%% Saving the eigenvalues
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.eps'),'-depsc2','-r600');
