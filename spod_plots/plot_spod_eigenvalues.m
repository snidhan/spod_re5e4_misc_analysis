%% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';

%% Plotting the eigenspectrum as a function of x/D

close all;
filename = './files/eigvalue_spectrum_diff_loc.mat';

load(filename);

% m = 1, x/D = 20
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);
figure;
for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  loglog(eigvalues_spectrum_plot(2).freq',eigvalues_spectrum_plot(2).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    hold on;
    grid on;
end

xlim([0 10]);
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2])

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15); %#ok<*NASGU>
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_20.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_20.eps'),'-depsc2','-r600');

% m = 1, x/D = 80
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);
figure;
for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  loglog(eigvalues_spectrum_plot(5).freq',eigvalues_spectrum_plot(5).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    hold on;
    grid on;
end

xlim([0 10]);
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2])

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_80.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_80.eps'),'-depsc2','-r600');

% m = 2, x/D = 20

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);
figure;
for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  loglog(eigvalues_spectrum_plot(3).freq',eigvalues_spectrum_plot(3).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    hold on;
    grid on;
end

xlim([0 10]);
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2])

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_20.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_20.eps'),'-depsc2','-r600');

% m = 2, x/D = 80

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);
figure;
for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  loglog(eigvalues_spectrum_plot(6).freq',eigvalues_spectrum_plot(6).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    hold on;
    grid on;
end

xlim([0 10]);
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2])

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_80.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_80.eps'),'-depsc2','-r600');


%% Contour plot of energy distribution

clearvars -except dirout;
close all;
filename = './files/eigvalue_contourf_diff_loc.mat';
load(filename);

% x/D = 20, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(1).FREQ(1:24,:), eigvalues_contourf_plot(1).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(1).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

xticks([0 1 2 3 4 5]);
colormap jet;
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$St$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',20),'_spod.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',20),'_spod.eps'),'-depsc','-r600');


% x/D = 40, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(2).FREQ(1:24,:), eigvalues_contourf_plot(2).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(2).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

xticks([0 1 2 3 4 5]);
colormap jet;
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$St$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',40),'_spod.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',40),'_spod.eps'),'-depsc','-r600');


% x/D = 40, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(3).FREQ(1:24,:), eigvalues_contourf_plot(3).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(3).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

xticks([0 1 2 3 4 5]);
colormap jet;
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$St$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',60),'_spod.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',60),'_spod.eps'),'-depsc','-r600');

% x/D = 80, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(4).FREQ(1:24,:), eigvalues_contourf_plot(4).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(4).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

xticks([0 1 2 3 4 5]);
colormap jet;
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$St$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',80),'_spod.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',80),'_spod.eps'),'-depsc','-r600');


%% Decay of leading SPOD eigenvalues as a function of x/D for m=1 and m=2

clearvars -except dirout;
close all;
filename = './files/eigvalue_decay_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(x_sample_closest(1:end), squeeze(eigenspectra_allm(6,1,2,1:end)), 'ko');
hold on;
h2 = loglog(x_sample_closest(1:end), squeeze(eigenspectra_allm(1,1,3,1:end)), 'ks');
y1 = 0.04*x_sample_closest(2:3,1).^(coeffs_m1_st0136_55(1));
h3 = loglog(x_sample_closest(2:3,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.10*x_sample_closest(14:20,1).^(coeffs_m1_st0136_70_120(1));
h4 = loglog(x_sample_closest(14:20,1),y2, 'r-', 'Linewidth', 2);
y3 = 0.005*x_sample_closest(4:6,1).^(coeffs_m2_st0(1));
h5 = loglog(x_sample_closest(4:6,1),y3, 'b-', 'Linewidth', 2);
h6 = text(12, 0.003,'$x^{-1.13}$','interpreter','latex','FontSize', 12);
h7 = text(55, 0.00015,'$x^{-1.42}$','interpreter','latex','FontSize', 12);
h8 = text(18, 0.00055,'$x^{-0.60}$','interpreter','latex','FontSize', 12);

xlim([1 125]);
xticks([1 10 100]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2], '$m=1, St=0.136$', '$m=2, St=0$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southwest';

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.eps'),'-depsc2','-r600');

%% Similarity analysis of the m=2, St=0 leading order SPOD mode

clearvars -except dirout;
close all;
filename = './files/eigvalues_similarity_diff_loc.mat';
load(filename);

C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

count = 1;

Legend = cell(10,1);
Legend{1} = 'x/D = 10';
Legend{2} = 'x/D = 20';
Legend{3} = 'x/D = 30';
Legend{4} = 'x/D = 40';
Legend{5} = 'x/D = 50';
Legend{6} = 'x/D = 60';
Legend{7} = 'x/D = 70';
Legend{8} = 'x/D = 80';
Legend{9} = 'x/D = 90';
Legend{10} = 'x/D = 100';


count = 1;
Nf_sampled = 100;
figure;
hold on;
for i = 2:2:20
   disp(i);
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,3,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count});
   count = count + 1;
end

xlim([0 0.5]);
ylim([0 0.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}(m=2, St=0)/(K^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.eps'),'-depsc2','-r600');


%% Similarity analysis of the m=1, St=0.136 leading order SPOD mode

clearvars -except dirout;
close all;
filename = './files/eigvalues_similarity_diff_loc.mat';
load(filename);

C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

count = 1;

Legend = cell(10,1);
Legend{1} = 'x/D = 10';
Legend{2} = 'x/D = 20';
Legend{3} = 'x/D = 30';
Legend{4} = 'x/D = 40';
Legend{5} = 'x/D = 50';
Legend{6} = 'x/D = 60';
Legend{7} = 'x/D = 70';
Legend{8} = 'x/D = 80';
Legend{9} = 'x/D = 90';
Legend{10} = 'x/D = 100';


count = 1;
Nf_sampled = 100;
figure;
hold on;
for i = 2:2:20
   disp(i);
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,2,i)/((ud_centerline_loc_planes(i,2))*LK_mean_loc_planes(i,2))^2, C{count});
   count = count + 1;
end

xlim([0 0.5]);
ylim([0 0.25]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$St$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}(m=1, St=0.136)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_m1_st0136_eigvalue_x_D_udld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_m1_st0136_eigvalue_x_D_udld', '.eps'),'-depsc2','-r600');

%% Bar chart plot for the dominance shift from m = 2 to m = 1 mode 

clearvars -except dirout;
close all;
filename = './files/eigvalue_barplot_diff_loc.mat';
load(filename);

% x/D = 20, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

color = [0 0 0; 1 0 0; 0 0 1];
colormap(color);

h1 = bar(eigvalues_bar_plot(1).mode, eigvalues_bar_plot(1).eigenvalue_fraction_contri(1:3,:)','grouped');

ylim([0 20]);

% set(h1(1),'FaceColor', 'c');
% set(h1(2),'FaceColor', 'b');
% set(h1(3),'FaceColor', 'm');

ax = gca;
ax.FontSize = 16; 

labels = {'SPOD Mode 1','SPOD Mode 2','SPOD Mode 3'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',20),'_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',20),'_spod.eps'),'-depsc2','-r600');

% x/D = 40, bar plot

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(2).mode, eigvalues_bar_plot(2).eigenvalue_fraction_contri(1:3,:)','grouped');

ylim([0 20]);

% set(h1(1),'FaceColor', 'c');
% set(h1(2),'FaceColor', 'b');
% set(h1(3),'FaceColor', 'm');

ax = gca;
ax.FontSize = 16; 

labels = {'SPOD Mode 1','SPOD Mode 2','SPOD Mode 3'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\xi^{(1)}$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',40),'_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',40),'_spod.eps'),'-depsc2','-r600');

% x/D = 60, bar plot

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(3).mode, eigvalues_bar_plot(3).eigenvalue_fraction_contri(1:3,:)','grouped');

ylim([0 20]);

% set(h1(1),'FaceColor', 'c');
% set(h1(2),'FaceColor', 'b');
% set(h1(3),'FaceColor', 'm');

ax = gca;
ax.FontSize = 16; 

labels = {'SPOD Mode 1','SPOD Mode 2','SPOD Mode 3'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',60),'_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',60),'_spod.eps'),'-depsc2','-r600');

% x/D = 80, Contourmap

figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(4).mode, eigvalues_bar_plot(4).eigenvalue_fraction_contri(1:3,:)','grouped');

ylim([0 20]);

% set(h1(1),'FaceColor', 'c');
% set(h1(2),'FaceColor', 'b');
% set(h1(3),'FaceColor', 'm');

ax = gca;
ax.FontSize = 16; 

labels = {'SPOD Mode 1','SPOD Mode 2','SPOD Mode 3'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',80),'_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',80),'_spod.eps'),'-depsc2','-r600');
