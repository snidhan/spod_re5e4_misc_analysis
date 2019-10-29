%% Name - Sheel Nidhan
%  Date - 11th September, 2019
%  Plots for the eigenvalues results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures_2.0/';
% dirout = 'C:\Users\snidh\Dropbox\research\sheel_papers\prf_spod_re5e4_frinf\template\figures_2.0\';
%% Plotting the eigenspectrum as a function of x/D

close all;
filename = './files/eigvalue_spectrum_diff_loc.mat';
load(filename);
close all;


x0=0;
y0=0;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,2,[.05, 0.05],[.15,.05],[.15 .05]);
%% m = 1, x/D = 20

axes(ha(1));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);

first_spod_mode  = eigvalues_spectrum_plot(2).eigvalue(:,2);
second_spod_mode = eigvalues_spectrum_plot(2).eigvalue(:,1);

hold on;
box on;

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    eigvalue = eigvalues_spectrum_plot(2).eigvalue(:,i);
    f = eigvalues_spectrum_plot(2).freq';
    idx = f>0 & eigvalue>0;
    h1(i) =  plot(f(idx), eigvalue(idx), 'LineWidth',2,'Color',grey);
    grid on;
end
patch([f(idx)' fliplr(f(idx)')], [second_spod_mode(idx)'  fliplr(first_spod_mode(idx)')], 'r');
hold off;
set(gca, 'XScale', 'log', 'YScale','log')

ax = gca;
ax.FontSize = 20;

xlim([2.7*10^-2 6]);
xticks([0.1 0.5 1 2 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2]);
yticklabels({'$10^{-12}$','$10^{-10}$','$10^{-8}$', '$10^{-6}$', '$10^{-4}$', '$10^{-2}$'});

% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15); %#ok<*NASGU>
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.99, 0]);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_20.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_20.eps'),'-depsc2','-r600');

%% m = 1, x/D = 80

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);

first_spod_mode  = eigvalues_spectrum_plot(5).eigvalue(:,2);
second_spod_mode = eigvalues_spectrum_plot(5).eigvalue(:,1);
hold on;
box on;

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  plot(eigvalues_spectrum_plot(5).freq',eigvalues_spectrum_plot(5).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    grid on;
end
patch([f(idx)' fliplr(f(idx)')], [second_spod_mode(idx)'  fliplr(first_spod_mode(idx)')], 'r');
hold off;
set(gca, 'XScale', 'log', 'YScale','log');

ax = gca;
ax.FontSize = 20;

xlim([2.7*10^-2 6]);
xticks([0.1 0.5 1 2 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2]);
set(gca, 'Yticklabel',[]);

% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15); %#ok<*NASGU>
%hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.99, 0]);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_80.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'eigenspectrum_m1_x_D_80.eps'),'-depsc2','-r600');

%% m = 2, x/D = 20

axes(ha(3))

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);

first_spod_mode  = eigvalues_spectrum_plot(3).eigvalue(:,2);
second_spod_mode = eigvalues_spectrum_plot(3).eigvalue(:,1);
hold on;
box on;

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  plot(eigvalues_spectrum_plot(3).freq',eigvalues_spectrum_plot(3).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    grid on;
end
patch([f(idx)' fliplr(f(idx)')], [second_spod_mode(idx)'  fliplr(first_spod_mode(idx)')], 'r');
hold off;
set(gca, 'XScale', 'log', 'YScale','log');

ax = gca;
ax.FontSize = 20;

xlim([2.7*10^-2 6]);
xticks([0.1 0.5 1 2 6]);
xticklabels({'0.1','0.5','1', '2', '6'});
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2]);
yticklabels({'$10^{-12}$','$10^{-10}$','$10^{-8}$', '$10^{-6}$', '$10^{-4}$', '$10^{-2}$'});

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.99, 0]);


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_20.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'eigenspectrum_m2_x_D_20.eps'),'-depsc2','-r600');

%% m = 2, x/D = 80

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

Nplot = 25;   % No. of modes to plot
C = repmat(linspace(1,0.1,Nplot).',1,3);

first_spod_mode  = eigvalues_spectrum_plot(6).eigvalue(:,2);
second_spod_mode = eigvalues_spectrum_plot(6).eigvalue(:,1);
hold on;
box on;

for i = 1:Nplot
    grey = C(Nplot-i+1,:);
    h1(i) =  plot(eigvalues_spectrum_plot(6).freq',eigvalues_spectrum_plot(6).eigvalue(:,i), 'LineWidth',2,'Color',grey);
    grid on;
end
patch([f(idx)' fliplr(f(idx)')], [second_spod_mode(idx)'  fliplr(first_spod_mode(idx)')], 'r');
hold off;
set(gca, 'XScale', 'log', 'YScale','log');

ax = gca;
ax.FontSize = 20;
xlim([2.710^-2 6]);
xticks([0.1 0.5 1 2 6]);
xticklabels({'0.1','0.5','1', '2', '6'});
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2]);
set(gca, 'Yticklabel',[]);

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
% hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.99, 0]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m12_x_D_20_80.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m12_x_D_20_80.eps'),'-depsc2','-r600');
%% Figure 7: Contour plot of energy distribution

clearvars -except dirout;
filename = './files/eigvalue_contourf_diff_loc.mat';
load(filename);

close all;
x0=0;
y0=0;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,2,[.05, 0.1],[.1,.05],[.1 .1]);

%% x/D = 20, Contourmap

axes(ha(1));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(1).FREQ(1:24,:), eigvalues_contourf_plot(1).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(1).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

ax = gca;
ax.FontSize = 20; 

xlim([0 5])
xticks([0 1 2 3 4 5]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
ylim([0 0.5]);
yticks([0 .1 .2 .3 .4 .5]);
yticklabels({'0','0.1','0.2', '0.3', '0.4', '0.5'});

colormap jet;
hBar = colorbar;

set(hBar, 'YTick', 1.5*10^-4:1.5*10^-4:9*10^-4);
set(hBar,'Position',[ax.Position(1)+0.36 ax.Position(2) 0.01 0.40])% To change size
set(hBar, 'TickLabelInterpreter', 'latex');

box on;

% hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.2, 0.95, 0]);

% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',20),'_spod.png'),'-dpng','-r600');
% print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',20),'_spod.eps'),'-depsc','-r600');

%% x/D = 40, Contourmap

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(2).FREQ(1:24,:), eigvalues_contourf_plot(2).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(2).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

ax = gca;
ax.FontSize = 20; 

xlim([0 5])
xticks([0 1 2 3 4 5]);
set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
ylim([0 0.5]);
yticks([0 .1 .2 .3 .4 .5]);
set(gca,'Yticklabel',[]);

colormap jet;
hBar = colorbar;
set(hBar,'Position',[ax.Position(1)+0.36 ax.Position(2) 0.01 0.40])% To change size
set(hBar, 'YTick', 1*10^-4:1*10^-4:6*10^-4);
set(hBar, 'TickLabelInterpreter', 'latex');

box on;

% hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.06, 0.95, 0]);

% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',40),'_spod.png'),'-dpng','-r600');
% print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',40),'_spod.eps'),'-depsc','-r600');

%% x/D = 80, Contourmap

axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(3).FREQ(1:24,:), eigvalues_contourf_plot(3).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(3).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

ax = gca;
ax.FontSize = 20; 
xlim([0 5])
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});
ylim([0 0.5]);
yticks([0 .1 .2 .3 .4 .5]);
yticklabels({'0','0.1','0.2', '0.3', '0.4', '0.5'});
        
colormap jet;
hBar = colorbar;
set(hBar,'Position',[ax.Position(1)+0.36 ax.Position(2) 0.01 0.40])% To change size
set(hBar, 'YTick', 0.6*10^-4:0.6*10^-4:3.6*10^-4);
set(hBar, 'TickLabelInterpreter', 'latex');

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
hYLabel = ylabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.2, 0.95, 0]);

% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'contourf_allmf_x_D_',sprintf('%03d',60),'_spod.png'),'-dpng','-r600');
% print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_',sprintf('%03d',60),'_spod.eps'),'-depsc','-r600');

%% x/D = 100, Contourmap

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(4).FREQ(1:24,:), eigvalues_contourf_plot(4).MODE(1:24,:), (squeeze(eigvalues_contourf_plot(4).eigvalue(1:24,1,1:6))), ...
            'Linestyle', 'none'); %#ok<*NASGU>

ax = gca;
ax.FontSize = 20; 
xlim([0 5])
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});
ylim([0 0.5]);
yticks([0 .1 .2 .3 .4 .5]);
set(gca,'Yticklabel',[]);
       
        
colormap jet;
hBar = colorbar;
set(hBar, 'YTick', 0.5*10^-4:0.5*10^-4:3*10^-4);
set(hBar,'Position',[ax.Position(1)+0.36 ax.Position(2) 0.01 0.40])% To change size
set(hBar, 'TickLabelInterpreter', 'latex');

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
% hYLabel = ylabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.06, 0.95, 0]);

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_spod','.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_spod','.eps'),'-depsc','-r600');
%% Decay of leading SPOD eigenvalues as a function of x/D for m=1 and m=2

clearvars -except dirout;
close all;
filename = './files/eigvalue_decay_diff_loc.mat';
load(filename);

figure;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(x_sample_closest(1:end-1), squeeze(eigenspectra_allm(6,1,2,1:end-1)), 'ko');
hold on;
h2 = loglog(x_sample_closest(1:end-1), squeeze(eigenspectra_allm(1,1,3,1:end-1)), 'ks');
y1 = 0.035*x_sample_closest(2:11,1).^(coeffs_m1_st0136_55(1));
h3 = loglog(x_sample_closest(2:11,1),y1, 'k--', 'Linewidth', 1.5);
y2 = 0.12*x_sample_closest(14:21,1).^(coeffs_m1_st0136_70_120(1));
h4 = loglog(x_sample_closest(14:21,1),y2, 'r--', 'Linewidth', 1.5);
y3 = 0.0054*x_sample_closest(2:21,1).^(coeffs_m2_st0(1));
h5 = loglog(x_sample_closest(2:21,1),y3, 'b--', 'Linewidth', 1.5);
%h6 = text(12, 0.003,'$x^{-1.13}$','interpreter','latex','FontSize', 12);
%h7 = text(55, 0.00015,'$x^{-1.42}$','interpreter','latex','FontSize', 12);
%h8 = text(18, 0.00055,'$x^{-0.60}$','interpreter','latex','FontSize', 12);

xlim([1 125]);
xticks([1 10 100]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\lambda^{(1)}$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4,h5], '$m=1, \mbox{\textit{St}}=0.136$ (VS)', '$m=2, \mbox{\textit{St}}=0$ (DH)', '$x^{-1.13}$','$x^{-1.45}$','$x^{-0.65}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southwest';

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'eigvaluedecay_x_D.eps'),'-depsc2','-r600');

%% Similarity analysis of the m=2, St=0 and m = 1, St = 0.135 leading order SPOD mode
close all;
figure;
x0=5;
y0=5;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.15 .02]);

%% Similarity of m = 2, St = 0

axes(ha(1));
filename = './files/eigvalues_similarity_diff_loc.mat';
load(filename);

C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

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


count = 1;
Nf_sampled = 100;
hold on;

for i = 4:2:20
   disp(i);
%    plot(f(1:Nf_sampled)*(LK_TKE_loc_planes(i,2)/TKE_centerline_loc_planes(i,2)^0.5), eigenspectra_allm(:,1,3,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count});
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,3,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count});

   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0  0.5])
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'0','0.1', '0.2', '0.3', '0.4', '0.5'});
ylim([0 0.2]);
yticks([0 0.05 0.1 0.15 0.20]);
yticklabels({'0','0.05', '0.10', '0.15', '0.20'});

box on;

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\lambda^{(1)}(m=2, \mbox{\textit{St}}=0)/(K_{o}^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.16, 0.98, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.eps'),'-depsc2','-r600');

%% Similarity analysis of the m=1, St=0.136 leading order SPOD mode

axes(ha(2));

C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

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


count = 1;
Nf_sampled = 100;
hold on;
for i = 4:2:20
   disp(i);
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,2,i)/((ud_centerline_loc_planes(i,2))*LK_mean_loc_planes(i,2))^2, C{count});
   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0  0.5])
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'0','0.1', '0.2', '0.3', '0.4', '0.5'});
ylim([0 0.2]);
yticks([0 0.05 0.1 0.15 0.20, 0.25]);
yticklabels({'0','0.05', '0.10', '0.15', '0.20', '0.25'});


box on;

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\lambda^{(1)}(m=1, \mbox{\textit{St}}=0.135)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.44, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.16, 0.98, 0]);

% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\lambda^{(1)}(m=1, \mbox{\textit{St}}=0.135)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_m1_st0136_m2_st0_eigvalue_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_m1_st0136_m2_st0_eigvalue_x_D', '.eps'),'-depsc2','-r600');
%% Unscaled analysis of the m=2, St=0 leading order SPOD mode

% clearvars -except dirout;
% close all;
% filename = './files/eigvalues_similarity_diff_loc.mat';
% load(filename);
% 
% C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle
% 
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
% 
% count = 1;
% Nf_sampled = 100;
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(f(1:Nf_sampled), eigenspectra_allm(:,1,3,i), C{count});
%    count = count + 1;
% end
% 
% 
% 
% box on;
% 
% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$\lambda^{(1)}(m=2, \mbox{\textit{St}}=0)/(K_{o}^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);
% 
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'unscaled_m2_st0_eigvalue_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'unscaled_m2_st0_eigvalue_x_D', '.eps'),'-depsc2','-r600');


%% Similarity analysis of the m=1, St=0.136 leading order SPOD mode

% clearvars -except dirout;
% close all;
% filename = './files/eigvalues_similarity_diff_loc.mat';
% load(filename);
% 
% C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle
% 
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
% 
% count = 1;
% Nf_sampled = 100;
% figure;
% hold on;
% for i = 4:2:20
%    disp(i);
%    plot(f(1:Nf_sampled), eigenspectra_allm(:,1,2,i), C{count});
%    count = count + 1;
% end
% 
% xlim([0 0.5]);
% ylim([0 0.002]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\lambda^{(1)}(m=1, \mbox{\textit{St}}=0.136)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',15);
% 
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'unscaled_m1_st0136_eigvalue_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'unscaled_m1_st0136_eigvalue_x_D', '.eps'),'-depsc2','-r600');

%% Bar chart plot for the dominance shift from m = 2 to m = 1 mode 

clearvars -except dirout;
close all;
filename = './files/eigvalue_barplot_diff_loc.mat';
load(filename);

close all;
x0=5;
y0=5;
width=15;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,2,[.05, 0.05],[.1,.05],[.1 .1]);
%% x/D = 20, Bar plot

axes(ha(1));

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
ax.FontSize = 20; 

% xlim([0  10])
xticks([0 1 2 3 4 5 6 7 8 9 10]);
set(gca,'Xticklabel',[]);
ylim([0 20]);
yticks([0 5 10 15 20]);
yticklabels({'0','5', '10', '15', '20'});


labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';

hYLabel = ylabel('$\xi^{(1)}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.15, 0.95, 0]);


% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',20),'_spod.png'),'-dpng2','-r600');
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',20),'_spod.eps'),'-depsc2','-r600');

%% x/D = 40, bar plot

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(2).mode, eigvalues_bar_plot(2).eigenvalue_fraction_contri(1:3,:)','grouped');

ax = gca;
ax.FontSize = 20; 

% xlim([0  10])
xticks([0 1 2 3 4 5 6 7 8 9 10]);
set(gca,'Xticklabel',[]);
ylim([0 20]);
yticks([0 5 10 15 20]);
set(gca,'Yticklabel',[]);


% labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
% hLegend = legend(labels,'Location','NorthEast');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 20;
% hLegend.FontWeight = 'bold';

% hXLabel = xlabel('$m$','interpreter','latex','fontsize',15);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',40),'_spod.png'),'-dpng2','-r600');
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',40),'_spod.eps'),'-depsc2','-r600');

%% x/D = 80, bar plot

axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(3).mode, eigvalues_bar_plot(3).eigenvalue_fraction_contri(1:3,:)','grouped');

ax = gca;
ax.FontSize = 20; 

% xlim([0  10])
xticks([0 1 2 3 4 5 6 7 8 9 10]);
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10'});
ylim([0 20]);
yticks([0 5 10 15 20]);
yticklabels({'0','5', '10', '15', '20'});


% labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
% hLegend = legend(labels,'Location','NorthEast');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 20;
% hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
hYLabel = ylabel('$\xi^{(1)}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.15, 0.95, 0]);

% set(gcf, 'PaperPositionMode', 'auto');  
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',80),'_spod.png'),'-dpng2','-r600');
% print(gcf,strcat(dirout, 'bar_allm_x_D_',sprintf('%03d',80),'_spod.eps'),'-depsc2','-r600');

%% x/D = 100, bar plot

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = bar(eigvalues_bar_plot(4).mode, eigvalues_bar_plot(4).eigenvalue_fraction_contri(1:3,:)','grouped');

ylim([0 20]);

% set(h1(1),'FaceColor', 'c');
% set(h1(2),'FaceColor', 'b');
% set(h1(3),'FaceColor', 'm');

ax = gca;
ax.FontSize = 20; 

% xlim([0  10])
xticks([0 1 2 3 4 5 6 7 8 9 10]);
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10'});
ylim([0 20]);
yticks([0 5 10 15 20]);
set(gca,'Yticklabel',[]);


% labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
% hLegend = legend(labels,'Location','NorthEast');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
% hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);


set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_204080100_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_204080100_spod.eps'),'-depsc2','-r600');
