%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting the contour map

clear;
addpath('./aux_plots/');
dirout = './';
%% Figure of contour plot of energy distribution

filename = './aux_plots/files/eigvalues_similarity_diff_loc.mat';
load(filename);

%% Setting up the figure environment
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

h1 = contourf(eigvalues_contourf_plot(1).FREQ(1:24,:), ...
              eigvalues_contourf_plot(1).MODE(1:24,:), ...
              (squeeze(eigvalues_contourf_plot(1).eigvalue(1:24,1,1:6))), ...
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
%% x/D = 40, Contourmap

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = contourf(eigvalues_contourf_plot(2).FREQ(1:24,:), ...
              eigvalues_contourf_plot(2).MODE(1:24,:), ...
             (squeeze(eigvalues_contourf_plot(2).eigvalue(1:24,1,1:6))), ...
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
set(hBar, 'YTick', 0.8*10^-4:0.8*10^-4:4.8*10^-4);
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

h1 = contourf(eigvalues_contourf_plot(3).FREQ(1:24,:),...
              eigvalues_contourf_plot(3).MODE(1:24,:),...
              (squeeze(eigvalues_contourf_plot(3).eigvalue(1:24,1,1:6))), ...
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
set(hBar, 'YTick', 0.5*10^-4:0.5*10^-4:3*10^-4);
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
set(hBar, 'YTick', 0.4*10^-4:0.4*10^-4:2.4*10^-4);
set(hBar,'Position',[ax.Position(1)+0.36 ax.Position(2) 0.01 0.40])% To change size
set(hBar, 'TickLabelInterpreter', 'latex');

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
% hYLabel = ylabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.06, 0.95, 0]);


%% Save the plot
set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'contourf_allmf_x_D_spod','.png'),'-dpng','-r600');
print(gcf,strcat(dirout, 'contourf_allmf_f_x_D_spod','.eps'),'-depsc','-r600');
