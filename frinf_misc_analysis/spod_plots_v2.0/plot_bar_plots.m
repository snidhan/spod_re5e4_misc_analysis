%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting the bar plots

clear;
addpath('./aux_plots/');
dirout = './';

%% Reading the filename

filename = './aux_plots/files/eigvalue_barplot_diff_loc.mat';
load(filename);

%% Setting up the figure environment

close all;
x0=0;
y0=0;
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
xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
set(gca,'Xticklabel',[]);
ylim([0 20]);
yticks([0 5 10 15 20]);
yticklabels({'0','5', '10', '15', '20'});


labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
hLegend = legend(labels,'Location','NorthEast');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';

hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
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
xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
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
xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'});
ylim([0 20]);
yticks([0 5 10 15 20]);
yticklabels({'0','5', '10', '15', '20'});


% labels = {'$\lambda^{(1)}$','$\lambda^{(2)}$','$\lambda^{(3)}$'};
% hLegend = legend(labels,'Location','NorthEast');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 20;
% hLegend.FontWeight = 'bold';

hXLabel = xlabel('$m$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
hYLabel = ylabel('$\xi^{(i)}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
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
xticks([0 1 2 3 4 5 6 7 8 9 10 11]);
xticklabels({'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'});
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

%% Saving the plot
set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'bar_allm_x_D_204080100_spod.png'),'-dpng2','-r600');
print(gcf,strcat(dirout, 'bar_allm_x_D_204080100_spod.eps'),'-depsc2','-r600');