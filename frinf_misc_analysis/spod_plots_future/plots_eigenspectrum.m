%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting the eigenspectrum

clear;
addpath('./aux_plots/');
dirout = './';
%% Figure of contour plot of energy distribution

filename = './aux_plots/files/eigvalue_spectrum_diff_loc.mat';
load(filename);

%% Setting up the figure environment
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
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.2, 0.45, 0]);
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

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.2, 0.45, 0]);
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
xlim([2.7*10^-2 6]);
xticks([0.1 0.5 1 2 6]);
xticklabels({'0.1','0.5','1', '2', '6'});
ylim([10^-12 1*10^-1]);
yticks([10^-12 10^-10 10^-8 10^-6 10^-4 10^-2]);
set(gca, 'Yticklabel',[]);

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
% hYLabel = ylabel('$\lambda$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.99, 0]);

%% Saving the plot
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'eigenspectrum_m12_x_D_20_80.png'),'-dpng','-r600');  
print(gcf,strcat(dirout, 'eigenspectrum_m12_x_D_20_80.eps'),'-depsc','-r600');

