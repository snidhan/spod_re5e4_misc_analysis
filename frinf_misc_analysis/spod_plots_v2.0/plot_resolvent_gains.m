%% Name - Sheel Nidhan
%  Date - 24th February, 2020
%  Plots for the resolvent values results for prf_spod_re5e4_frinf paper

clear;
dirout = './';
% dirout = 'C:\Users\snidh\Dropbox\research\sheel_papers\prf_spod_re5e4_frinf\template\figures_2.0\';
%% Plotting the eigenspectrum as a function of x/D

close all;

x0=0;
y0=0;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,2,[.05, 0.08],[.15,.05],[.15 .05]);

%% m = 1, x/D = 40, kx ~ St = 0.135

axes(ha(1));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m1_kx0.8482_x_D_40.mat';
load(filename);

hold on;
box on;

h1 = semilogy(St, SS(1,:,1).^2, 'ko', 'MarkerSize',10);
hold on;
h2 = semilogy(St, SS(2,:,1).^2, 'rs', 'MarkerSize',10);
h3 = semilogy(St, SS(3,:,1).^2, 'bd', 'MarkerSize',10);

hLegend = legend([h1,h2,h3], '$\sigma_{1}^{2}$', '$\sigma_{2}^{2}$', '$\sigma_{3}^{2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
hLegend.Location = 'northeast';


set(gca, 'XScale', 'linear', 'YScale','log')

ax = gca;
ax.FontSize = 20;

xlim([0 0.5]);
xticks([0 0.1 0.2 0.3 0.4 0.5]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([10^-2 10^8]);
yticks([10^0 10^2 10^4 10^6 10^8]);
yticklabels({'$10^{0}$','$10^{2}$','$10^{4}$', '$10^{6}$', '$10^{8}$'});

% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15); %#ok<*NASGU>
hYLabel = ylabel('$\sigma^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.25, 0.99, 0]);

%% m = 2, x/D = 40, kx ~ St = 0

axes(ha(3));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m2_kx0_x_D_40.mat';
load(filename);

hold on;
box on;

h1 = semilogy(St, SS(1,:,1).^2, 'ko', 'MarkerSize',10);
hold on;
h2 = semilogy(St, SS(2,:,1).^2, 'rs', 'MarkerSize',10);
h3 = semilogy(St, SS(3,:,1).^2, 'bd', 'MarkerSize',10);

set(gca, 'XScale', 'linear', 'YScale','log')

ax = gca;
ax.FontSize = 20;

xlim([0 0.5]);
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'$0$','$0.1$','$0.2$', '$0.3$', '$0.4$', '$0.5$'});
ylim([10^-2 10^16]);
yticks([10^0 10^4 10^8 10^12 10^16]);
yticklabels({'$10^{0}$','$10^{4}$','$10^{8}$', '$10^{12}$', '$10^{16}$'});
   
hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20); %#ok<*NASGU>
hYLabel = ylabel('$\sigma^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.25, 0.99, 0]);

%% m = 1, x/D = 100, kx ~ St = 0.135

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m1_kx0.8482_x_D_100.mat';
load(filename);

hold on;
box on;

h1 = semilogy(St, SS(1,:,1).^2, 'ko', 'MarkerSize',10);
hold on;
h2 = semilogy(St, SS(2,:,1).^2, 'rs', 'MarkerSize',10);
h3 = semilogy(St, SS(3,:,1).^2, 'bd', 'MarkerSize',10);

set(gca, 'XScale', 'linear', 'YScale','log')

ax = gca;
ax.FontSize = 20;

xlim([0 0.5]);
xticks([0 0.1 0.2 0.3 0.4 0.5]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([10^-2 10^8]);
yticks([10^0 10^2 10^4 10^6 10^8]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.
    
% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15); %#ok<*NASGU>
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.1, 0.99, 0]);

%% m = 2, x/D = 40, kx ~ St = 0

axes(ha(4));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

filename = './files/m2_kx0_x_D_100.mat';
load(filename);

hold on;
box on;

h1 = semilogy(St, SS(1,:,1).^2, 'ko', 'MarkerSize',10);
hold on;
h2 = semilogy(St, SS(2,:,1).^2, 'rs', 'MarkerSize',10);
h3 = semilogy(St, SS(3,:,1).^2, 'bd', 'MarkerSize',10);

set(gca, 'XScale', 'linear', 'YScale','log')

ax = gca;
ax.FontSize = 20;

xlim([0 0.5]);
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'$0$','$0.1$','$0.2$', '$0.3$', '$0.4$', '$0.5$'});
ylim([10^-2 10^16]);
yticks([10^0 10^4 10^8 10^12 10^16]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.
    
hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20); %#ok<*NASGU>
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.1, 0.99, 0]);

%% Saving the plots

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'resolvent_spectrum_m12_x_D_40_100.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'resolvent_spectrum_m12_x_D_40_100.eps'),'-depsc2','-r600');