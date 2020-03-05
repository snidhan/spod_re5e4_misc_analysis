%% Name - Sheel Nidhan
%  Date - 15 December, 2019

% Code for plotting the centerline statistical quantities from re5e4 frinf
% data

clear;
addpath('./aux_plots/');
dirout = './';

%% Figure 1 - Plotting the defect velocity and wake width based on that as a function of x/D

close all;
figure;
x0=0;
y0=0;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);

% Plotting the defect velocity as a function of x/D

filename = './aux_plots/files/Defect_centerline.dat';
udefect  = importdata(filename);

loc = udefect(33:end,1);                         
log_max_uw = log(udefect(33:end,2));
log_loc =  log(loc);
[coeffs, S] = polyfit(log_loc, log_max_uw, 1);
y_fitted = polyval(coeffs, log_loc);

axes(ha(1));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(udefect(:,1), udefect(:,2), 'ko');
hold on;
y1 = 0.5*udefect(10:25,1).^(-0.9);
h2 = loglog(udefect(10:25,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.25*udefect(30:43,1).^(-0.6);
h3 = loglog(udefect(30:43,1),y2, 'r-', 'Linewidth', 2);
h4=text(10, 0.028,'$x^{-0.9}$','interpreter','latex','FontSize', 20);
h5=text(65, 0.025,'$x^{-0.6}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([1 10 100]);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$U_{d}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

ax = gca;
ax.FontSize = 20; 

% Plotting the wake width as a function of x/D

filename = './aux_plots/files/Half_length_zwhazi_WMEAN.dat';
wake_width  = importdata(filename);

loc = wake_width(33:end,1);                         
log_max_uw = log(wake_width(33:end,2));
log_loc =  log(loc);
[coeffs, S] = polyfit(log_loc, log_max_uw, 1);
y_fitted = polyval(coeffs, log_loc);

axes(ha(2));

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(wake_width(:,1), wake_width(:,4), 'ko');
hold on;
y1 = 0.35*wake_width(10:25,1).^(0.45);
h2 = loglog(wake_width(10:25,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.68*wake_width(30:43,1).^(0.3);
h3 = loglog(wake_width(30:43,1),y2, 'r-', 'Linewidth', 2);
h4=text(20, 1.3,'$x^{0.45}$','interpreter','latex','FontSize', 20);
h5=text(75, 2.35,'$x^{0.3}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([1 10 100]);
ylim([0.5 4]);
yticks([0.5 1 2  3 4]);

ax = gca;
ax.FontSize = 20; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$L_{d}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'udefect_ldefect_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'udefect_ldefect_x_D.eps'),'-depsc2','-r600');

%% Figure 2 - Plotting centerline turbulent quantities and L_{k} as a function of x/D

close all;
figure;
x0=0;
y0=0;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);


filename = './aux_plots/files/UZ_rms_centerline.dat';
uz_rms   = importdata(filename);

filename = './aux_plots/files/UX_rms_centerline.dat';
ux_rms   = importdata(filename);

filename = './aux_plots/files/UY_rms_centerline.dat';
uy_rms   = importdata(filename);

filename = './aux_plots/files/TKE_centerline.dat';
tke_centerline   = importdata(filename);


set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

axes(ha(1));

h1 = loglog(ux_rms(:,1), ux_rms(:,2), 'bs');
hold on;
h2 = loglog(uy_rms(:,1), uy_rms(:,2), 'rd');
h3 = loglog(uz_rms(:,1), uz_rms(:,2), 'gv');
h4 = loglog(tke_centerline(:,1), tke_centerline(:,2).^0.5, 'ko');
y1 = 0.3*uy_rms(10:25,1).^(-2/3);
h5 = loglog(uy_rms(10:25,1),y1, 'k-', 'Linewidth', 2);
h6 = text(12, 0.030,'$x^{-2/3}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([1 10 100]);
ax = gca;
ax.FontSize = 20; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K_{o}^{1/2}$','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0] );
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

hLegend = legend([h1,h2,h3,h4], '$\langle u_{x}''u_{x}'' \rangle^{1/2}_{r=0}$', '$\langle u_{y}''u_{y}'' \rangle^{1/2}_{r=0}$', ...
    '$\langle u_{z}''u_{z}'' \rangle^{1/2}_{r=0}$', '$K_{o}^{1/2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
hLegend.Location = 'southwest';

% Plotting the TKE wake width as a function of x/D

axes(ha(2));

filename = './aux_plots/files/Half_length_zwhazi_TKE.dat';

tke_wake_width  = importdata(filename);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(tke_wake_width(:,1), tke_wake_width(:,4), 'ko');
hold on;
y1 = 0.68*tke_wake_width(15:40,1).^(1/3);
h2 = loglog(tke_wake_width(15:40,1),y1, 'k-', 'Linewidth', 2);
h3=text(50, 2.4,'$x^{1/3}$','interpreter','latex','FontSize', 20);

ax = gca;
ax.FontSize = 20; 
xlim([0 125]);
xticks([1 10 100]);
ylim([0.5 4]);
yticks([0.5 1 2  3 4]);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20); %#ok<*NASGU>
hYLabel = ylabel('$L_{k}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'k_ltke_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'k_ltke_x_D.eps'),'-depsc2','-r600');

%% Figure 3 - Plotting the ratio of K^1/2 and Udefect as a function of x/D

close;

filename = './aux_plots/files/UZ_rms_centerline.dat';
uz_rms   = importdata(filename);

filename = './aux_plots/files/UX_rms_centerline.dat';
ux_rms   = importdata(filename);

filename = './aux_plots/files/UY_rms_centerline.dat';
uy_rms   = importdata(filename);

filename = './aux_plots/files/TKE_centerline.dat';
tke_centerline   = importdata(filename);

filename = './aux_plots/files/Defect_centerline.dat';
udefect  = importdata(filename);

ratio_tke = tke_centerline(:,2).^(0.5)./udefect(:,2);
ratio_ux  = ux_rms(:,2)./udefect(:,2);
ratio_uy  = uy_rms(:,2)./udefect(:,2);
ratio_uz  = uz_rms(:,2)./udefect(:,2);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = plot(udefect(:,1), ratio_ux, 'bs');
hold on;
h2 = plot(udefect(:,1), ratio_uy, 'rd');
h3 = plot(udefect(:,1), ratio_uz, 'gv');
h4 = plot(udefect(:,1), ratio_tke, 'ko');
y1 = 1.22*udefect(:,1).^(0);
h5 = plot(udefect(:,1),y1, 'k--', 'Linewidth', 2);
y1 = 1.5*udefect(:,1).^(0);
h6 = plot(udefect(:,1),y1, 'r--', 'Linewidth', 2);

xlim([0 125]);
xticks([0 20 40 60 80 100 120]);
ylim([0 2.1]);
yticks([0 0.3 0.6 0.9 1.2 1.5 1.8 2.1])

ax = gca;
ax.FontSize = 15; 

box on;

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$K_{o}^{1/2}/U_{d}$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4], '$\langle u_{x}^{''}u_{x}^{''} \rangle^{1/2}_{r=0}/U_{d}$', ...
    '$\langle u_{y}^{''}u_{y}^{''} \rangle^{1/2}_{r=0}/U_{d}$', ...
    '$\langle u_{z}^{''}u_{z}^{''} \rangle^{1/2}_{r=0}/U_{d}$', ...
    '$K_{o}^{1/2}/U_{d}$');

hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southeast';

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ratioturbud_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ratioturbud_x_D.eps'),'-depsc2','-r600');
