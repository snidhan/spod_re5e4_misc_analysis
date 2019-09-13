%% Name - Sheel Nidhan
%  Date - 10th September, 2019
%  Plots for the statistical results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';


%% Plotting the defect velocity as a function of x/D

close;
filename = './files/Defect_centerline.dat';

udefect  = importdata(filename);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(udefect(:,1), udefect(:,2), 'ko');
hold on;
y1 = 0.5*udefect(10:25,1).^(-0.9);
h2 = loglog(udefect(10:25,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.25*udefect(30:43,1).^(-0.6);
h3 = loglog(udefect(30:43,1),y2, 'r-', 'Linewidth', 2);
h4=text(10, 0.03,'$x^{-0.9}$','interpreter','latex','FontSize', 12);
h5=text(70, 0.025,'$x^{-0.6}$','interpreter','latex','FontSize', 12);


xlim([0 125]);
xticks([1 10 100]);

ax = gca;
ax.FontSize = 16; 



hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$U_{d}$','interpreter','latex','fontsize',15);


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'udefect_x_D.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'udefect_x_D.eps'),'-depsc2','-r600');

%% Plotting the wake width as a function of x/D

close;

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
filename = './files/Half_length_zwhazi_WMEAN.dat';

wake_width  = importdata(filename);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(wake_width(:,1), wake_width(:,4), 'ko');
hold on;
y1 = 0.35*wake_width(10:25,1).^(0.45);
h2 = loglog(wake_width(10:25,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.68*wake_width(30:43,1).^(0.3);
h3 = loglog(wake_width(30:43,1),y2, 'r-', 'Linewidth', 2);
h4=text(20, 1.3,'$x^{0.45}$','interpreter','latex','FontSize', 12);
h5=text(75, 2.4,'$x^{0.3}$','interpreter','latex','FontSize', 12);

xlim([0 125]);
xticks([1 10 100]);
ylim([0.5 4]);
yticks([0.5 1 2  3 4]);

ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$L_{d}$','interpreter','latex','fontsize',15);


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ldefect_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ldefect_x_D.eps'),'-depsc2','-r600');

%% Plotting centerline turbulent quantities as a function of x/D

close;

filename = './files/UZ_rms_centerline.dat';
uz_rms   = importdata(filename);

filename = './files/UX_rms_centerline.dat';
ux_rms   = importdata(filename);

filename = './files/UY_rms_centerline.dat';
uy_rms   = importdata(filename);

filename = './files/TKE_centerline.dat';
tke_centerline   = importdata(filename);


set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


h1 = loglog(ux_rms(:,1), ux_rms(:,2), 'bs');
hold on;
h2 = loglog(uy_rms(:,1), uy_rms(:,2), 'rd');
h3 = loglog(uz_rms(:,1), uz_rms(:,2), 'gv');
h4 = loglog(tke_centerline(:,1), tke_centerline(:,2).^0.5, 'ko');
y1 = 0.3*wake_width(10:25,1).^(-2/3);
h5 = loglog(wake_width(10:25,1),y1, 'k-', 'Linewidth', 2);
h6 = text(15, 0.030,'$x^{-2/3}$','interpreter','latex','FontSize', 12);

xlim([0 125]);
xticks([1 10 100]);
ax = gca;
ax.FontSize = 16; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$K^{1/2}$, $u_{x}^{''}$, $u_{y}^{''}$, $u_{z}^{''}$','interpreter','latex','fontsize',15);


hLegend = legend([h1,h2,h3,h4], '$u_{x}^{''}$', '$u_{y}^{''}$', '$u_{z}^{''}$', '$K^{1/2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'turbcenterline_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'turbcenterline_x_D.eps'),'-depsc2','-r600');

%% Plotting the TKE wake width as a function of x/D

close;

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
filename = './files/Half_length_zwhazi_TKE.dat';

tke_wake_width  = importdata(filename);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

h1 = loglog(tke_wake_width(:,1), tke_wake_width(:,4), 'ko');
hold on;
y1 = 0.68*wake_width(15:40,1).^(1/3);
h2 = loglog(wake_width(15:40,1),y1, 'k-', 'Linewidth', 2);
h3=text(50, 2.4,'$x^{1/3}$','interpreter','latex','FontSize', 12);

xlim([0 125]);
xticks([1 10 100]);
ylim([0.5 4]);
yticks([0.5 1 2  3 4]);

ax = gca;
ax.FontSize = 16; 


hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15); %#ok<*NASGU>
hYLabel = ylabel('$L_{k}$','interpreter','latex','fontsize',15);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ltke_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ltke_x_D.eps'),'-depsc2','-r600');

%% Plotting the ratio of K^1/2 and Udefect as a function of x/D

close;

filename = './files/UZ_rms_centerline.dat';
uz_rms   = importdata(filename);

filename = './files/UX_rms_centerline.dat';
ux_rms   = importdata(filename);

filename = './files/UY_rms_centerline.dat';
uy_rms   = importdata(filename);

filename = './files/TKE_centerline.dat';
tke_centerline   = importdata(filename);

filename = './files/Defect_centerline.dat';
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
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$K^{1/2}/U_{d}$','interpreter','latex','fontsize',15);


hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('($K^{1/2}$, $u_{x}^{''}$, $u_{y}^{''}$, $u_{z}^{''}$)/$(U_{d})$','interpreter','latex','fontsize',15);


hLegend = legend([h1,h2,h3,h4], '$u_{x}^{''}/U_{d}$', '$u_{y}^{''}/U_{d}$', '$u_{z}^{''}/U_{d}$', '$K^{1/2}/U_{d}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Location = 'southeast';


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'ratioturbud_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ratioturbud_x_D.eps'),'-depsc2','-r600');

%% Similarity profiles of <u_{x}u_{r}> at different x/D locations using defect velocity wake width

close;

load('./files/similarity_uxur.mat');

figure; 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

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
for i = 2:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 3]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/max($<u_{x}u_{r}>_{r}$)','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld', '.eps'),'-depsc2','-r600');

close;

% Unnormalized profiles of <u_{x}u_{r}> at different x/D locations

figure; 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

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
for i = 2:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 3]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/max($<u_{x}u_{r}>_{r}$)','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk', '.eps'),'-depsc2','-r600');

close;

% Similarity profiles of <u_{x}u_{r}> at different x/D locations using tke wake width

figure; 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

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
for i = 2:2:20
   disp(i);
   plot(rc, -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 8]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/max($<u_{x}u_{r}>_{r}$)','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'unnormalized_uxur_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'unnormalized_uxur_x_D', '.eps'),'-depsc2','-r600');

close;



%% Scaling of <uxur> with different quantities

close all;

load('./files/scaling_uxur.mat');

% Scaling with defect velocity
figure; 
hold on;
C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

count = 1;

Legend = cell(8,1);
Legend{1} = 'x/D = 30';
Legend{2} = 'x/D = 40';
Legend{3} = 'x/D = 50';
Legend{4} = 'x/D = 60';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 80';
Legend{7} = 'x/D = 90';
Legend{8} = 'x/D = 100';

count = 1;
for i = 6:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_Ud(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 4]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{d}^{2}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedud_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedud_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with TKE
figure; 
hold on;
C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

count = 1;

Legend = cell(8,1);
Legend{1} = 'x/D = 30';
Legend{2} = 'x/D = 40';
Legend{3} = 'x/D = 50';
Legend{4} = 'x/D = 60';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 80';
Legend{7} = 'x/D = 90';
Legend{8} = 'x/D = 100';

count = 1;
for i = 6:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_tke(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 4]);
ylim([0 0.4]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/$K$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with Kddelta/dx
figure; 
hold on;
C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

count = 1;

Legend = cell(8,1);
Legend{1} = 'x/D = 30';
Legend{2} = 'x/D = 40';
Legend{3} = 'x/D = 50';
Legend{4} = 'x/D = 60';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 80';
Legend{7} = 'x/D = 90';
Legend{8} = 'x/D = 100';

count = 1;
for i = 6:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_tkeddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 4]);
ylim([0 0.6]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}KdL_{k}/dx$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedk_dlkdx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedk_dlkdx_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with udddelta/dx
figure; 
hold on;
C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

count = 1;

Legend = cell(8,1);
Legend{1} = 'x/D = 30';
Legend{2} = 'x/D = 40';
Legend{3} = 'x/D = 50';
Legend{4} = 'x/D = 60';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 80';
Legend{7} = 'x/D = 90';
Legend{8} = 'x/D = 100';

count = 1;
for i = 6:2:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_udddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 4]);
ylim([0 1]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}U_{d}dL_{d}/dx$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedud_dlddx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedud_dlddx_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with ud(lk/ld)dld/dx
figure; 
hold on;
C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

count = 1;

Legend = cell(8,1);
Legend{1} = 'x/D = 30';
Legend{2} = 'x/D = 40';
Legend{3} = 'x/D = 50';
Legend{4} = 'x/D = 60';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 80';
Legend{7} = 'x/D = 90';
Legend{8} = 'x/D = 100';

count = 1;
for i = 6:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_udlk_ld_ddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
   count = count + 1;
end

xlim([0 4]);
ylim([0 0.4]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}U_{d}(L_{k}/L_{d})dL_{d}/dx$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedudlk_ld_dlddx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedudlk_ld_dlddx_uxur_x_D', '.eps'),'-depsc2','-r600');