%% Name - Sheel Nidhan
%  Date - 10th September, 2019
%  Plots for the statistical results for prf_spod_re5e4_frinf paper

dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures/';
% dirout = '/home/sheel/Dropbox/research/sheel_papers/prf_spod_re5e4_frinf/template/figures_2.0/';
% dirout = 'C:\Users\snidh\Dropbox\research\sheel_papers\prf_spod_re5e4_frinf\template\figures_2.0\';
%% Figure 1 - Plotting the defect velocity and wake width based on that as a function of x/D

close all;
figure;
x0=5;
y0=5;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);

% Plotting the defect velocity as a function of x/D

filename = './files/Defect_centerline.dat';
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
h4=text(10, 0.029,'$x^{-0.9}$','interpreter','latex','FontSize', 20);
h5=text(65, 0.025,'$x^{-0.6}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([1 10 100]);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$U_{d}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

ax = gca;
ax.FontSize = 20; 

% Plotting the wake width as a function of x/D

filename = './files/Half_length_zwhazi_WMEAN.dat';
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
x0=5;
y0=5;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);


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

axes(ha(1));

h1 = loglog(ux_rms(:,1), ux_rms(:,2), 'bs');
hold on;
h2 = loglog(uy_rms(:,1), uy_rms(:,2), 'rd');
h3 = loglog(uz_rms(:,1), uz_rms(:,2), 'gv');
h4 = loglog(tke_centerline(:,1), tke_centerline(:,2).^0.5, 'ko');
y1 = 0.3*uy_rms(10:25,1).^(-2/3);
h5 = loglog(uy_rms(10:25,1),y1, 'k-', 'Linewidth', 2);
h6 = text(15, 0.030,'$x^{-2/3}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([1 10 100]);
ax = gca;
ax.FontSize = 20; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K_{o}^{1/2}$','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0] );
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

hLegend = legend([h1,h2,h3,h4], '$\langle u_{x}^{''}u_{x}^{''} \rangle^{1/2}_{r=0}$', '$\langle u_{y}^{''}u_{y}^{''} \rangle^{1/2}_{r=0}$', ...
    '$\langle u_{z}^{''}u_{z}^{''} \rangle^{1/2}_{r=0}$', '$K_{o}^{1/2}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
hLegend.Location = 'southwest';

% Plotting the TKE wake width as a function of x/D

axes(ha(2));

filename = './files/Half_length_zwhazi_TKE.dat';

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

%% Plotting the ratio of L_{k} and L_{d} as a function of x/D

% close;
% 
% filename = './files/Half_length_zwhazi_WMEAN.dat';
% LK_mean  = importdata(filename);
% 
% filename = './files/Half_length_zwhazi_TKE.dat';
% LK_tke   = importdata(filename);
% 
% ratio = LK_tke(:,2)./LK_mean(:,2);

% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% h1 = plot(LK_mean(:,1), ratio, 'ko');
% 
% xlim([0 125]);
% xticks([0 20 40 60 80 100 120]);
% ylim([1 1.5]);
% yticks([1 1.1 1.2 1.3 1.4 1.5])
% 
% ax = gca;
% ax.FontSize = 16; 

% box on;
% 
% hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$K^{1/2}/U_{d}$','interpreter','latex','fontsize',15);
% 
% 
% hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$L_{k}/L_{d}$','interpreter','latex','fontsize',15);


% hLegend = legend([h1,h2,h3,h4], '$u_{x}^{''}/U_{d}$', '$u_{y}^{''}/U_{d}$', '$u_{z}^{''}/U_{d}$', '$K^{1/2}/U_{d}$');
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Location = 'southeast';


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ratiolkld_x_D.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ratiolkld_x_D.eps'),'-depsc2','-r600');


%% Similarity profiles of U at different x/D locations

close;

load('./files/similarity_w.mat');
figure;

hold on;

lineStyles = distinguishable_colors(11);
count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 50';
Legend{2} = 'x/D = 55';
Legend{3} = 'x/D = 60';
Legend{4} = 'x/D = 65';
Legend{5} = 'x/D = 70';
Legend{6} = 'x/D = 75';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 85';
Legend{9} = 'x/D = 90';
Legend{10} = 'x/D = 95';
Legend{11} = 'x/D = 100';
% Legend{10} = 'x/D = 55';
% Legend{11} = 'x/D = 60';
% Legend{12} = 'x/D = 65';
% Legend{13} = 'x/D = 70';
% Legend{14} = 'x/D = 75';
% Legend{15} = 'x/D = 80';
% Legend{16} = 'x/D = 85';
% Legend{17} = 'x/D = 90';
% Legend{18} = 'x/D = 95';
% Legend{19} = 'x/D = 100';
% Legend{11} = 'x/D = 110';
% Legend{11} = 'x/D = 120';

count = 1;
for i = 10:1:20
   disp(i);
   plot(rc(1:nr)/LK_mean_loc_planes(i,2), (1-mean_w_1d(1:nr,i))/max(1-mean_w_1d(:,i)), 'Color', lineStyles(count,:), 'Linewidth',2);

   count = count + 1;
end

xlim([0 3]);
ylim([0 1.2]);

ax = gca;
ax.FontSize = 16; 

box on;

hXLabel = xlabel('$\eta = r/L_{d}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$f(\eta)$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 10;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_w_x_D_ld_50_5_100', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_w_x_D_ld_50_5_100', '.eps'),'-depsc2','-r600');

% Unnormalized profiles of u at different x/D locations

% figure; 
% hold on;
% 
% lineStyles = distinguishable_colors(10);
% count = 1;
% 
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 25';
% Legend{3} = 'x/D = 30';
% Legend{4} = 'x/D = 35';
% Legend{5} = 'x/D = 40';
% Legend{6} = 'x/D = 45';
% Legend{7} = 'x/D = 50';
% Legend{8} = 'x/D = 55';
% Legend{9} = 'x/D = 60';
% 
% 
% count = 1;
% for i = 2:2:12
%    disp(i);
%    plot(rc(1:nr)/LK_mean_loc_planes(i,2), (1-mean_w_1d(1:nr,i))/max(1-mean_w_1d(:,i)), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 10]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$f(\eta)$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'unnormalized_w_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'unnormalized_w_x_D', '.eps'),'-depsc2','-r600');


%% Figure 5: Similarity profiles of <u_{x}u_{r}> at different x/D locations using defect velocity wake width       

close all;

load('./files/similarity_uxur.mat');
reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess',2);
reystress_uw_1d = reystress_uw_1d_smooth;

x0=0;
y0=0;
width=15;
height=15;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);

[ha, pos] = tight_subplot(3,2,[.05, 0.05],[.1,.05],[.1 .02]);
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(ha(1),'Position',[new,pos{1}(2:end)])
delete(ha(2));
%% Plotting <u_x u_r> as a function of r/D from x/D = 20 to 120

axes(ha(1));

%C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
%color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(11,@srgb_to_Jab);
count = 1;

Legend = cell(11,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 30$';
Legend{3} = '$x/D = 40$';
Legend{4} = '$x/D = 50$';
Legend{5} = '$x/D = 60$';
Legend{6} = '$x/D = 70$';
Legend{7} = '$x/D = 80$';
Legend{8} = '$x/D = 90$';
Legend{9} = '$x/D = 100$';
Legend{10} = '$x/D = 110$';
Legend{11} = '$x/D = 120$';

hold on;
for i =  4:2:20
   disp(i);
   plot(rc, -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc, -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 15;

xlim([0 10]);
xticks([0 2 4 6 8 10]);
xticklabels({'0','2','4', '6', '8', '10'});

ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%% Plotting <u_xu_r> as a function of r/Ld from x/D = 20 to 60

axes(ha(3));

hold on;
C = {'-','-','-','-','-','-','-','-','-','-'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

%lineStyles = linspecer(10);
lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 25$';
Legend{3} = '$x/D = 30$';
Legend{4} = '$x/D = 35$';
Legend{5} = '$x/D = 40$';
Legend{6} = '$x/D = 45$';
Legend{7} = '$x/D = 50$';
Legend{8} = '$x/D = 55$';
Legend{9} = '$x/D = 60$';

count = 1;

for i = 4:12
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

% for i =  21:22
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end

ax = gca;
ax.FontSize = 15;

xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca, 'Xticklabel', []);
% xticks([0 1 2 3 4]);
% xticklabels({'0','1','2', '3', '4'});
ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

%hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$-\langle u_{x}u_{r}\rangle$/max$(-\langle u_{x}u_{r}\rangle)_{r}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_20_60', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_20_60', '.eps'),'-depsc2','-r600');


%% Plotting <u_xu_r> as a function of r/Lk from x/D = 20 to 60

axes(ha(4));

hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 25$';
Legend{3} = '$x/D = 30$';
Legend{4} = '$x/D = 35$';
Legend{5} = '$x/D = 40$';
Legend{6} = '$x/D = 45$';
Legend{7} = '$x/D = 50$';
Legend{8} = '$x/D = 55$';
Legend{9} = '$x/D = 60$';

count = 1;
for i = 4:12
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

% for i =  21:22
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end

xlim([0 4]);
xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca, 'Xticklabel', []);
ylim([0 1.2]);

% ax = gca;
% ax.FontSize = 16; 

box on;

% hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$-\langle u_{x}u_{r}\rangle$/max$(-\langle u_{x}u_{r}\rangle)_{r}$','interpreter','latex','fontsize',15);
hTitle  = title('(c)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk_20_60', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk_20_60', '.eps'),'-depsc2','-r600');

%% Plot of <u_{x}u_{r}> using r/Ld from x/D = 70 to 120

axes(ha(5));

hold on;
C = {'-','-','-','-','-','-','-','-','-','-'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

%lineStyles = linspecer(10);
lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 70$';
Legend{2} = '$x/D = 75$';
Legend{3} = '$x/D = 80$';
Legend{4} = '$x/D = 85$';
Legend{5} = '$x/D = 90$';
Legend{6} = '$x/D = 95$';
Legend{7} = '$x/D = 100$';
Legend{8} = '$x/D = 110$';
Legend{9} = '$x/D = 120$';

count = 1;

for i = 14:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end


ax = gca;
ax.FontSize = 15;

xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4'});

ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_70_120', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_70_120', '.eps'),'-depsc2','-r600');

%% Plot of <u_{x}u_{r}> using r/Lk from x/D = 70 to 120

axes(ha(6));
 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(9,@srgb_to_Jab);

count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 70$';
Legend{2} = '$x/D = 75$';
Legend{3} = '$x/D = 80$';
Legend{4} = '$x/D = 85$';
Legend{5} = '$x/D = 90$';
Legend{6} = '$x/D = 95$';
Legend{7} = '$x/D = 100$';
Legend{8} = '$x/D = 110$';
Legend{9} = '$x/D = 120$';

count = 1;
for i = 14:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 15;
xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4'});
ylim([0 1.2]);
% yticks([0 .2 .4 .6 .8 1]);
% yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
% hYLabel = ylabel('$-\langle u_{x}^{''}u_{r}^{''}\rangle$/max$(-\langle u_{x}^{''}u_{r}^{''}\rangle)_{r}$','interpreter','latex','fontsize',10,...
%     'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(e)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_complete', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_complete', '.eps'),'-depsc2','-r600');

%% Figure 4 - Similarity profiles of <K> at different x/D locations using defect velocity wake width

close all;
load('./files/similarity_tke.mat');

% figure; 
% hold on;
% C = {'-','-','-','-','-','-','-','-','-','-'}; % Cell array of linestyle
% color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle
% 
% %lineStyles = linspecer(10);
% lineStyles = maxdistcolor(11,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(11,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% Legend{10} = 'x/D = 110';
% Legend{11} = 'x/D = 120';
% 
% count = 1;
% for i = 4:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% for i =  21:22
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0.05 3]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$K$/$K_{o}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_k_x_D_ld', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_k_x_D_ld', '.eps'),'-depsc2','-r600');


% K profiles normalized by L_{k} at different x/D locations

close all;
figure;
x0=5;
y0=5;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);

axes(ha(2));
hold on;

C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(11,@srgb_to_Jab);
count = 1;

% Legend = cell(11,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 30';
% Legend{3} = 'x/D = 40';
% Legend{4} = 'x/D = 50';
% Legend{5} = 'x/D = 60';
% Legend{6} = 'x/D = 70';
% Legend{7} = 'x/D = 80';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 100';
% Legend{10} = 'x/D = 110';
% Legend{11} = 'x/D = 120';

count = 1;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 16; 

ylim([0 1.2]);
xlim([0.05 5]);
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});


box on;

hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',20);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_k_x_D_lk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_k_x_D_lk', '.eps'),'-depsc2','-r600');


% Similarity profiles of K at different x/D locations using tke wake width

axes(ha(1)); 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(11,@srgb_to_Jab);
count = 1;

Legend = cell(11,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 30$';
Legend{3} = '$x/D = 40$';
Legend{4} = '$x/D = 50$';
Legend{5} = '$x/D = 60$';
Legend{6} = '$x/D = 70$';
Legend{7} = '$x/D = 80$';
Legend{8} = '$x/D = 90$';
Legend{9} = '$x/D = 100$';
Legend{10} = '$x/D = 110$';
Legend{11} = '$x/D = 120$';


for i = 4:2:20
   disp(i);
   plot(rc, tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc, tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0.05 10]);
xticks([0 2 4 6 8 10]);
xticklabels({'0','2','4', '6', '8', '10'});

ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K$/$K_{o}$','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Units = 'normalized';
hLegend.Position = [-0.1 0.1 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'k_func_r_D_r_Lk_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'k_func_r_D_r_Lk_x_D', '.eps'),'-depsc2','-r600');

%% Figure 6 ??: Scaling of max(-<uxur>) as a function of x/D

close all;

load('./files/similarity_uxur.mat');
reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');
reystress_uw_1d = reystress_uw_1d_smooth;

max_uw = zeros(size(reystress_uw_1d,2),1);
for i = 1:size(reystress_uw_1d,2)
    max_uw(i,1) = max(abs(reystress_uw_1d(:,i)));
end

loc = LK_mean_loc_planes(4:22,1);                         
log_max_uw = log(max_uw(4:22,1));
log_loc =  log(loc);
[coeffs, S] = polyfit(log_loc, log_max_uw, 1);
y_fitted = polyval(coeffs, log_loc);

h1 = loglog(LK_mean_loc_planes(:,1), max_uw, 'ko');
hold on;
y1 = 0.15*LK_mean_loc_planes(6:12,1).^(coeffs(1));
h5 = loglog(LK_mean_loc_planes(6:12,1),y1, 'k-', 'Linewidth', 2);
h6 = text(45, 0.0008,'$x^{-1.52}$','interpreter','latex','FontSize', 20);

xlim([0 125]);
xticks([0 1 10 100]);

ax = gca;
ax.FontSize = 20; 

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('max$(-\langle u''_{x}u''_{r} \rangle)_{r}$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'scaling_uxur_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'scaling_uxur_x_D', '.eps'),'-depsc2','-r600');

%% Scaling of spread of <uxur> as a function of x/D

close all;

load('./files/similarity_uxur.mat');
reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');
% reystress_uw_1d_smooth = reystress_uw_1d;

figure;
hold on;
for i = 1:12
    %plot(rc, -reystress_uw_1d(:,i), 'k', 'Linewidth',2);
    plot(rc, -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i))), 'Linewidth',2);
end

ylim([0 1.2]);
%%
reystress_uw_1d = reystress_uw_1d_smooth;
final_increasing_region = zeros(size(reystress_uw_1d,2),1);
idx_final_increasing_region = zeros(size(reystress_uw_1d,2),1);

figure;
hold on;
for i = 1:22
   rs = -reystress_uw_1d(:,i)./max(-reystress_uw_1d(:,i));
   count = 1;
   for j = 2:size(rs,1)-1
    disp(j);
    rs_forward = rs(j+1,1) - rs(j,1);
    rs_backward = rs(j,1) - rs(j-1,1);
    if (rs_forward < 0) || (abs(rs_forward/rs_backward) < 10^-4)
        final_increasing_region(i,1) = rc(j,1);
        idx_final_increasing_region(i,1) = j;
        break;
    end
   end
   
   plot(rc, rs, '-', 'Linewidth',2)
   plot(rc(idx_final_increasing_region(i,1)), rs(idx_final_increasing_region(i,1)), 'k*');
end

%xlim([0 6]);
ylim([-0.01 1.2]);

final_increasing_region(16) = [];
loc = LK_mean_loc_planes(:,1); 
loc(16) = [];
log_max = log(final_increasing_region(4:end,1));
log_loc =  log(loc(4:end,1));
[coeffs, S] = polyfit(log_loc, log_max, 1);
y_fitted = polyval(coeffs, log_loc);

figure;
h1 = loglog(loc, final_increasing_region, 'ko');
hold on;
y1 = 0.44*loc(4:18,1).^(coeffs(1));
h5 = loglog(loc(4:18,1),y1, 'k-', 'Linewidth', 2);
h6 = text(35, 2.1,'$x^{0.38}$','interpreter','latex','FontSize', 20);

xlim([1 125]);
xticks([0 1 10 100]);

ax = gca;
ax.FontSize = 20; 


hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\delta_{max(-\langle u''_{x}u''_{r} \rangle)_{r}}$','interpreter','latex','fontsize',20);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'width_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'width_uxur_x_D', '.eps'),'-depsc2','-r600');  
%% Find the approximation between the slope at r=0 with the scaling of max(uxur)/delta(max(uxur)

close all;

load('./files/similarity_uxur.mat');
reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');

% Finding the maxima
max_uw = zeros(size(reystress_uw_1d,2),1);
for i = 1:size(reystress_uw_1d,2)
    [max_uw(i,1), max_idx(i,1)] = max(abs(reystress_uw_1d(:,i)));
end

for i = 1:size(reystress_uw_1d,2)
   Slope(i).slope = abs(reystress_uw_1d(3:max_idx(i,1),i) - reystress_uw_1d(1:max_idx(i,1)-2,i))./rc(3:max_idx(i,1) - rc(1:max_idx(i,1)-2))    
   Slope(i).rc = rc(2:max_idx(i,1)-1);
end

%% Finding the location of maxima
final_increasing_region = zeros(size(reystress_uw_1d,2),1);
idx_final_increasing_region = zeros(size(reystress_uw_1d,2),1);
for i = 1:size(reystress_uw_1d,2)
   rs = -reystress_uw_1d(:,i)./max(-reystress_uw_1d(:,i));
   count = 1;
   for j = 2:size(rs,1)-1
    disp(j);
    rs_forward = rs(j+1,1) - rs(j,1);
    rs_backward = rs(j,1) - rs(j-1,1);
    if (rs_forward < 0) || (abs(rs_forward/rs_backward) < 10^-4)
        final_increasing_region(i,1) = rc(j,1);
        idx_final_increasing_region(i,1) = j;
        break;
    end
   end
end

slope_scaling = max_uw./final_increasing_region;

% Finding the slope at r = 0;

slope_r0 = zeros(size(reystress_uw_1d,2),1);
for i = 1:size(reystress_uw_1d,2)
    
    slope_r0(i,1) = mean(abs((reystress_uw_1d(5:10,i) - reystress_uw_1d(3:8,i))./(rc(5:10,1) - rc(3:8,1))));

end


figure;

h1 = loglog(LK_mean_loc_planes(:,1),  slope_scaling, 'ro');
hold on;
h2 = plot(LK_mean_loc_planes(:,1),  slope_r0, 'bo');

% percentage_error = abs((slope_r0-slope_scaling))./slope_r0*100;
%% Scaling of <uxur> with different quantities

% close all;
% 
% load('./files/scaling_uxur.mat');

%% Scaling with defect velocity
% figure; 
% hold on;
% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle
% 
% lineStyles = maxdistcolor(11,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(8,1);
% Legend{1} = 'x/D = 30';
% Legend{2} = 'x/D = 40';
% Legend{3} = 'x/D = 50';
% Legend{4} = 'x/D = 60';
% Legend{5} = 'x/D = 70';
% Legend{6} = 'x/D = 80';
% Legend{7} = 'x/D = 90';
% Legend{8} = 'x/D = 100';
% 
% count = 1;
% for i = 6:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_Ud(1:nr,i), 'Color', color{count}, 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 4]);
% ylim([0 1.2]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{d}^{2}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedud_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedud_uxur_x_D', '.eps'),'-depsc2','-r600');

%% Scaling with TKE
% figure; 
% hold on;

% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 70';
% Legend{2} = 'x/D = 75';
% Legend{3} = 'x/D = 80';
% Legend{4} = 'x/D = 85';
% Legend{5} = 'x/D = 90';
% Legend{6} = 'x/D = 95';
% Legend{7} = 'x/D = 100';
% Legend{8} = 'x/D = 110';
% Legend{9} = 'x/D = 120';
% 
% count = 1;
% for i = 14:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_tke(1:nr,i), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% for i = 21:22
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_tke(1:nr,i), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 4]);
% ylim([0 0.4]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$K_{o}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.eps'),'-depsc2','-r600');

%% Scaling with Udddelta/dx
% figure; 
% hold on;

% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle

% lineStyles = maxdistcolor(9,@srgb_to_Jab);
% count = 1;
% 
% Legend = cell(9,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 25';
% Legend{3} = 'x/D = 30';
% Legend{4} = 'x/D = 35';
% Legend{5} = 'x/D = 40';
% Legend{6} = 'x/D = 45';
% Legend{7} = 'x/D = 50';
% Legend{8} = 'x/D = 55';
% Legend{9} = 'x/D = 60';
% 
% count = 1;
% for i = 4:12
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_udddelta_dx(1:nr,i), 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end
% 
% 
% xlim([0 4]);
% ylim([0 0.8]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$K_{o}$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedtke_uxur_x_D', '.eps'),'-depsc2','-r600');

%% Scaling with Kddelta/dx
% figure; 
% hold on;
% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle
% 
% count = 1;
% 
% Legend = cell(8,1);
% Legend{1} = 'x/D = 30';
% Legend{2} = 'x/D = 40';
% Legend{3} = 'x/D = 50';
% Legend{4} = 'x/D = 60';
% Legend{5} = 'x/D = 70';
% Legend{6} = 'x/D = 80';
% Legend{7} = 'x/D = 90';
% Legend{8} = 'x/D = 100';
% 
% count = 1;
% for i = 6:2:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_tkeddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 4]);
% ylim([0 0.6]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}KdL_{k}/dx$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedk_dlkdx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedk_dlkdx_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with udddelta/dx
% figure; 
% hold on;
% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle
% 
% count = 1;
% 
% Legend = cell(8,1);
% Legend{1} = 'x/D = 30';
% Legend{2} = 'x/D = 40';
% Legend{3} = 'x/D = 50';
% Legend{4} = 'x/D = 60';
% Legend{5} = 'x/D = 70';
% Legend{6} = 'x/D = 80';
% Legend{7} = 'x/D = 90';
% Legend{8} = 'x/D = 100';
% 
% count = 1;
% for i = 6:2:20
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_udddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 4]);
% ylim([0 1]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}U_{d}dL_{d}/dx$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedud_dlddx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedud_dlddx_uxur_x_D', '.eps'),'-depsc2','-r600');

% Scaling with ud(lk/ld)dld/dx
% figure; 
% hold on;
% C = {'-','-','-','-','--','--','--','--'}; % Cell array of linestyle
% color = {'k','b','r','g','g','r','b','k'}; % Cell array of linestyle
% 
% count = 1;
% 
% Legend = cell(8,1);
% Legend{1} = 'x/D = 30';
% Legend{2} = 'x/D = 40';
% Legend{3} = 'x/D = 50';
% Legend{4} = 'x/D = 60';
% Legend{5} = 'x/D = 70';
% Legend{6} = 'x/D = 80';
% Legend{7} = 'x/D = 90';
% Legend{8} = 'x/D = 100';
% 
% count = 1;
% for i = 6:2:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d_udlk_ld_ddelta_dx(1:nr,i), C{count}, 'Color', color{count}, 'Linewidth',2);
%    count = count + 1;
% end
% 
% xlim([0 4]);
% ylim([0 0.4]);
% 
% ax = gca;
% ax.FontSize = 16; 
% 
% box on;
% 
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/$U_{\infty}U_{d}(L_{k}/L_{d})dL_{d}/dx$','interpreter','latex','fontsize',15);
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 15;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'normalizedudlk_ld_dlddx_uxur_x_D', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'normalizedudlk_ld_dlddx_uxur_x_D', '.eps'),'-depsc2','-r600');