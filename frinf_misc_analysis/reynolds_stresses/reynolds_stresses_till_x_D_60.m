%% Name: Sheel Nidhan
%  To calculate the Reynolds stresses and its similarity in the 'non-equilibrium' region
clear; clc; close all;
%% Parameters

loc_planes = [20; 25; 30; 35; 40; 45; 50; 55];
nstart = 1892600;
nend   = 2613200;
stride  = 100;
dir_in_planes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/reystresses/run_2.0/';
nr = 354;
ntheta = 256;
N = (nend-nstart)/stride + 1;

reystress_uw_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_uu_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_ww_2d = zeros(nr,ntheta,size(loc_planes,1));
reystress_vv_2d = zeros(nr,ntheta,size(loc_planes,1));


reystress_uw_1d = zeros(nr,size(loc_planes,1));
reystress_uu_1d = zeros(nr,size(loc_planes,1));
reystress_ww_1d = zeros(nr,size(loc_planes,1));
reystress_vv_1d = zeros(nr,size(loc_planes,1));
tke_1d          = zeros(nr,size(loc_planes,1));
%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);
LK_TKE_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1))); %#ok<*NCOMMA>
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);
LK_mean_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_centerline.dat';
TKE_Centerline   = importdata(filename);
TKE_centerline_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-TKE_Centerline(:,1)));
    TKE_centerline_loc_planes(i,2) = TKE_Centerline(idx,2);
    TKE_centerline_loc_planes(i,1) = TKE_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/UX_rms_centerline.dat';
ux_Centerline   = importdata(filename);
ux_centerline_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ux_Centerline(:,1)));
    ux_centerline_loc_planes(i,2) = ux_Centerline(idx,2);
    ux_centerline_loc_planes(i,1) = ux_Centerline(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);
ud_centerline_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/MKE_areaI_FINF.dat';
mke_area   = importdata(filename);
mke_area_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-mke_area(:,1)));
    mke_area_loc_planes(i,2) = mke_area(idx,2);
    mke_area_loc_planes(i,1) = mke_area(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_areaI_FINF.dat';
tke_area   = importdata(filename);
tke_area_loc_planes = zeros(size(loc_planes,1),2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-tke_area(:,1)));
    tke_area_loc_planes(i,2) = tke_area(idx,2);
    tke_area_loc_planes(i,1) = tke_area(idx,1);
end

%% Reading the Reynolds stresses
for x_loc_planes = 1:size(loc_planes,1)

    filename = strcat(dir_in_planes, 'reystress_x_D_', int2str(loc_planes(x_loc_planes,1)), '.mat');
    disp(filename);
    load(filename);
    
    reystress_uw_2d(:,:,x_loc_planes) =  reystress_uw_av;
    reystress_ww_2d(:,:,x_loc_planes) =  reystress_ww_av;
    reystress_uu_2d(:,:,x_loc_planes) =  reystress_uu_av;
    reystress_vv_2d(:,:,x_loc_planes) =  reystress_vv_av;

    
    reystress_uw_1d(:,x_loc_planes) =  reystress_uw_av_th;
    reystress_ww_1d(:,x_loc_planes) =  reystress_ww_av_th;
    reystress_uu_1d(:,x_loc_planes) =  reystress_uu_av_th;
    reystress_vv_1d(:,x_loc_planes) =  reystress_vv_av_th;
    tke_1d(:,x_loc_planes)          =  reystress_vv_1d(:,x_loc_planes) + reystress_uu_1d(:,x_loc_planes) + reystress_ww_1d(:,x_loc_planes); 
end

reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');
tke_1d_smooth          = smoothdata(tke_1d,'loess');
%% Fit a function to the Ud profile

log_x_loc = log(loc_planes);
log_ud = log(ud_centerline_loc_planes(:,2));
[coeffs_ud, S] = polyfit(log_x_loc, log_ud, 1)

%% Fit a function to the TKE profile

log_x_loc = log(loc_planes);
log_tke = log(TKE_centerline_loc_planes(:,2));
[coeffs_tke, S] = polyfit(log_x_loc, log_tke, 1)
%% Fit a function to the LK_TKE profile

log_x_loc = log(loc_planes);
log_Lk = log(LK_TKE_loc_planes(:,2));
[coeffs_lk, S] = polyfit(log_x_loc, log_Lk, 1)

%% Fit a function to the LK_mean profile

log_x_loc = log(loc_planes);
log_Ld = log(LK_mean_loc_planes(:,2));
[coeffs_ld, S] = polyfit(log_x_loc, log_Ld, 1)

%% Slope of scaling factor Uddld/dx

log_x_D = log(loc_planes);
log_ld  = log(LK_mean_loc_planes(:,2));
[coeffs_ud_dlddx, S] = polyfit(log_x_D, log_ld, 1)

dld_dx = coeffs_ud_dlddx(1)*LK_mean_loc_planes(:,2)./loc_planes
scaling_factor_ud_dld_dx = ud_centerline_loc_planes(:,2).*dld_dx;

%% Slope of max of RS

max_rs = max(abs(reystress_uw_1d_smooth));
log_max_rs = log(max_rs');
[coeffs_maxrs, S] = polyfit(log_x_D, log_max_rs, 1)

%% Slope of max of TKE
max_tke = max(abs(tke_1d_smooth));
log_max_tke = log(max_tke');
[coeffs_tke_max, S] = polyfit(log_x_D, log_max_tke, 1)
%% Reading the grid files

theta = linspace(0,2*pi,ntheta)';
numvar = 3;   
fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Plot the unscale reynolds stresses

figure;
hold on;

C = {'-','-','-','-','-','-','-','--','--','--','--','--','--','--'}; % Cell array of linestyle
color = {'y','c','k','b','r','g','m','m','g','r','b','k','c','y'}; % Cell array of linestyle

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

for i = 1:size(loc_planes)
    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(:,i)/scaling_factor_ud_dld_dx(i,1), ...
        'Linestyle', C{i}, 'color', color{i}, 'Linewidth',2);
    %plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i))), 'o', 'color', color{i}, 'Linewidth',2);
end

ylim([0, 0.7]);
xlim([0, 3]);

hXLabel = xlabel('$r/L_{d}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$-\langle u''_{x}u''_{r} \rangle/(U_{\infty}U_{d}dL_{d}/dx)$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$-\langle u''_{x}u''_{r} \rangle/K_{o}$','interpreter','latex','fontsize',20);

ax = gca;
ax.FontSize = 20; 

% xlim([0 100]);
% xticks([10 20 30 40 60 80 100]);

Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 25$';
Legend{3} = '$x/D = 30$';
Legend{4} = '$x/D = 35$';
Legend{5} = '$x/D = 40$';
Legend{6} = '$x/D = 45$';
Legend{7} = '$x/D = 50$';
Legend{8} = '$x/D = 55$';

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'rxr_normalized_uinfty_ud_dld.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'rxr_normalized_uinfty_ud_dld.eps'),'-depsc2','-r600');

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'rxr_normalized_ko.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'rxr_normalized_ko.eps'),'-depsc2','-r600');
%% Plotting the Reynolds stresses

figure;
hold on;

C = {'-','-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','g','r','b','k'}; % Cell array of linestyle


for i = 1:size(loc_planes,1)
    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(:,i)/max(tke_1d(:,i)), ...
        'Linestyle', C{i}, 'color', color{i}, 'Linewidth',2);
end

ylim([0, 0.2]);
xlim([0, 3]);
%% Plotting the Reynolds stresses 

figure;
hold on;

C = {'-','-','-','-','-','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','g','r','b','k'}; % Cell array of linestyle


reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');

for i = 1:size(loc_planes,1)
    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(:,i)/max(tke_1d(:,i)), ...
        'Linestyle', C{i}, 'color', color{i}, 'Linewidth',2);
end

ylim([0, 0.2]);
xlim([0, 3]);

%% Fit a function to single profile at x/D = 40

% yr = -reystress_uw_1d_smooth(:,8)/max(abs(reystress_uw_1d(:,8)));
% yr = -reystress_uw_1d_smooth(:,8)/TKE_centerline_loc_planes(8,2);
% 
% % eta = rc/LK_mean_loc_planes(8,2);
% eta = rc/LK_TKE_loc_planes(8,2);
% 
% figure;
% plot(eta, yr, 'k-', 'Linewidth',2);
% 
% yr_zero_padded = [0; yr];
% eta_zero_padded = [0; eta];
% 
% hold on;
% plot(eta_zero_padded, yr_zero_padded, 'k--', 'Linewidth',2);
% 
% alpha = 10;
% beta  = 0.3;
% gamma = 0.1;
% geta  = alpha*eta_zero_padded.*exp(-beta*eta_zero_padded.^2 - gamma*eta_zero_padded.^4);
% plot(eta_zero_padded, geta, 'b--', 'Linewidth',2);
% 
% % Fitting to g(eta) = alpha*eta*exp(-beta*(eta)^2))
% geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
% 
% x = eta_zero_padded;
% y = yr_zero_padded;
% [f1, gof] = fit(x,y,geta_fit, 'Start', [alpha beta gamma])
% 
% geta  = (f1.a)*eta_zero_padded.*exp(-(f1.b)*eta_zero_padded.^2 - (f1.c)*eta_zero_padded.^4);
% plot(eta_zero_padded, geta, 'r--', 'Linewidth',2);


%% Fit a function to all the profiles from x/D = 5 to x/D = 70 by normalizing by max(RS)

% close all;
% 
% alpha_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% beta_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% gamma_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% sse_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% rmse_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% rsquare_fit_maxrs_ld = zeros(size(loc_planes,1),1);
% 
% for i = 1:size(loc_planes,1)
%     yr = -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i)));
%     eta = rc/LK_mean_loc_planes(i,2);
%     disp(max(abs(reystress_uw_1d(:,i))));
%     disp(LK_mean_loc_planes(i,2));
%     yr_zero_padded = [0; yr];
%     y = yr_zero_padded;
%     eta_zero_padded = [0; eta];
%     x = eta_zero_padded;
%     [f1, gof] = fit(x,y,geta_fit, 'Start', [alpha beta gamma])
%     alpha_fit_maxrs_ld(i,1) = f1.a;
%     beta_fit_maxrs_ld(i,1) = f1.b;
%     gamma_fit_maxrs_ld(i,1) = f1.c;
%     rmse_fit_maxrs_ld(i,1) = gof.rmse;
%     sse_fit_maxrs_ld(i,1) = gof.sse;
%     rsquare_fit_maxrs_ld(i,1) = gof.rsquare;
% end

%% Plotting the Reynolds stresses with the fitted function

% close all;
% 
% figure;
% hold on;
% 
% C = {'-','-','-','-','-','-','-','--','--','--','--','--','--','--'}; % Cell array of linestyle
% color = {'y','c','k','b','r','g','m','m','g','r','b','k','c','y'}; % Cell array of linestyle
% 
% reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');
% 
% for i = 4:size(loc_planes,1)-2
%     plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i))), ...
%         'Linestyle', C{i}, 'color', color{i}, 'Linewidth',2);
% end
% geta  = (f1.a)*eta_zero_padded.*exp(-(f1.b)*eta_zero_padded.^2 - (f1.c)*eta_zero_padded.^4);
% plot(eta_zero_padded, geta, 'c-');
% ylim([0, 1.2]);
% xlim([0, 3]);

%% Fit a function to all the profiles from x/D = 5 to x/D = 70 by normalizing by K_centerline

% close all;
% 
% alpha_fit_k_lk = zeros(size(loc_planes,1),1);
% beta_fit_k_lk = zeros(size(loc_planes,1),1);
% gamma_fit_k_lk = zeros(size(loc_planes,1),1);
% sse_fit_k_lk = zeros(size(loc_planes,1),1);
% rmse_fit_k_lk = zeros(size(loc_planes,1),1);
% rsquare_fit_k_lk = zeros(size(loc_planes,1),1);
% 
% for i = 1:size(loc_planes,1)
%     yr = -reystress_uw_1d_smooth(:,i)/TKE_centerline_loc_planes(i,2);
%     eta = rc/LK_TKE_loc_planes(i,2);
%     disp(TKE_centerline_loc_planes(i,2));
%     disp(LK_TKE_loc_planes(i,2));
%     yr_zero_padded = [0; yr];
%     y = yr_zero_padded;
%     eta_zero_padded = [0; eta];
%     x = eta_zero_padded;
%     [f1, gof] = fit(x,y,geta_fit, 'Start', [alpha beta gamma])
%     alpha_fit_k_lk(i,1) = f1.a;
%     beta_fit_k_lk(i,1) = f1.b;
%     gamma_fit_k_lk(i,1) = f1.c;
%     rmse_fit_k_lk(i,1) = gof.rmse;
%     sse_fit_k_lk(i,1) = gof.sse;
%     rsquare_fit_k_lk(i,1) = gof.rsquare;
% end

%% Plotting the Reynolds stresses with the fitted function

% close all;
% 
% figure;
% hold on;
% 
% C = {'-','-','-','-','-','-','-','--','--','--','--','--','--','--'}; % Cell array of linestyle
% color = {'y','c','k','b','r','g','m','m','g','r','b','k','c','y'}; % Cell array of linestyle
% 
% reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess');
% 
% for i = 4:size(loc_planes,1)-2
%     plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i))), ...
%         'Linestyle', C{i}, 'color', color{i}, 'Linewidth',2);
% end
% geta  = (f1.a)*eta_zero_padded.*exp(-(f1.b)*eta_zero_padded.^2 - (f1.c)*eta_zero_padded.^4);
% plot(eta_zero_padded, geta, 'c-');
% ylim([0, 1.2]);
% xlim([0, 3]);

%% Plotting the RMS errors of two fits as a function of x/D

% figure;
% hold on;
% h1 = plot(loc_planes, rmse_fit_maxrs_ld(:,1), 'bo');
% h2 = plot(loc_planes, rmse_fit_k_lk(:,1), 'ro');
% 
% figure;
% hold on;
% h1 = plot(loc_planes, sse_fit_maxrs_ld(:,1), 'bo');
% h2 = plot(loc_planes, sse_fit_k_lk(:,1), 'ro');
% 
% figure;
% hold on;
% h1 = plot(loc_planes, rsquare_fit_maxrs_ld(:,1), 'b*');
% h2 = plot(loc_planes, rsquare_fit_k_lk(:,1), 'ro');


%% Fitting a single function to the whole set of data using inbuilt MATLAB tool by using RS

ytotal = [];
etotal = [];
for i = 1:size(loc_planes,1)
    yr = -reystress_uw_1d_smooth(:,i)/scaling_factor_ud_dld_dx(i);
    eta = rc/LK_mean_loc_planes(i,2);
    yr_zero_padded = [0; yr];
    eta_zero_padded = [0; eta];
    ytotal = [ytotal; yr_zero_padded];
    etotal = [etotal; eta_zero_padded];
end


alpha = 10;
beta  = 0.3;
gamma = 0.1;
geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
[f1_ud_dld_dx, gof_ud_dld_dx] = fit(etotal,ytotal,geta_fit, 'Start', [alpha beta gamma])
%% Fitting a single function to the whole set of data using inbuilt MATLAB tool by using TKE

ytotal = [];
etotal = [];
for i = 1:size(loc_planes,1)
    yr = -reystress_uw_1d_smooth(:,i)/max(tke_1d_smooth(:,i));
    eta = rc/LK_TKE_loc_planes(i,2);
    yr_zero_padded = [0; yr];
    eta_zero_padded = [0; eta];
    ytotal = [ytotal; yr_zero_padded];
    etotal = [etotal; eta_zero_padded];
end


alpha = 10;
beta  = 0.3;
gamma = 0.1;
geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
[f1_k_lk, gof_k_lk] = fit(etotal,ytotal,geta_fit, 'Start', [alpha beta gamma])

%% Fitting a single function to the whole set of data using inbuilt MATLAB tool by using TKE

ytotal = [];
etotal = [];
for i = 1:size(loc_planes,1)
    yr = -reystress_uw_1d_smooth(:,i)/max(tke_1d_smooth(:,i));
    eta = rc/LK_mean_loc_planes(i,2);
    yr_zero_padded = [0; yr];
    eta_zero_padded = [0; eta];
    ytotal = [ytotal; yr_zero_padded];
    etotal = [etotal; eta_zero_padded];
end


alpha = 10;
beta  = 0.3;
gamma = 0.1;
geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
[f1_k_ld, gof_k_ld] = fit(etotal,ytotal,geta_fit, 'Start', [alpha beta gamma])
%% Fitting a single function to the whole set of data using inbuilt MATLAB tool by using RS

ytotal = [];
etotal = [];
for i = 1:size(loc_planes,1)
    yr = -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i)));
    eta = rc/LK_mean_loc_planes(i,2);
    yr_zero_padded = [0; yr];
    eta_zero_padded = [0; eta];
    ytotal = [ytotal; yr_zero_padded];
    etotal = [etotal; eta_zero_padded];
end


alpha = 10;
beta  = 0.3;
gamma = 0.1;
geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
[f1_maxrs_ld, gof_maxrs_ld] = fit(etotal,ytotal,geta_fit, 'Start', [alpha beta gamma])
%% Fitting a single function to the whole set of data using inbuilt MATLAB tool by using RS

ytotal = [];
etotal = [];
for i = 1:size(loc_planes,1)
    yr = -reystress_uw_1d_smooth(:,i)/max(abs(reystress_uw_1d_smooth(:,i)));
    eta = rc/LK_TKE_loc_planes(i,2);
    yr_zero_padded = [0; yr];
    eta_zero_padded = [0; eta];
    ytotal = [ytotal; yr_zero_padded];
    etotal = [etotal; eta_zero_padded];
end


alpha = 10;
beta  = 0.3;
gamma = 0.1;
geta_fit = 'a*x*exp(-b*x^2-c*x^4)';
[f1_maxrs_lk, gof_maxrs_lk] = fit(etotal,ytotal,geta_fit, 'Start', [alpha beta gamma])

%% Plotting the following
close all;
h1 = loglog(LK_TKE(:,1), ud_Centerline(:,2), 'ko');
hold on;
% Slope of max of Ko/Ld
max_tke_ld = max(abs(tke_1d_smooth))'./LK_mean_loc_planes(:,2);
log_max_tke_ld = log(max_tke_ld);
[coeffs_tke_ld, S] = polyfit(log_x_D, log_max_tke_ld, 1)
coeffs_ud_calculated = coeffs_tke_ld+1;
y1 = 0.6*LK_mean_loc_planes(1:end,1).^(coeffs_ud_calculated(1));
h2 = loglog(LK_mean_loc_planes(1:end,1),y1, 'k--', 'Linewidth', 2);
