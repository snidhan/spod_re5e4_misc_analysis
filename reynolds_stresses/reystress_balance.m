%% Name: Sheel Nidhan
%  To calculate the Reynolds stresses and its similarity in the 'non-equilibrium' region
clear; clc; close all;

%% Parameters

loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
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
%% Loading the mean velocity profile
 
close all;

load('./similarity_w.mat');
mean_w_1d_smooth       = smoothdata(mean_w_1d, 'loess');
ud_1d_smooth = 1 - mean_w_1d_smooth;

max_ud_1d_smooth = max(ud_1d_smooth);
x_D_ud  = linspace(5,100,20)';
x_D_ud(21,1) = 110;
x_D_ud(22,1) = 120;


log_x_loc = log(x_D_ud(4:12));
log_max_ud = log(max_ud_1d_smooth(4:12))';
[coeffs_max_ud_60, S] = polyfit(log_x_loc, log_max_ud, 1)

log_x_loc = log(x_D_ud(16:end));
log_max_ud = log(max_ud_1d_smooth(16:end))';
[coeffs_max_ud_120, S] = polyfit(log_x_loc, log_max_ud, 1)

log_x_loc = log(x_D_ud(4:22));
log_max_ud = log(max_ud_1d_smooth(4:22))';
[coeffs_max_ud_total, S] = polyfit(log_x_loc, log_max_ud, 1)


h1 = loglog(x_D_ud, max_ud_1d_smooth(1:22)', 'ko');
hold on;
y1 = 0.8*x_D_ud(4:12,1).^(coeffs_max_ud_60(1));
h2 = loglog(x_D_ud(4:12,1),y1, 'k-', 'Linewidth', 2);
y2 = 0.6*x_D_ud(16:22,1).^(coeffs_max_ud_120(1));
h3 = loglog(x_D_ud(16:22,1),y2, 'r-', 'Linewidth', 2);
y3 = 0.5*x_D_ud(4:22,1).^(coeffs_max_ud_total(1));
h3 = loglog(x_D_ud(4:22,1),y3, 'k--', 'Linewidth', 2);
h4=text(30, 0.05,'$x^{-0.85}$','interpreter','latex','FontSize', 15);
h5=text(100, 0.02,'$x^{-0.80}$','interpreter','latex','FontSize', 15);


hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$U_{d}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$-\langle u''_{x}u''_{r} \rangle/K_{o}$','interpreter','latex','fontsize',20);

ax = gca;
ax.FontSize = 20; 

xlim([0 120]);
% xticks([10 20 30 40 60 80 100]);

% Legend{1} = '$max(U_{defect})$';
% Legend{2} = '$U_{defect}(r=0)$';
% 
% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 20;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'ud_max_variation.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'ud_max_variation.eps'),'-depsc2','-r600');




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

%% Plotting the balance at a given location

w_x_D_50 = mean_w_1d_smooth(:,10);
uw_x_D_50 = reystress_uw_1d_smooth(:,7);

figure;
plot(rc, uw_x_D_50, 'ko')
figure;
plot(rc, w_x_D_50, 'ro')

%% LHS U_{\infty}*dUd/dx

U_inf = 1;

ld = LK_mean_loc_planes(7,2);
ud = ud_centerline_loc_planes(7,2);
x_loc = 50;
alpha = -0.91;
beta  = 0.39;

Ud = U_inf - w_x_D_50;
Uo = U_inf - min(w_x_D_50(:,1));

f_eta = Ud/Uo;
eta  = rc/LK_mean_loc_planes(7,2);

count = 1;
for i = 2:size(f_eta,1)-1
    df = f_eta(i+1,1) - f_eta(i-1,1);
    d_eta = eta(i+1,1)   - eta(i-1,1);
    df_deta(count,1) = df/d_eta;
    count = count + 1;
end

% LHS can be simplified into: \alpha*Uo/x_loc*f_eta - \beta*Uo/x*\eta*df/d\eta

lhs_term1 = alpha*ud/x_loc.*f_eta(2:end-1);
lhs_term2 = beta*ud/x_loc.*eta(2:end-1).*df_deta;

lhs_total = lhs_term1 - lhs_term2;

%% RHS : 1/r*d(r*uxur)/dr

r_times_uxur = rc.*uw_x_D_50;

count = 1;

for i = 2:size(uw_x_D_50,1)-1
    dnum = r_times_uxur(i+1,1) - r_times_uxur(i-1,1);
    dden = rc(i+1,1) - rc(i-1,1);
    rhs_total(count,1) = (1/rc(i,1))*dnum/dden;
    count = count + 1;
end

%% Plotting the balance of LHS and RHS

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure;
h1 = plot(rc(2:end-1,1), abs(rhs_total), 'ko');
hold on;
h2 = plot(rc(2:end-1,1), abs(lhs_total), 'ro');


hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('Equation Term','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$-\langle u''_{x}u''_{r} \rangle/K_{o}$','interpreter','latex','fontsize',20);

ax = gca;
ax.FontSize = 20; 

% xlim([0 100]);
% xticks([10 20 30 40 60 80 100]);

Legend{1} = '$\frac{1}{r}\frac{\partial}{\partial r}(r\langle u''_{x}u''_{r} \rangle)$';
Legend{2} = '$U_{\infty}\frac{\partial}{\partial x}(U_{\infty}-U)$';

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'balance_rans_v1.0.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'balance_rans_v1.0.eps'),'-depsc2','-r600');


%% Decay of MKE as a function of x/D

log_x_loc = log(mke_area(17:end,1));
log_mke = log(mke_area(17:end,2));
[coeffs_mke_total, S] = polyfit(log_x_loc, log_mke, 1)

product = (ud_Centerline(:,2).^2).*(LK_mean(:,2).^2);
log_product = log(product(17:end,1));
[coeffs_product_total, S] = polyfit(log_x_loc, log_product, 1)

close all;
figure;
h1 = loglog(mke_area(:,1), mke_area(:,2), 'ko');
hold on;
h2 = loglog(mke_area(:,1), product(:,1), 'ro');
y1 = 0.20*mke_area(17:end,1).^(coeffs_mke_total(1));
h3 = loglog(mke_area(17:end,1), y1, 'b--', 'Linewidth', 2);


Legend{1} = 'Area-integrated MKE';
Legend{2} = '$U_{d}^{2}L_{d}^{2}$';
Legend{3} = '$x^{-0.89}$';


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

ax = gca;
ax.FontSize = 20;

xlim([0 120]);
xticks([10 60 120]);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'mke_decay_rate.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'mke_decay_rate.eps'),'-depsc2','-r600');
