%% Name: Sheel Nidhan
%  To calculate the Reynolds stresses and its similarity
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

reystress_uw_1d = zeros(nr,size(loc_planes,1));
reystress_uu_1d = zeros(nr,size(loc_planes,1));
reystress_ww_1d = zeros(nr,size(loc_planes,1));

%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);
LK_TKE_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1))); %#ok<*NCOMMA>
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);
LK_mean_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_centerline.dat';
TKE_Centerline   = importdata(filename);
TKE_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-TKE_Centerline(:,1)));
    TKE_centerline_loc_planes(i,2) = TKE_Centerline(idx,2);
    TKE_centerline_loc_planes(i,1) = TKE_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/UX_rms_centerline.dat';
ux_Centerline   = importdata(filename);
ux_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ux_Centerline(:,1)));
    ux_centerline_loc_planes(i,2) = ux_Centerline(idx,2);
    ux_centerline_loc_planes(i,1) = ux_Centerline(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Defect_centerline.dat';
ud_Centerline   = importdata(filename);
ud_centerline_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-ud_Centerline(:,1)));
    ud_centerline_loc_planes(i,2) = ud_Centerline(idx,2);
    ud_centerline_loc_planes(i,1) = ud_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/MKE_areaI_FINF.dat';
mke_area   = importdata(filename);
mke_area_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-mke_area(:,1)));
    mke_area_loc_planes(i,2) = mke_area(idx,2);
    mke_area_loc_planes(i,1) = mke_area(idx,1);
end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_areaI_FINF.dat';
tke_area   = importdata(filename);
tke_area_loc_planes = zeros(22,2);
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
    
    reystress_uw_1d(:,x_loc_planes) =  reystress_uw_av_th;
    reystress_ww_1d(:,x_loc_planes) =  reystress_ww_av_th;
    reystress_uu_1d(:,x_loc_planes) =  reystress_uu_av_th;
end

%% Reading the grid files

theta = linspace(0,2*pi,ntheta)';
numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Similarity analysis the Reynolds stresses

dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
nr = 354; ntheta = 256;

save(strcat(dirout, 'similarity_uxur.mat'), 'LK_TKE_loc_planes', 'reystress_uw_1d', 'LK_mean_loc_planes','rc', 'nr');

%figure; 
% hold on;
% C = {'k','b','r','g','m','c','k--','b--','r--','g--','m--','c--','k-.','b-.','r-.','g-.','m-.','c-.'}; % Cell array of colros.
% count = 1;

% Legend = cell(7,1);
% Legend{1} = 'x/D = 70';
% Legend{2} = 'x/D = 75';
% Legend{3} = 'x/D = 80';
% Legend{4} = 'x/D = 85';
% Legend{5} = 'x/D = 90';
% Legend{6} = 'x/D = 95';
% Legend{7} = 'x/D = 100';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 95';
% Legend{10} = 'x/D = 100';

%for i = 14:20
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), C{count}, 'Linewidth',2);
%    count = count + 1;
%end



% ylim([0 1.2])
% xlim([0 3]);

% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/max($<u_{x}u_{r}>_{r}$)','interpreter','latex','fontsize',15);
% %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('unnormalized_uw_x_D_70_100_lk', '.png'),'-dpng','-r600');  
% print(gcf,strcat('unnormalized_uw_x_D_70_100_lk', '.eps'),'-depsc','-r600');
% close;

%% Scaling for Reynolds stresses in streamwise direction

% figure; 
% hold on;
% C = {'k','b','r','g','m','c','k--','b--','r--','g--','m--','c--','k-.','b-.','r-.','g-.','m-.','c-.'}; % Cell array of colros.
% Legend = cell(7,1);
% Legend{1} = 'x/D = 20';
% Legend{2} = 'x/D = 25';
% Legend{3} = 'x/D = 30';
% Legend{4} = 'x/D = 35';
% Legend{5} = 'x/D = 40';
% Legend{6} = 'x/D = 45';
% Legend{7} = 'x/D = 50';
% Legend{8} = 'x/D = 90';
% Legend{9} = 'x/D = 95';
% Legend{10} = 'x/D = 100';


%% Scaling Reynolds stresses with Ud^2
reystress_uw_1d_Ud = zeros(size(rc,1),size(reystress_uw_1d,2));

for i = 1:size(reystress_uw_1d,2)
    reystress_uw_1d_Ud(:,i) = reystress_uw_1d(:,i)/(ud_centerline_loc_planes(i,2)^2);
end

%% Scaling Reynolds stress with K
reystress_uw_1d_tke = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 1:size(reystress_uw_1d,2)
    reystress_uw_1d_tke(:,i) = reystress_uw_1d(:,i)/(TKE_centerline_loc_planes(i,2));
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with TKE

wake_width_tke = LK_TKE_loc_planes(5:20,2);                       
loc = LK_TKE_loc_planes(5:20,1);                         
log_wake_width = log(wake_width_tke);
log_loc =  log(loc);
[coeffs, S] = polyfit(log_loc, log_wake_width, 1);
y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_tke./loc;
scaling_factor = (TKE_centerline_loc_planes(5:20,2).^0.5).*ddelta_dx; 
count = 1;

reystress_uw_1d_tkeddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 5:20
    reystress_uw_1d_tkeddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
    count = count + 1;
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with mean
wake_width_mean_20_60 = LK_mean_loc_planes(5:12,2);                       
loc_20_60 = LK_mean_loc_planes(5:12,1);                         
log_wake_width = log(wake_width_mean_20_60);
log_loc =  log(loc_20_60);
[coeffs_20_60, S] = polyfit(log_loc, log_wake_width, 1);
y_fitted = polyval(coeffs, log_wake_width); %#ok<*NASGU>

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_mean_20_60./loc_20_60;
scaling_factor_20_60 = (ud_centerline_loc_planes(5:12,2)).*ddelta_dx; 
count = 1;

wake_width_mean_65_100 = LK_mean_loc_planes(13:20,2);                       
loc_65_100 = LK_mean_loc_planes(13:20,1);                         
log_wake_width = log(wake_width_mean_65_100);
log_loc =  log(loc_65_100);
[coeffs_65_100, S] = polyfit(log_loc, log_wake_width, 1); %#ok<*ASGLU>
y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_mean_65_100./loc_65_100;
scaling_factor_65_100 = (ud_centerline_loc_planes(13:20,2)).*ddelta_dx; 
count = 1;

scaling_factor = [scaling_factor_20_60; scaling_factor_65_100];

reystress_uw_1d_udddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 5:20
    reystress_uw_1d_udddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
    count = count + 1;
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with ud(Lk/Ld)dLd/dx
wake_width_mean_20_60 = LK_mean_loc_planes(5:12,2);                       
loc_20_60 = LK_mean_loc_planes(5:12,1);                         
log_wake_width = log(wake_width_mean_20_60);
log_loc =  log(loc_20_60);
[coeffs_20_60, S] = polyfit(log_loc, log_wake_width, 1);
y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_mean_20_60./loc_20_60;
scaling_factor_20_60 = (ud_centerline_loc_planes(5:12,2).*(LK_TKE_loc_planes(5:12,2)./LK_mean_loc_planes(5:12,2))).*ddelta_dx; 
count = 1;

wake_width_mean_65_100 = LK_mean_loc_planes(13:20,2);                       
loc_65_100 = LK_mean_loc_planes(13:20,1);                         
log_wake_width = log(wake_width_mean_65_100);
log_loc =  log(loc_65_100);
[coeffs_65_100, S] = polyfit(log_loc, log_wake_width, 1);
y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_mean_65_100./loc_65_100;
scaling_factor_65_100 = (ud_centerline_loc_planes(13:20,2).*(LK_TKE_loc_planes(13:20,2)./LK_mean_loc_planes(13:20,2))).*ddelta_dx; 
count = 1;

scaling_factor = [scaling_factor_20_60; scaling_factor_65_100];

reystress_uw_1d_udlk_ld_ddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 5:20
    reystress_uw_1d_udlk_ld_ddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
    count = count + 1;
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 mixed
wake_width_tke = LK_TKE_loc_planes(5:20,2);                       
loc = LK_TKE_loc_planes(5:20,1);                         
log_wake_width = log(wake_width_tke);
log_loc =  log(loc);
[coeffs, S] = polyfit(log_loc, log_wake_width, 1);
y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs(1)*wake_width_tke./loc;
scaling_factor = (ud_centerline_loc_planes(5:20,2)).*ddelta_dx; 
count = 1;

reystress_uw_1d_udddeltak_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 5:20
    reystress_uw_1d_udddeltak_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
    count = count + 1;
end


dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
nr = 354; ntheta = 256;

save(strcat(dirout, 'scaling_uxur.mat'), 'reystress_uw_1d_Ud', 'reystress_uw_1d_tke', 'reystress_uw_1d_tkeddelta_dx', ...
                                         'reystress_uw_1d_udddeltak_dx', 'reystress_uw_1d_udddelta_dx', 'reystress_uw_1d_udlk_ld_ddelta_dx', ...
                                          'LK_TKE_loc_planes', 'LK_mean_loc_planes', 'rc', 'nr');


%% Plot parameters
% hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('-$<u_{x}u_{r}>$/($u_{x}^{''})^{2}$','interpreter','latex','fontsize',15); %%check
% %hTitle = title('Variation of $C_{p}$ vs $\theta$','interpreter','latex','fontsize',15);

% hLegend = legend(Legend);
% hLegend.Interpreter = 'Latex';
% hLegend.FontSize = 10;
% hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('ux_uw_x_D_20_50_lk', '.png'),'-dpng','-r600');  
% print(gcf,strcat('ux_uw_x_D_20_50_lk', '.eps'),'-depsc','-r600');

%% Creating polar field
% viridis=viridis();
% x_D = 20;
% disp(loc_planes(x_D));
% 
% figure;
% [C1,h1] = polarcont(rc, theta, -reystress_uw_2d(1:nr-2,:,x_D), 10);
% axis equal
% xlim([-8 8]);
% ylim([-8 8]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
%   
% figure;
% [C2,h2] = polarcont(rc, theta, reystress_uu_2d(1:nr-2,:,x_D), 10);
% axis equal
% xlim([-8 8]);
% ylim([-8 8]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{r}u_{r}>$','interpreter','latex','fontsize',15);
% 
% figure;
% [C3,h3] = polarcont(rc, theta, reystress_ww_2d(1:nr-2,:,x_D), 10);
% axis equal
% xlim([-8 8]);
% ylim([-8 8]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);