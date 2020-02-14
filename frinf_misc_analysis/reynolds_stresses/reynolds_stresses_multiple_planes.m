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
reystress_vv_2d = zeros(nr,ntheta,size(loc_planes,1));


reystress_uw_1d = zeros(nr,size(loc_planes,1));
reystress_uu_1d = zeros(nr,size(loc_planes,1));
reystress_ww_1d = zeros(nr,size(loc_planes,1));
reystress_vv_1d = zeros(nr,size(loc_planes,1));
tke_1d          = zeros(nr,size(loc_planes,1));
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
    reystress_vv_2d(:,:,x_loc_planes) =  reystress_vv_av;

    
    reystress_uw_1d(:,x_loc_planes) =  reystress_uw_av_th;
    reystress_ww_1d(:,x_loc_planes) =  reystress_ww_av_th;
    reystress_uu_1d(:,x_loc_planes) =  reystress_uu_av_th;
    reystress_vv_1d(:,x_loc_planes) =  reystress_vv_av_th;
    tke_1d(:,x_loc_planes)          =  reystress_vv_1d(:,x_loc_planes) + reystress_uu_1d(:,x_loc_planes) + reystress_ww_1d(:,x_loc_planes); 
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



%% Evaluating the value of RHS of RANS equation using the data

% count = 1;
% for x_loc_planes = 1:size(loc_planes,1)
%     reystress_uw_x_D = reystress_uw_1d(:,x_loc_planes);
%     reystress_uw_x_D_smooth = smoothdata(reystress_uw_x_D,'loess');
%     max_uw(count) = max(abs(reystress_uw_x_D));
%     term_1(count) = 0.5*(reystress_uw_x_D_smooth(3)/rc(3) + reystress_uw_x_D_smooth(2)/rc(2));
%     term_2(count) = (reystress_uw_x_D_smooth(3)-reystress_uw_x_D(2))/(rc(3)-rc(2));
%     count = count + 1;
% end
% 
% figure;
% h1 = loglog(loc_planes, -term_1, 'ko');
% hold on;
% h2 = loglog(loc_planes, max_uw, 'ko');

                       
% loc = LK_TKE_loc_planes(16:22,1)';                         
% log_max_uw = log(max_uw(16:22));
% log_loc =  log(loc);
% [coeffs, S] = polyfit(log_loc, log_max_uw, 1);
% y_fitted = polyval(coeffs, log_loc);
% 
% 
% count = 1;
% figure;
% hold on;
% for x_loc_planes = 1:12
%     reystress_uw_x_D = reystress_uw_1d(:,x_loc_planes);
%     f = fit(rc/LK_mean_loc_planes(x_loc_planes,2),-reystress_uw_x_D,'smoothingspline');
%     val(count) = (f(rc(74))-f(rc(64)))/(rc(74)-rc(64));
%     plot(rc/LK_mean_loc_planes(x_loc_planes,2), -reystress_uw_x_D)
%     plot(f,rc/LK_mean_loc_planes(x_loc_planes,2),-reystress_uw_x_D);
%   [d1, d2] = differentiate(f,rc);
%   rhs = d1./rc;
%   frhs  = fit(rc(3:end), rhs(3:end), 'smoothingspline');
%   plot(frhs,rc,rhs);
%     count = count + 1;
% end
% 
% 
% figure;
% loglog(loc_planes, val, 'ko')
% 
% loc = LK_TKE_loc_planes(4:12,1)';                         
% log_max_uw = log(val(4:12));
% log_loc =  log(loc);
% [coeffs, S] = polyfit(log_loc, log_max_uw, 1);
% y_fitted = polyval(coeffs, log_loc);


%% Fitting curve to Reynolds stress data 

close all;

% eta = rc/LK_mean_loc_planes(6,2);
% plot(eta,-reystress_uw_1d(:,6),'k-','Linewidth',2);
% hold on;
% 
F = @(x,xdata)x(5)*xdata.*(1+x(1)*xdata.^2+x(2)*xdata.^4).*exp(-x(3)*xdata.^2-x(4)*xdata.^4);

x0 = [0.1; 0.1; 0.1; 0.1; 0.4];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,eta,-reystress_uw_1d(:,6));
% 
% plot(eta,F(x,eta),'k--', 'Linewidth',2);
% 
% F = @(x,xdata)x(2)*xdata.*exp(-x(1)*xdata.^2);
% 
% x0 = [0.1,1];
% [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,eta,-reystress_uw_1d(:,4)/max(-reystress_uw_1d(:,4)));
% 
% plot(eta,F(x,eta),'r--', 'Linewidth',2);
% 
% hold off



count = 1;
for i = 4:12
    reystress_20_60(354*(count-1)+1:354*count,1) = -reystress_uw_1d(:,i)./max(-reystress_uw_1d(:,i));
    count  = count + 1;
end

reystress_combined = reystress_20_60(:);

count = 1;
for i = 4:12
    eta_combined(354*(count-1)+1:354*count,1) = rc./LK_mean_loc_planes(i,2);
    count  = count + 1;
end

F =@(x,xdata)xdata.*(1+x(1)*xdata.^2+x(2)*xdata.^4).*exp(-x(3)*xdata.^2-x(4)*xdata.^4); 
x0 = [-0.5; 0.08; -0.25; 0.33];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,eta_combined,reystress_combined);
fitdata = F(x,eta_combined);

% plot(eta_combined(1:354),reystress_combined(1:354),'k-','Linewidth',2);
hold on;
plot(eta_combined(355:708),fitdata(355:708),'r-','Linewidth',2);



%% Similarity analysis the Reynolds stresses

dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
nr = 354; ntheta = 256;

save(strcat(dirout, 'similarity_uxur.mat'), 'LK_TKE_loc_planes', 'reystress_uw_1d', 'LK_mean_loc_planes','rc', 'nr');
save(strcat(dirout, 'similarity_tke.mat'), 'LK_TKE_loc_planes',  'tke_1d', 'LK_mean_loc_planes','rc', 'nr');

% 
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
% reystress_uw_1d_Ud = zeros(size(rc,1),size(reystress_uw_1d,2));
% 
% for i = 1:size(reystress_uw_1d,2)
%     reystress_uw_1d_Ud(:,i) = reystress_uw_1d(:,i)/(ud_centerline_loc_planes(i,2)^2);
% end

%% Scaling Reynolds stress with K
reystress_uw_1d_tke = zeros(size(rc,1),size(reystress_uw_1d,2));
for i = 1:size(reystress_uw_1d,2)
    reystress_uw_1d_tke(:,i) = reystress_uw_1d(:,i)/(TKE_centerline_loc_planes(i,2));
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with TKE

% wake_width_tke = LK_TKE_loc_planes(5:20,2);                       
% loc = LK_TKE_loc_planes(5:20,1);                         
% log_wake_width = log(wake_width_tke);
% log_loc =  log(loc);
% [coeffs, S] = polyfit(log_loc, log_wake_width, 1);
% y_fitted = polyval(coeffs, log_wake_width);
% 
% % Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
% ddelta_dx = coeffs(1)*wake_width_tke./loc;
% scaling_factor = (TKE_centerline_loc_planes(5:20,2).^0.5).*ddelta_dx; 
% count = 1;
% 
% reystress_uw_1d_tkeddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
% for i = 5:20
%     reystress_uw_1d_tkeddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
%     count = count + 1;
% end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with mean
wake_width_mean_20_60 = LK_mean_loc_planes(4:12,2);                       
loc_20_60 = LK_mean_loc_planes(4:12,1);                         
log_wake_width = log(wake_width_mean_20_60);
log_loc =  log(loc_20_60);
[coeffs_20_60, S] = polyfit(log_loc, log_wake_width, 1);
% y_fitted = polyval(coeffs, log_wake_width); %#ok<*NASGU>

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs_20_60(1)*wake_width_mean_20_60./loc_20_60;
scaling_factor_20_60 = (ud_centerline_loc_planes(4:12,2)).*ddelta_dx; 
count = 1;

wake_width_mean_65_100 = LK_mean_loc_planes(13:20,2);                       
loc_65_100 = LK_mean_loc_planes(13:20,1);                         
log_wake_width = log(wake_width_mean_65_100);
log_loc =  log(loc_65_100);
[coeffs_65_100, S] = polyfit(log_loc, log_wake_width, 1); %#ok<*ASGLU>
% y_fitted = polyval(coeffs, log_wake_width);

% Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
ddelta_dx = coeffs_65_100(1)*wake_width_mean_65_100./loc_65_100;
scaling_factor_65_100 = (ud_centerline_loc_planes(13:20,2)).*ddelta_dx; 
count = 1;

scaling_factor = [scaling_factor_20_60; scaling_factor_65_100];
reystress_uw_1d_udddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));

for i = 4:20
    reystress_uw_1d_udddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
    count = count + 1;
end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 but with ud(Lk/Ld)dLd/dx
% wake_width_mean_20_60 = LK_mean_loc_planes(5:12,2);                       
% loc_20_60 = LK_mean_loc_planes(5:12,1);                         
% log_wake_width = log(wake_width_mean_20_60);
% log_loc =  log(loc_20_60);
% [coeffs_20_60, S] = polyfit(log_loc, log_wake_width, 1);
% y_fitted = polyval(coeffs, log_wake_width);
% 
% % Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
% ddelta_dx = coeffs(1)*wake_width_mean_20_60./loc_20_60;
% scaling_factor_20_60 = (ud_centerline_loc_planes(5:12,2).*(LK_TKE_loc_planes(5:12,2)./LK_mean_loc_planes(5:12,2))).*ddelta_dx; 
% count = 1;
% 
% wake_width_mean_65_100 = LK_mean_loc_planes(13:20,2);                       
% loc_65_100 = LK_mean_loc_planes(13:20,1);                         
% log_wake_width = log(wake_width_mean_65_100);
% log_loc =  log(loc_65_100);
% [coeffs_65_100, S] = polyfit(log_loc, log_wake_width, 1);
% y_fitted = polyval(coeffs, log_wake_width);
% 
% % Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
% ddelta_dx = coeffs(1)*wake_width_mean_65_100./loc_65_100;
% scaling_factor_65_100 = (ud_centerline_loc_planes(13:20,2).*(LK_TKE_loc_planes(13:20,2)./LK_mean_loc_planes(13:20,2))).*ddelta_dx; 
% count = 1;
% 
% scaling_factor = [scaling_factor_20_60; scaling_factor_65_100];
% 
% reystress_uw_1d_udlk_ld_ddelta_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
% for i = 5:20
%     reystress_uw_1d_udlk_ld_ddelta_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
%     count = count + 1;
% end

%% Scaling of Reynolds stress with Eq. 7.10 of Dairay et al. 2015 mixed
% wake_width_tke = LK_TKE_loc_planes(5:20,2);                       
% loc = LK_TKE_loc_planes(5:20,1);                         
% log_wake_width = log(wake_width_tke);
% log_loc =  log(loc);
% [coeffs, S] = polyfit(log_loc, log_wake_width, 1);
% y_fitted = polyval(coeffs, log_wake_width);
% 
% % Calculating d\delta/dx = \delta/x * d(log \delta)/d (log x)
% ddelta_dx = coeffs(1)*wake_width_tke./loc;
% scaling_factor = (ud_centerline_loc_planes(5:20,2)).*ddelta_dx; 
% count = 1;
% 
% reystress_uw_1d_udddeltak_dx = zeros(size(rc,1),size(reystress_uw_1d,2));
% for i = 5:20
%     reystress_uw_1d_udddeltak_dx(:,i) =  reystress_uw_1d(:,i)/(scaling_factor(count));  
%     count = count + 1;
% end
% 
% 
% dirout = '/home/sheel/Work/codes/spod_re5e4_misc_analysis/spod_plots/files/';
% nr = 354; ntheta = 256;

% save(strcat(dirout, 'scaling_uxur.mat'), 'reystress_uw_1d_Ud', 'reystress_uw_1d_tke', 'reystress_uw_1d_tkeddelta_dx', ...
%                                          'reystress_uw_1d_udddeltak_dx', 'reystress_uw_1d_udddelta_dx', 'reystress_uw_1d_udlk_ld_ddelta_dx', ...
%                                           'LK_TKE_loc_planes', 'LK_mean_loc_planes', 'rc', 'nr');

save(strcat(dirout, 'scaling_uxur.mat'), 'reystress_uw_1d_tke', ...
                                       'reystress_uw_1d_udddelta_dx', ...
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