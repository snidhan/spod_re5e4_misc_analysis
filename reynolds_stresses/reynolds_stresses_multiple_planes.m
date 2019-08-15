clear; clc; close all;
%% Parameters

var1 = 'up';
var2 = 'wp';
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
%loc = [50];
nstart = 1892600;
nend   = 2613200;
stride  = 100;
dir_in_planes = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/reystresses/';
nr = 356;
ntheta = 258;
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

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1)));
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/TKE_centerline.dat';
TKE_Centerline   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-TKE_Centerline(:,1)));
    TKE_centerline_loc_planes(i,2) = TKE_Centerline(idx,2);
    TKE_centerline_loc_planes(i,1) = TKE_Centerline(idx,1);
end


filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);

for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
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

%% Plotting the Reynolds stresses

figure; 
hold on

for i = 8:size(loc_planes,1)
    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr-2,i)/max(abs(reystress_uw_1d(1:nr-2,i))) ...
           ,'k--','Linewidth',2);
end

%% Creating polar field
viridis=viridis();
x_D = 11;
disp(loc_planes(x_D));

figure;
[C1,h1] = polarcont(rc, theta, -reystress_uw_2d(1:nr-2,:,x_D), 10);
axis equal
xlim([-8 8]);
ylim([-8 8]);    
colormap(viridis);
colorbar;
hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
hTitle = title('$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
  
figure;
[C2,h2] = polarcont(rc, theta, reystress_uu_2d(1:nr-2,:,x_D), 10);
axis equal
xlim([-8 8]);
ylim([-8 8]);    
colormap(viridis);
colorbar;
hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
hTitle = title('$<u_{r}u_{r}>$','interpreter','latex','fontsize',15);

figure;
[C3,h3] = polarcont(rc, theta, reystress_ww_2d(1:nr-2,:,x_D), 10);
axis equal
xlim([-8 8]);
ylim([-8 8]);    
colormap(viridis);
colorbar;
hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
hTitle = title('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);
