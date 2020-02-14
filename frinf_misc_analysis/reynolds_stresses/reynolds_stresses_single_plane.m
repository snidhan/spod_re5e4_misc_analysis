%% Written by Sheel Nidhan 
clear; clc;
%% Parameters

var1 = 'up';
var2 = 'wp';
var3 = 'vp';
loc = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
nstart = 1892600;
nend   = 2613200;
stride  = 100;
dir_in = '/mnt/RAID5/sheel/spod_re5e4/frinf/data_files_uniform/';
nr = 356;
ntheta = 258;
N = (nend-nstart)/stride + 1;

for x_loc = 1:size(loc,1)

% Reading the data files of radial velocity
disp(dir_in);
dir = strcat(dir_in, 'x_D_', int2str(loc(x_loc,1)),'/');

u = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var1, '/', var1, '_', num2str(num,'%08.f'), '_', int2str(loc(x_loc,1)), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    disp(fid);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            u(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

% Reading the data files of streamwise velocity

w = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var2, '/', var2, '_', num2str(num,'%08.f'), '_', int2str(loc(x_loc,1)), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            w(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

% Reading the files of azimuthal velocity

v = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(loc(x_loc,1)), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            v(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

%% Centering the velocity values of respective velocities
disp('centering the velocities');
u_centered = 0.5*(u(2:nr-1,2:ntheta-1,:) + u(1:nr-2,2:ntheta-1,:));
v_centered = 0.5*(v(2:nr-1,2:ntheta-1,:) + v(2:nr-1,1:ntheta-2,:));
w_centered = w(2:nr-1,2:ntheta-1,:);
disp('centered the velocities');
clear u; clear v; clear w;
%% Average of the velocity fields in time 

w_mean = squeeze(mean(w_centered,3));
v_mean = squeeze(mean(v_centered,3));
u_mean = squeeze(mean(u_centered,3));
disp('removed the mean velocities');
% w_mean_th_time = squeeze(mean(w_mean_theta,2));
% u_mean_th_time = squeeze(mean(u_mean_theta,2));
% v_mean_th_time = squeeze(mean(v_mean_theta,2));


%% Fluctuations of radial, streamwise and azimuthal velocity

u_fluc = u_centered - u_mean;
v_fluc = v_centered - v_mean;
w_fluc = w_centered - w_mean;

%% Calculating the averaged Reynolds stress over time

reystress_ww = w_fluc.*w_fluc;
reystress_ww_av = mean(reystress_ww,3);
reystress_ww_av_th = squeeze(mean(reystress_ww_av,2));

reystress_uw = u_fluc.*w_fluc;
reystress_uw_av = mean(reystress_uw,3);
reystress_uw_av_th = squeeze(mean(reystress_uw_av,2));


reystress_uu = u_fluc.*u_fluc;
reystress_uu_av = mean(reystress_uu,3);
reystress_uu_av_th = squeeze(mean(reystress_uu_av,2));


reystress_vv = v_fluc.*v_fluc;
reystress_vv_av = mean(reystress_vv,3);
reystress_vv_av_th = squeeze(mean(reystress_vv_av,2));

reystress_wv = w_fluc.*v_fluc;
reystress_wv_av = mean(reystress_wv,3);
reystress_wv_av_th = squeeze(mean(reystress_wv_av,2));

reystress_uv = u_fluc.*v_fluc;
reystress_uv_av = mean(reystress_uv,3);
reystress_uv_av_th = squeeze(mean(reystress_uv_av,2));

disp('Calculated the reynolds stresses');

clear u_fluc; clear v_fluc; clear w_fluc;
clear u_centered; clear v_centered; clear w_centered;
clear reystress_uu; clear reystress_vv; clear reystress_ww;
clear reystress_uv; clear reystress_wv; clear reystress_uw;
%% Reading the grid files

% theta = linspace(0,2*pi,ntheta)';
% numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
%               % numvar = 4 (velocity and density is used for kernel)
% 
% fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
% %fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid
% 
% D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
% 
% r = D(1:end-9,2);
% 
% for i = 1:size(r,1)-2
%     rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
% end


%% Saving the mat files

% dir_out = '/mnt/RAID5/sheel/spod_re5e4/frinf/reystresses/';
dir_out = './';
file_out = strcat(dir_out, 'reystress_x_D_', int2str(loc(x_loc,1)), '.mat');
save(file_out, '-v7.3');


% %% Creating polar field
% viridis=viridis();
% 
% figure;
% [C1,h1] = polarcont(rc, theta, reystress_uw_av(1:nr-2,:), 10);
% axis equal
% xlim([-15 15]);
% ylim([-15 15]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{x}u_{r}>$','interpreter','latex','fontsize',15);
%   
% 
% figure;
% [C2,h2] = polarcont(rc, theta, reystress_uu_av(1:nr-2,:), 10);
% axis equal
% xlim([-15 15]);
% ylim([-15 15]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{r}u_{r}>$','interpreter','latex','fontsize',15);
% 
% figure;
% [C3,h3] = polarcont(rc, theta, reystress_ww_av(1:nr-2,:), 10);
% axis equal
% xlim([-8 8]);
% ylim([-8 8]);    
% colormap(viridis);
% colorbar;
% hXLabel = xlabel('$y$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$z$','interpreter','latex','fontsize',15);
% hTitle = title('$<u_{x}u_{x}>$','interpreter','latex','fontsize',15);
% 
% 
% %% Creating the line plot of the reynolds stress
% 
% figure;
% hold on;
% h4 = plot(rc, reystress_uw_av_th(1:nr-2,1), 'k','Linewidth',2);
% h5 = plot(rc, reystress_uu_av_th(1:nr-2,1), 'r','Linewidth',2);
% h6 = plot(rc, reystress_ww_av_th(1:nr-2,1), 'b','Linewidth',2);


end
