%% %% Written by Sheel Nidhan
%  This code is used to calculate the TKE of the field that is obtained from actual data

clear; clc;
%% Reading the datafiles of u_x velocity field

var = 'wp';
loc = '50';
nstart = 1892600;
nend   = 1892600;
stride  = 100;
dir = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_50/';
nr = 356;
ntheta = 258;
N = (nend-nstart)/stride + 1;
vel_ux = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var, '/', var, '_', num2str(num,'%08.f'), '_', loc, '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            vel_ux(i,j,n) = a((j-1)*nr + i, 1);
        end
    end
    
end

% var = 'up';
% loc = '50';
% nstart = 1892600;
% nend   = 2613200;
% stride  = 100;
% dir = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_50/';
% nr = 356;
% ntheta = 258;
% N = (nend-nstart)/stride + 1;
% vel_ur = zeros(nr,ntheta,N);
% 
% for n = 1:N
%     num = (n-1)*stride + nstart;
%     filename = strcat(dir, var, '/', var, '_', num2str(num,'%08.f'), '_', loc, '_', 'uniform_pchip.res');
%     disp(filename);
%     fid = fopen(filename);
%     h = fread(fid,0,'*uint64'); % May need adjusting
%     a = fread(fid, nr*ntheta, '*double');
%     fclose(fid);
%    for j = 1:ntheta
%         for i = 1:nr
%             vel_ur(i,j,n) = a((j-1)*nr + i, 1);
%         end
%     end
%     
% end
% 
% var = 'vp';
% loc = '50';
% nstart = 1892600;
% nend   = 2613200;
% stride  = 100;
% dir = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_50/';
% nr = 356;
% ntheta = 258;
% N = (nend-nstart)/stride + 1;
% vel_utheta = zeros(nr,ntheta,N);
% 
% for n = 1:N
%     num = (n-1)*stride + nstart;
%     filename = strcat(dir, var, '/', var, '_', num2str(num,'%08.f'), '_', loc, '_', 'uniform_pchip.res');
%     disp(filename);
%     fid = fopen(filename);
%     h = fread(fid,0,'*uint64'); % May need adjusting
%     a = fread(fid, nr*ntheta, '*double');
%     fclose(fid);
%    for j = 1:ntheta
%         for i = 1:nr
%             vel_utheta(i,j,n) = a((j-1)*nr + i, 1);
%         end
%     end
%     
% end

%% Subtracting the temporal mean from the field

vel_mean = mean(vel_ux,3);
vel_ux = vel_ux - vel_mean;
nr_tot = 356; ntheta_tot = 258;
% Calculating the TKE in at each location

tke_ux = 0;
for i = 1:N
    disp(i);
    tke_ux = tke_ux + vel_ux(:,:,i) .*vel_ux(:,:,i);
end

tke_ux = tke_ux/N;
tke_ux = tke_ux(2:nr_tot-1, 2:ntheta_tot-1);

% vel_ur = vel_ur - mean(vel_ur,3);
% nr_tot = 356; ntheta_tot = 258;
% Calculating the TKE in at each location
% 
% tke_ur = 0;
% for i = 1:N
%     disp(i);
%     tke_ur = tke_ur + vel_ur(:,:,i) .*vel_ur(:,:,i);
% end
% 
% tke_ur = tke_ur/N;
% tke_ur = tke_ur(2:nr_tot-1, 2:ntheta_tot-1);
% 
% vel_utheta = vel_utheta - mean(vel_utheta,3);
% nr_tot = 356; ntheta_tot = 258;
% % Calculating the TKE in at each location
% 
% tke_utheta = 0;
% for i = 1:N
%     disp(i);
%     tke_utheta = tke_utheta + vel_utheta(:,:,i) .*vel_utheta(:,:,i);
% end
% 
% tke_utheta = tke_utheta/N;
% tke_utheta = tke_utheta(2:nr_tot-1, 2:ntheta_tot-1);

%% Calculating the weights

numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end
 
nr = size(rc,1);
ntheta = 256;

% Weights of radial direction
nothetar = length(rc);
weight_thetar = zeros(nr,1);
weight_thetar(1) = pi*( rc(1) + (rc(2)-rc(1))/2)^2 - pi*(rc(1))^2;

for i=2:nothetar-1
    weight_thetar(i) = pi*( rc(i) + (rc(i+1)-rc(i))/2 )^2 - pi*( rc(i) - (rc(i)-rc(i-1))/2 )^2;
end

weight_thetar(nothetar) = pi*rc(end)^2 - pi*( rc(end) - (rc(end)-rc(end-1))/2 )^2;

% Weights in azimuthal direction
weight_theta = (2*pi/ntheta)*ones(ntheta,1);   % Check once again SHEEL NIDHAN

if numvar == 4
    weight_rtheta = weight_thetar*weight_theta';
    weight_rtheta_column = weight_rtheta(:);
elseif numvar == 3
    weight_rtheta = weight_thetar;
    weight_rtheta_column = weight_rtheta(:);
end


%% Integrating over the whole area 

total_integrated_tke = (2*pi/256)*trapz(trapz(rc,rc.*tke_ux,1)); %+ (2*pi/256)*trapz(trapz(rc,tke_ur,1)) + (2*pi/256)*trapz(trapz(rc,tke_utheta,1));
%total_integrated_tke = (2*pi/256)*trapz(trapz(rc(1:187,1),tke_ux(1:187,:),1)) + (2*pi/256)*trapz(trapz(rc(1:187,1),tke_ur(1:187,:),1)) ...
%    + (2*pi/256)*trapz(trapz(rc(1:187,1),tke_utheta(1:187,:),1));

save('tke_ux_x_D_50.mat');
