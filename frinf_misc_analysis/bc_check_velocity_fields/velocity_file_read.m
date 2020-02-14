%% Sheel Nidhan
%  Date - 13th February 2020

clear; clc; close all;
%% Read in the velocity field at one time instant

x = 40;
% dir_vel   = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
dir_vel = './';

%% Read the grid file

nr  = 354;
ntheta = 256;

numvar = 3;

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

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

%% Read the radial velocity field

nr = 356; ntheta = 258;

var1 = 'up';
dir = dir_vel;
u = zeros(nr,ntheta,1);
num = 2600000;
filename = strcat(dir, var1, '_', num2str(num,'%08.f'), '_', int2str(x), '.res');
disp(filename);
fid = fopen(filename);
h = fread(fid,0,'*uint64'); % May need adjusting
a = fread(fid, nr*ntheta, '*double');
fclose(fid);

for j = 1:ntheta
    for i = 1:nr
        u(i,j,1) = a((j-1)*nr + i, 1);
    end
end

nr = 356; ntheta = 258;

%%
var1 = 'vp';
dir = dir_vel;
v = zeros(nr,ntheta,1);
num = 2600000;
filename = strcat(dir, var1, '_', num2str(num,'%08.f'), '_', int2str(x), '.res');
disp(filename);
fid = fopen(filename);
h = fread(fid,0,'*uint64'); % May need adjusting
a = fread(fid, nr*ntheta, '*double');
fclose(fid);

for j = 1:ntheta
    for i = 1:nr
        v(i,j,1) = a((j-1)*nr + i, 1);
    end
end

%%
var1 = 'wp';
dir = dir_vel;
w = zeros(nr,ntheta,1);
num = 2600000;
filename = strcat(dir, var1, '_', num2str(num,'%08.f'), '_', int2str(x), '.res');
disp(filename);
fid = fopen(filename);
h = fread(fid,0,'*uint64'); % May need adjusting
a = fread(fid, nr*ntheta, '*double');
fclose(fid);

for j = 1:ntheta
    for i = 1:nr
        w(i,j,1) = a((j-1)*nr + i, 1);
    end
end

%% Plotting the value at centerline for all the velocity fields

figure;
plot(u(1,:), 'k-', 'Linewidth',2);

hXLabel = xlabel('$N_{\theta}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$u_{r}$','interpreter','latex','fontsize',15);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'ur_axis_time_instant_2600000_complete.png','-dpng','-r300');  


figure;
plot(v(1,:), 'k-', 'Linewidth',2);

hXLabel = xlabel('$N_{\theta}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$u_{\theta}$','interpreter','latex','fontsize',15);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'utheta_axis_time_instant_2600000_complete.png','-dpng','-r300');  

figure;
plot(w(1,:), 'k-', 'Linewidth',2);

hXLabel = xlabel('$N_{\theta}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$u_{x}$','interpreter','latex','fontsize',15);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'ux_axis_time_instant_2600000_complete.png','-dpng','-r300');  