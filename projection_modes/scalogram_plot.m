clear; clc;
load('coefficients_projection_x_D_80.mat')

f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);

[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);

%%
h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');

colormap hot;
colorbar;

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;
x = 80;
dir_modes = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_', int2str(x), '/eigenmodes/');
dir_spectrum = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/spectrum/';
dir_vel   = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');


%% Read w velocity field 

nr = 356; ntheta = 258;

var3 = 'wp';
dir = dir_vel;

w = zeros(nr,ntheta,N);
for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var3, '/', var3, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
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

disp('centering the velocities');
w_centered = w(2:nr-1,2:ntheta-1,:);
w_mean = squeeze(mean(w_centered,3));
w_fluc = w_centered - w_mean;
w_fluc_sampled = w_fluc(:,:,1:end);
clear w_fluc; clear w_centered; clear w;

%% Reading the grid files

ntheta = 256;
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

figure;
%% Plotting the turbulent fluctuation


t = 111;

[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);
axis equal;
set(h,'Linecolor','none');
colormap hot;
colorbar;

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',15);

ylim([-8 8]);
xlim([-8 8]);
%% set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_80.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_80.eps'),'-depsc2','-r600');