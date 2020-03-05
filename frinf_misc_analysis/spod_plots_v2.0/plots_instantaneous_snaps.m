%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting instantaneous snapshots containing structures of m=1 and m=2 

clear;
addpath('./aux_plots/');
dirout = './';

%% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;


%% Generating the plot object
% close all;
x0=0;
y0=0;
width=12;
height=8;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,3,[.05, 0.05],[.1,.15],[.1 .05]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Read w velocity field at x/D = 10

nr = 356; ntheta = 258;
  var3 = 'wp';
x = 10;
% dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
dir  = strcat('../../../../../spod_re5e4_frinf_work/frinf_rana/data_files_uniform/x_D_', int2str(x), '/');

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

save('w_fluc_sampled_x_D_10.mat', 'w_fluc_sampled', '-v7.3');

%% Reading the grid files

ntheta = 256;
theta = linspace(0,2*pi,ntheta)';
numvar = 3;   %#ok<*NASGU> % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

% fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
fid = fopen('./x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end

%% Plotting the m=2 mode at x/D = 10

load('./aux_plots/files/coefficients_projection_x_D_10.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;

maximum = max(max(abs(squeeze(c(1,3,1:50,:)))));
[x,y] = find(abs(squeeze(c(1,3,1:50,:))) == maximum);

axes(ha(1));
t = 4663;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>

set(h,'Linecolor','none');
colormap hot;
caxis([-0.3 0.3]);

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(1).Position(1) ha(1).Position(2)+ 0.36 ha(1).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.3:0.15:0.3);
set(hBar, 'TickLabelInterpreter', 'latex');

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',15);

ax = gca;
ax.FontSize = 20;

xlim([-3 3]);
xticks([-3 -2 -1 0 1 2 3]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-3 3]);
yticks([-3 -2 -1 0 1 2 3]);
yticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});

%print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_10.png'),'-dpng2','-r600');  
%print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_10.eps'),'-depsc2','-r600');
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

box on;

%% Plotting the m=1 mode at x/D = 10


load('./aux_plots/files/coefficients_projection_x_D_10.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;z
% colorbar;
maximum = max(max(abs(squeeze(c(1,2,1:50,:)))));
[x,y] = find(abs(squeeze(c(1,2,1:50,:))) == maximum);


axes(ha(4));
t = 2650;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);

set(h,'Linecolor','none');
colormap hot;
caxis([-0.2 0.2]);
%colorbar;

ax = gca;
ax.FontSize = 20;

xlim([-3 3]);
xticks([-3 -2 -1 0 1 2 3]);
xticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});
ylim([-3 3]);
yticks([-3 -2 -1 0 1 2 3]);
yticklabels({'$-3$', '$-2$', '$-1$', '$0$', '$1$', '$2$', '$3$'});

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20,  'Position', [-0.14, -3.9, 0]);
hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_10.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_10.eps'),'-depsc2','-r600');

%% Plotting the snapshots representative of m=1, 2 mode at x/D = 40;

clearvars w_mean w_fluc_sampled;
% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;

% Read w velocity field 

nr = 356; ntheta = 258;
var3 = 'wp';
x = 40;
%dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
dir  = strcat('../../../../../spod_re5e4_frinf_work/frinf_rana/data_files_uniform/x_D_', int2str(x), '/');
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

%% Plotting the m=2 mode at x/D = 40

axes(ha(2));

load('./aux_plots/files/coefficients_projection_x_D_40.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;


maximum = max(max(abs(squeeze(c(1,2,1:50,:)))));
[x,y] = find(abs(squeeze(c(1,2,1:50,:))) == maximum)

t = 4556;
% t = 1280;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>


set(h,'Linecolor','none');
colormap hot;
caxis([-0.10 0.10]);

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(2).Position(1) ha(2).Position(2)+ 0.36 ha(2).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.1:0.05:0.1);
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 20;


xlim([-6 6]);
xticks([-6 -4 -2 0 2 4 6]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-6 6]);
yticks([-6 -4 -2 0 2 4 6]);
yticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});

% hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

% print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_40.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m2st0_', int2str(t), '_x_D_40.eps'),'-depsc2','-r600');


%% Plotting the m=1 mode at x/D = 40

axes(ha(5));

load('./aux_plots/files/coefficients_projection_x_D_40.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;

t = 761;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);

set(h,'Linecolor','none');
colormap hot;
caxis([-0.10 0.10]);
%colorbar;

ax = gca;
ax.FontSize = 20;

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(2).Position(1) ha(2).Position(2)+ 0.36 ha(2).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.1:0.05:0.1);
set(hBar, 'TickLabelInterpreter', 'latex');


xlim([-6 6]);
xticks([-6 -4 -2 0 2 4 6]);
xticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});
ylim([-6 6]);
yticks([-6 -4 -2 0 2 4 6]);
yticklabels({'$-6$', '$-4$', '$-2$', '$0$', '$2$', '$4$', '$6$'});


hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20,'Position', [-0.14, -7.82, 0]);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(e)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_40.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ux_m1st0135_', int2str(t), '_x_D_40.eps'),'-depsc2','-r600');

%% Plotting the snapshots representative of m=1, 2 mode at x/D = 80;

clearvars w_mean w_fluc_sampled;
% Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;

% Read w velocity field 

nr = 356; ntheta = 258;
var3 = 'wp';
x = 80;
% dir  = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');
dir  = strcat('../../../../../spod_re5e4_frinf_work/frinf_rana/data_files_uniform/x_D_', int2str(x), '/');
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

%% Plotting the m=2 mode at x/D = 80

axes(ha(3));

load('./aux_plots/files/coefficients_projection_x_D_80.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);
% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,3,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;

t = 5133;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10); %#ok<*ASGLU>

set(h,'Linecolor','none');
colormap hot;
caxis([-0.05 0.05]);

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(3).Position(1) ha(3).Position(2)+ 0.36 ha(3).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.05:0.025:0.05);
set(hBar, 'TickLabelInterpreter', 'latex');

ax = gca;
ax.FontSize = 20;

xlim([-9 9]);
xticks([-9 -6 -3 0 3 6 9]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
ylim([-9 9]);
yticks([-9 -6 -3 0 3 6 9]);
yticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});

hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);
%% Plotting the m=1 mode at x/D = 80

axes(ha(6));

load('./aux_plots/files/coefficients_projection_x_D_80.mat')
f_sampled = f(1:50);
N_snaps   = linspace(1,7200,7200);
[N_SNAPS, F_SAMPLED] = meshgrid(f_sampled, N_snaps);

% h1 = contourf(N_SNAPS, F_SAMPLED, abs(squeeze(c(1,2,1:50,:))'),'LineStyle','none');
% colormap hot;
% colorbar;


t = 1325;
[C,h,x,y] = polarcont(rc, theta, w_fluc_sampled(:,:,t), 10);


set(h,'Linecolor','none');
colormap hot;
caxis([-0.05 0.05]);

ax = gca;
ax.FontSize = 20;

xlim([-9 9]);
xticks([-9 -6 -3 0 3 6 9]);
xticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});
ylim([-9 9]);
yticks([-9 -6 -3 0 3 6 9]);
yticklabels({'$-9$', '$-6$', '$-3$', '$0$', '$3$', '$6$', '$9$'});

hXLabel = xlabel('$z/D$','interpreter','latex','fontsize',20,'Position', [-0.14, -11.75, 0]);
% hYLabel = ylabel('$y/D$','interpreter','latex','fontsize',20);
hTitle  = title('(f)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);
% latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.14, 0.95, 0]);

%% Saving the plots
print(gcf,strcat(dirout, 'ux_m12_', '_instant_realizations.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'ux_m12_', '_instant_realizations.eps'),'-depsc2','-r600');

