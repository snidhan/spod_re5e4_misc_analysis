%% Written by Sheel Nidhan
%  Date - 24th February, 2020
clear; clc;
close all;
%% Reading the datafiles of u_x velocity field

var = 'wp';
loc = '80';
nstart = 1892600;
nend   = 2613200;
stride  = 100;
dir = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_80/';
nr = 356;
ntheta = 258;
N = (nend-nstart)/stride + 1;
vel = zeros(nr,ntheta,N);

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
            vel(i,j,n) = a((j-1)*nr + i, 1);
        end
    end
    
end
%% Getting fluctuations from the time series data

vel_centered = vel(2:nr-1,2:ntheta-1,:);
disp('centered the velocities');
clear vel;

%% Getting the mean

vel_mean = squeeze(mean(vel_centered,3));
disp('removed the mean velocities');

%% Fluctuations of radial, streamwise and azimuthal velocity

vel_fluc = vel_centered - vel_mean;

%% Performing the azimuthal FFT of the velocity field 

nr_truncated = nr-2; ntheta_truncated = ntheta-2;
% vel_fluc = vel_fluc(2:nr-1, 2:ntheta-1,:); %% Removing the ghost cells

for i = 1:N
    for j = 1:nr_truncated
            vel_fluc(j,:,i) = (1/ntheta_truncated)*fft(vel_fluc(j,:,i)); %% Azimuthal FFT at each radial location
%         vel(j,:,i) = fft(vel(j,:,i)); %% Azimuthal FFT at each radial location
 
    end
end

%% Loading the grid file in radial direction

nr  = 354;
ntheta = 256;

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*SAGROW> % Centered the grid faces to grid centers
end

%% Contour plot of velocity fluctuations as a function of r/D and t

nt = 7200;
nr_sampled = 300;
vel_m1 = squeeze(vel_fluc(1:nr_sampled,2,1:nt));
vel_m2 = squeeze(vel_fluc(1:nr_sampled,3,1:nt));

%% Save .mat file for plots

save('x_D_80_fft.mat', 'vel_m1', 'vel_m2', 'nr_sampled', 'nt', 'r');

%% Plotting the contour of fluctuations for m = 1

close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% close all;
% figure
dt = 0.07;
t = dt*linspace(0,nt-1,nt)';
[R T] = meshgrid(rc(1:nr_sampled), t);

% [C,h] = contourf(T,R,real(vel_m1)', 20);
% colormap jet;
% set(h, 'edgecolor', 'none')
% colorbar;
% caxis([-0.05, 0.05]);
% hold on;
% hline = plot(t,0.8*ones(size(t,1),1),'k--', 'Linewidth',3);
%% Setting the figure environment and saving the contours

% x0=0;
% y0=0;
% width=30;
% height=10;
% set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
% 
% ax = gca;
% ax.FontSize = 15;
% 
% hXLabel = xlabel('$tU_{\infty}/D$','interpreter','latex','fontsize',20); %#ok<*NASGU>
% hYLabel = ylabel('$r/D$','interpreter','latex','fontsize',20); %#ok<*NASGU>
% hTitle  = title('$Re(\tilde{u''}_{x}(m=1,r,t))$','interpreter','latex','fontsize', 20);
% 
% dirout = './';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_fluc_x_D_', loc, '_m1.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_fluc_x_D_', int2str(loc), '_m1.eps'),'-depsc2','-r600');

%% FFT of m = 1 mode to uncover the vortex dynamics

dt = 0.0721;
data = (vel_m1)';
Fs   = 1/dt;
L = size(t,1);
W = hamming(L);

% Plotting data at r/D \approx 0.8;
% h_fft_time_series = plot(t, data(:,60), 'k-', 'Linewidth',2);

FFT_time_series = fft(W.*data(:,60));

P2 = abs(FFT_time_series/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

close all;
loglog(f,P1, 'k-', 'Linewidth',2); 
hold on;
loglog(f(500:end), 10^-3*f(500:end).^(-5/3), 'k--', 'Linewidth',2);
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');
%% Plotting the contour of fluctuations for m = 2

close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% close all;
% figure
dt = 0.07;
t = dt*linspace(0,nt-1,nt)';
[R T] = meshgrid(rc(1:nr_sampled), t);

% [C,h] = contourf(T,R,real(vel_m2)', 20);
% colormap jet;
% set(h, 'edgecolor', 'none')
% colorbar;
% caxis([-0.05, 0.05]);

%% Setting the figure environment

% x0=0;
% y0=0;
% width=30;
% height=10;
% set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
% 
% ax = gca;
% ax.FontSize = 15;
% 
% hXLabel = xlabel('$tU_{\infty}/D$','interpreter','latex','fontsize',20); %#ok<*NASGU>
% hYLabel = ylabel('$r/D$','interpreter','latex','fontsize',20); %#ok<*NASGU>
% hTitle  = title('$Re(\tilde{u''}_{x}(m=2,r,t))$','interpreter','latex','fontsize', 20);
% 
% % Saving the contours
% 
% dirout = './';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'ux_fluc_x_D_', loc, '_m2.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'ux_fluc_x_D_', int2str(loc), '_m2.eps'),'-depsc2','-r600');

%% FFT of m = 1 mode to uncover the vortex dynamics

dt = 0.0721;
data = (vel_m2)';
Fs   = 1/dt;
L = size(t,1);
W = hamming(L);

% Plotting data at r/D \approx 0.8;
% h_fft_time_series = plot(t, data(:,60), 'k-', 'Linewidth',2);

FFT_time_series = fft(W.*data(:,200));

P2 = abs(FFT_time_series/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

close all;
loglog(f,P1, 'k-', 'Linewidth',2); 
hold on;
loglog(f(500:end), 10^-3*f(500:end).^(-5/3), 'k--', 'Linewidth',2);
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');