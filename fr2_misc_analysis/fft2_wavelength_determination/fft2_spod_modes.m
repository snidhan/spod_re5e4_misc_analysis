function [max_kx, max_ky] = fft2_spod_modes(x_D, f_idx)

%% Name - Sheel Nidhan
%  Date - 20th February 2020
%  Code to determine the most dominant angle of emanation of IWs from the SPOD mode using FFT2
%% SPOD Parameters

close all;

Nfreq = 512;
Novlp = 256;
N     = 7000;
stride = 100;
nstart = 2329600;
nend = nstart + (N-1)*stride;
nr = 333;
ntheta = 256;
numvar = 4;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nblk_sampled = 3;
Nfreq_sampled = 25;
Nrows = numvar*nr*ntheta;

filename = strcat('/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_3.0/matlab_files/eigenmodes/eigenmodes_x_D_', int2str(x_D), '.mat');
%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*AGROW,*SAGROW> % Centered the grid faces to grid centers
end

rc = rc(1:333);

theta = linspace(0,2*pi,ntheta)';
%% dt for calculating frequency

dt = 0.0905441280000332;

%% Fixing the frequency of SPOD spectrum

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

f = f';

%% Reading SPOD eigenvalues files

load(filename); %#ok<*LOAD>

%% FFT of leading eigenmode

% fft_mode = real(eigmodes(:,:,4,1,6));
fft_mode = eigmodes(:,:,4,1,f_idx); %#ok<*IDISVAR,*NODEF>

[C,h,x,y] = polarcont(rc, theta, real(fft_mode), 10); %#ok<*ASGLU>
axis equal;
set(h,'Linecolor','none');
set(h,'LevelList',-1:0.01:1)
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20); 
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);

%% Map the eigenmode on a cartesian coordinate using polarcont function

% truncated_domain_radius = 7;
% 
% [val nr_truncated] = min(abs(truncated_domain_radius-rc));
% 
% x_truncated = zeros(nr_truncated,ntheta);
% y_truncated = zeros(nr_truncated,ntheta);
% 
% for j = 1:nr_truncated
%     for k = 1:ntheta
%         x_truncated(j,k) = rc(j)*cos(theta(k));
%         y_truncated(j,k) = rc(j)*sin(theta(k));
%     end
% end
% 
% [C,h] = contourf(x_truncated,y_truncated,fft_mode(1:nr_truncated,:),10);
% axis equal;
% set(h,'Linecolor','none');
% set(h,'LevelList',-1:0.01:1)
% colorbar;

%% Using scattered interpolant to interpolate data

[X_cart, Y_cart] = meshgrid(linspace(-7,7,600));

[T, R] = meshgrid(theta,rc);
[X, Y] = pol2cart(T, R);

S = scatteredInterpolant(X(:), Y(:), fft_mode(:));
interpolated_fft_mode = S(X_cart, Y_cart);

figure;
% interpolated_fft_mode = interp2(X,Y,fft_mode,X_cart,Y_cart);
[C,h] = contourf(X_cart, Y_cart, real(interpolated_fft_mode), 10);
axis equal;
set(h,'Linecolor','none');
set(h,'LevelList',-1:0.01:1)

hold on;

x = linspace(-6,6,100)';
y = tand(85.91)*(x+5.387) + 3.307;
xlim([-7 7]);
ylim([-7 7]);
plot(x,y,'k-','linewidth',2)
%% FFT2 on the uniform data set 

% interpolated_fft_mode_padded = padarray(interpolated_fft_mode, [2 2], 0, 'both');

interpolated_fft_mode(1:2,:) = 0;
interpolated_fft_mode(end:end-1,:) = 0;

interpolated_fft_mode(:,1:2) = 0;
interpolated_fft_mode(:,end:end-1) = 0;

FFTed_mode = fft2(interpolated_fft_mode);
abs_FFTed_mode = abs(FFTed_mode);

%% FFT parameters 

L = size(FFTed_mode,1);

sampling_frequency = 1/((X_cart(1,end) - X_cart(end,1))/L);
sampling_period = 1/sampling_frequency;
kx = sampling_frequency*(-L/2+1:L/2)/L;

sampling_frequency = 1/((Y_cart(end,1) - Y_cart(1,end))/L);
sampling_period = 1/sampling_frequency;
ky = sampling_frequency*(-L/2+1:L/2)/L;

[KX KY] = meshgrid(kx, ky);

%% Manipulating the FFT mode calculation

abs_FFTed_mode_rearranged(:,1:L/2-1)   = abs_FFTed_mode(:,L/2+2:L);
abs_FFTed_mode_rearranged(:,L/2:L)     = abs_FFTed_mode(:,1:L/2+1);

abs_FFTed_mode_rearranged_2(1:L/2-1,:) = abs_FFTed_mode_rearranged(L/2+2:L,:);
abs_FFTed_mode_rearranged_2(L/2:L,:)   = abs_FFTed_mode_rearranged(1:L/2+1,:);

%% Save the wavevector in the 
[val idx] = max(abs_FFTed_mode_rearranged_2(:));
[idx1 idx2] = find(abs_FFTed_mode_rearranged_2 == val);

max_kx = KX(idx1,idx2);
max_ky = KY(idx1,idx2);
%% FFT parameters 

figure;
[C,h] = contourf(2*pi*KX, 2*pi*KY, abs_FFTed_mode_rearranged_2, 10);
xlim([-12 12]);
ylim([-12 12]);
% axis equal;
set(h,'Linecolor','none');
% set(h,'LevelList',-1:0.01:1)
colorbar
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('fft2_eigenmode_','x_D_',int2str(x_D), '_fidx_', int2str(f_idx), '.png'),'-dpng','-r600');  

end