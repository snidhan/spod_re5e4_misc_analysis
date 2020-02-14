%% Written by Sheel Nidhan
% Used to project the modes on actual flow and to find the flow instant matching with the mode
clear; clc; close all;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
numvar = 3;
x = 10;
dir_modes = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_', int2str(x), '/eigenmodes/');
dir_spectrum = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/matlab_files/spectrum/';
dir_vel   = strcat('/home/sheel/Work2/projects_data/spod_re5e4/frinf/data_files_uniform/x_D_', int2str(x), '/');


%% Reading the time file
        
time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);

%% Reading SPOD modes files

Nblk = floor((N-Novlp)/(Nfreq-Novlp));
nr  = 354;
numvar = 3;
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;
Nf_sampled = 50;
Nblk_sampled = 3;
mode = [0; 1; 2];
eigmodes = zeros(Nrows, Nf_sampled, size(mode,1));

for i = 1:Nf_sampled
    for j = 1:size(mode,1)
        freq = sprintf('%04d',i);   
        m = sprintf('%03d',mode(j,1));        
        filename = strcat(dir_modes, 'eigenmode_freq_',freq,'_',m,'.mod');
        disp(filename);
        A = importdata(filename);
        eigmodes(:,i,j) = A(:,1) + sqrt(-1).*A(:,2);
    end
end

save ('eigmodes_m012_nf50_x_D_10.mat', 'eigmodes', 'mode', 'Nf_sampled');


%load('eigmodes_m012_nf512_x_D_80.mat');

%% Rearranging the eigmodes for visualization

for i = 1:Nf_sampled
    disp(i);
    for j = 1:size(mode,1)
        eigenmodes_arranged(:,:,i,j) = reshape(eigmodes(:,i,j), [Nrows_permode, Nblk]); %#ok<*SAGROW,*AGROW>
    end
end
eigenmodes_proper_ranked = flip(eigenmodes_arranged, 2);
clear eigenmodes_arranged;
disp('Ranked the eigenmodes');

for i = 1:Nf_sampled
    disp(i);
    for j = 1:Nblk_sampled
        for k = 1:size(mode,1)
            eigenmodes_separated(:,:,j,i,k) = reshape(eigenmodes_proper_ranked(:,j,i,k), [nr, numvar]);
        end
    end
end

disp('Separated the eigenmodes')
clear eigenmodes_proper_ranked;
clear eigmodes;

u_eigenmode = squeeze(eigenmodes_separated(:,1,:,:,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:,:,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:,:,:));

%% Loading the grid file in radial direction

nr  = 354;
ntheta = 256;

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

%% Loading the matlabfiles for SPOD modes

mode_sampled = [0 ;1; 2];
eigenspectra_allm = zeros(Nf_sampled, Nblk_sampled, size(mode_sampled,1));
for i_mode = 1:size(mode_sampled,1)
        filename = strcat(dir_spectrum, 'spectrum_x_D_', int2str(x), '.mat');
        disp(filename);
        load(filename);
        eigenspectra_allm(:,:,i_mode) = eigvalue(1:Nf_sampled,1:Nblk_sampled,i_mode);
end

%% Renormalizing the eigenmodes by multiplying with their respective sqrt(lambda)

for i_mode = 1:size(mode_sampled,1)
    for f_mode = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
                u_eigenmode(:,Nb,f_mode,i_mode) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode))*u_eigenmode(:,Nb,f_mode,i_mode);
                v_eigenmode(:,Nb,f_mode,i_mode) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode))*v_eigenmode(:,Nb,f_mode,i_mode);
                w_eigenmode(:,Nb,f_mode,i_mode) = sqrt(eigenspectra_allm(f_mode,Nb,i_mode))*w_eigenmode(:,Nb,f_mode,i_mode);
        end
    end
end
%% Normalizing eigenmodes as per Moin and Moser eqn. 3.3

mode  = [0; 1; 2];
for m = 1:size(mode,1)
    for fn = 1:Nf_sampled
        for Nb = 1:Nblk_sampled
             spod_mode_u = u_eigenmode(:,Nb,fn,m);
             spod_mode_v = v_eigenmode(:,Nb,fn,m);
             spod_mode_w = w_eigenmode(:,Nb,fn,m);
             spod_mode_mag = spod_mode_u.*conj(spod_mode_u) + spod_mode_v.*conj(spod_mode_v) + spod_mode_w.*conj(spod_mode_w);
             alpha = dot(weight_thetar,spod_mode_mag);
             u_eigenmode(:,Nb,fn,m) = spod_mode_u/sqrt(alpha);
             v_eigenmode(:,Nb,fn,m) = spod_mode_v/sqrt(alpha);
             w_eigenmode(:,Nb,fn,m) = spod_mode_w/sqrt(alpha);
        end
    end
end

%% Defining c
Nf_sampled_scalogram = 50;
c = zeros(1, size(mode,1), Nf_sampled_scalogram, N);

%% Read u velocity field 

nr = 356; ntheta = 258;

var1 = 'up';
dir = dir_vel;
u = zeros(nr,ntheta,N);
for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var1, '/', var1, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
    disp(filename);
    fid = fopen(filename);
    h = fread(fid,0,'*uint64'); % May need adjusting
    a = fread(fid, nr*ntheta, '*double');
    fclose(fid);
   for j = 1:ntheta
        for i = 1:nr
            u(i,j,n) = a((j-1)*nr + i, 1);
        end
   end
end

disp('centering the velocities');
u_centered = 0.5*(u(2:nr-1,2:ntheta-1,:) + u(1:nr-2,2:ntheta-1,:));
u_mean = squeeze(mean(u_centered,3));
u_fluc = u_centered - u_mean;
u_fluc_sampled = u_fluc(:,:,1:end); %#ok<*NASGU>
clear u_fluc; clear u_centered; clear u;

% Fourier decomposition of the velocity fields in azimuthal

nr = 354;
ntheta = 256;
for i = 1:N
    for j = 1:nr
            u_fluc_fft(j,:,i) = (1/ntheta)*fft(u_fluc_sampled(j,:,i)); %% Azimuthal FFT at each radial location
    end
end

u_fluc_fft_m012 = u_fluc_fft(:,1:3,:);

clear u_fluc_sampled u_fluc_fft;
% Calculating the coefficients 


for i = 1:1
    for j = 1:size(mode,1)
        disp(j);
        for k = 1:Nf_sampled_scalogram
            for l = 1:N
            product_u = u_fluc_fft_m012(:,j,l).*conj(u_eigenmode(:,i,k,j));
            %product_v = v_fluc_sampled(:,j,k).*conj(v_eigenmode(:,i,k,j));
            %product_w = w_fluc_sampled(:,j,k).*conj(w_eigenmode(:,i,k,j));
            total_product = product_u;
            c(i,j,k,l) = dot(weight_thetar, total_product) + c(i,j,k,l);
            end
        end
    end
end

clear u_fluc_fft_m012;

%% Read v velocity field 

nr = 356; ntheta = 258;

var2 = 'vp';
dir = dir_vel;

v = zeros(nr,ntheta,N);

for n = 1:N
    num = (n-1)*stride + nstart;
    filename = strcat(dir, var2, '/', var2, '_', num2str(num,'%08.f'), '_', int2str(x), '_', 'uniform_pchip.res');
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

disp('centering the velocities');
v_centered = 0.5*(v(2:nr-1,2:ntheta-1,:) + v(2:nr-1,1:ntheta-2,:));
v_mean = squeeze(mean(v_centered,3));
v_fluc = v_centered - v_mean;
v_fluc_sampled = v_fluc(:,:,1:end);
clear v_fluc; clear v_centered; clear v;

% Fourier decomposition of the velocity fields in azimuthal

nr = 354;
ntheta = 256;
for i = 1:N
    for j = 1:nr
            v_fluc_fft(j,:,i) = (1/ntheta)*fft(v_fluc_sampled(j,:,i)); %% Azimuthal FFT at each radial location
    end
end

v_fluc_fft_m012 = v_fluc_fft(:,1:3,:);

clear v_fluc_sampled v_fluc_fft;
% Calculating the coefficients 


for i = 1:1
    for j = 1:size(mode,1)
        disp(j);
        for k = 1:Nf_sampled_scalogram
            for l = 1:N
            product_v = v_fluc_fft_m012(:,j,l).*conj(v_eigenmode(:,i,k,j));
            %product_v = v_fluc_sampled(:,j,k).*conj(v_eigenmode(:,i,k,j));
            %product_w = w_fluc_sampled(:,j,k).*conj(w_eigenmode(:,i,k,j));
            total_product = product_v;
            c(i,j,k,l) = dot(weight_thetar, total_product) + c(i,j,k,l);
            end
        end
    end
end

clear v_fluc_fft_m012;

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

% Fourier decomposition of the velocity fields in azimuthal

nr = 354;
ntheta = 256;
for i = 1:N
    for j = 1:nr
            w_fluc_fft(j,:,i) = (1/ntheta)*fft(w_fluc_sampled(j,:,i)); %% Azimuthal FFT at each radial location
    end
end

w_fluc_fft_m012 = w_fluc_fft(:,1:3,:);

clear w_fluc_sampled w_fluc_fft;
% Calculating the coefficients 


for i = 1:1
    for j = 1:size(mode,1)
        disp(j);
        for k = 1:Nf_sampled_scalogram
            for l = 1:N
            product_w = w_fluc_fft_m012(:,j,l).*conj(w_eigenmode(:,i,k,j));
            %product_v = v_fluc_sampled(:,j,k).*conj(v_eigenmode(:,i,k,j));
            %product_w = w_fluc_sampled(:,j,k).*conj(w_eigenmode(:,i,k,j));
            total_product = product_w;
            c(i,j,k,l) = dot(weight_thetar, total_product) + c(i,j,k,l);
            end
        end
    end
end

clear w_fluc_fft_m012;
%% Saving coefficients file
save('coefficients_projection_x_D_10.mat', 'c', 'Nf_sampled', 'f');
