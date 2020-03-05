%% Written by Sheel Nidhan
%  Calculates the ratio of energy residing between TKE and TPE for different modes

clear; clc; close all;
%% SPOD Parameters
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

%x_sampled =  [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
x_sampled = [20; 30; 40; 50; 60; 70; 80; 90; 100]; %#ok<*NBRAK>
% loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];
loc_planes = [20; 30; 40; 50; 60; 70; 80; 90; 100];

dir_modes = '/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_3.0/matlab_files/eigenmodes/';
dir_spec  = '/home/sheel/Work2/projects_data/spod_re5e4/fr2/spod_data/run_3.0/matlab_files/';

disp(dir_modes);
disp(dir_spec);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*AGROW,*SAGROW> % Centered the grid faces to grid centers
end

rc = rc(1:333);

theta = linspace(0,2*pi,ntheta)';

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
%     weight_rtheta = weight_thetar*weight_theta';            % When implementing serial spod version
    weight_rtheta = weight_thetar;                            % When implementing parallel spod version
    weight_rtheta_column = weight_rtheta(:)/(2*pi);
elseif numvar == 3
    weight_rtheta = weight_thetar;
    weight_rtheta_column = weight_rtheta(:);
end
weight_theta = repmat(weight_rtheta_column,1,ntheta);

%% Reading the time file

dt = 0.0905441280000332;

%%  Partitioning the time data in different blocks

Nblk = floor((N-Novlp)/(Nfreq-Novlp));

qstart = zeros(Nblk,1);
qend  = zeros(Nblk,1);

for i = 1:Nblk
    qstart(i,1) = (i-1)*(Nfreq-Novlp) + 1;
    qend(i,1) = qstart(i) + Nfreq - 1;
end

%% Fixing the frequency of SPOD spectrum

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

%% Loading the matlabfiles for SPOD modes

u_eigenmode_allm   = zeros(nr, ntheta, Nblk_sampled, Nfreq_sampled, size(x_sampled,1));
w_eigenmode_allm   = zeros(nr, ntheta, Nblk_sampled, Nfreq_sampled, size(x_sampled,1));
v_eigenmode_allm   = zeros(nr, ntheta, Nblk_sampled, Nfreq_sampled, size(x_sampled,1));
rho_eigenmode_allm = zeros(nr, ntheta, Nblk_sampled, Nfreq_sampled, size(x_sampled,1));

for x_d = 1:size(x_sampled,1)
        filename = strcat(dir_modes, 'eigenmodes_x_D_', int2str(x_sampled(x_d,1)), '.mat');
        disp(filename);
        load(filename);
        u_eigenmode_allm(:,:,:,:,x_d)   = eigmodes(:,:,1,:,:);
        v_eigenmode_allm(:,:,:,:,x_d)   = eigmodes(:,:,2,:,:);
        w_eigenmode_allm(:,:,:,:,x_d)   = eigmodes(:,:,3,:,:);
        rho_eigenmode_allm(:,:,:,:,x_d) = eigmodes(:,:,4,:,:);

end

%% Loading the eigenvalues for SPOD modes

eigenspectra_allm = zeros(Nfreq_sampled, Nblk_sampled, size(x_sampled,1));
filename = strcat(dir_spec,'fr2_eigvalue_spectrum_slice.mat');
load(filename);

for x_d = 1:size(x_sampled,1)
        for j = 1:size(eigvalues_spectrum_plot,1)
        if (eigvalues_spectrum_plot(j).loc == x_sampled(x_d,1))
            eigenspectra_allm(:,:,x_d) = eigvalues_spectrum_plot(x_d).eigvalue(1:Nfreq_sampled, 1:Nblk_sampled);
        end
        end
end
%% Renormalizing the eigenmodes by multiplying with their respective sqrt(lambda)

for x_d = 1:size(x_sampled,1)
    for f_mode = 1:Nfreq_sampled
        for Nb = 1:Nblk_sampled
            u_eigenmode_allm(:,:,Nb,f_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,x_d))*u_eigenmode_allm(:,:,Nb,f_mode,x_d);
            v_eigenmode_allm(:,:,Nb,f_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,x_d))*v_eigenmode_allm(:,:,Nb,f_mode,x_d);
            w_eigenmode_allm(:,:,Nb,f_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,x_d))*w_eigenmode_allm(:,:,Nb,f_mode,x_d);
            rho_eigenmode_allm(:,:,Nb,f_mode,x_d) = sqrt(eigenspectra_allm(f_mode,Nb,x_d))*rho_eigenmode_allm(:,:,Nb,f_mode,x_d);
        end
    end
end
%% Calculating the ratio of TKE to TPE

energy_partition = zeros(size(x_sampled,1), Nfreq_sampled, Nblk_sampled);
energy_rho       = zeros(size(x_sampled,1), Nfreq_sampled, Nblk_sampled);

for x = 1:size(x_sampled,1)
    for fn = 1:Nfreq_sampled
        for Nb = 1:Nblk_sampled
             spod_mode_u =   u_eigenmode_allm(:,:,Nb,fn,x);
             spod_mode_v =   v_eigenmode_allm(:,:,Nb,fn,x);
             spod_mode_w =   w_eigenmode_allm(:,:,Nb,fn,x);
             spod_mode_rho = rho_eigenmode_allm(:,:,Nb,fn,x);
%             full_spod_mode = [spod_mode_u; spod_mode_v; spod_mode_w];
%             spod_mode_u = spod_mode_u/norm(full_spod_mode);
%             spod_mode_v = spod_mode_v/norm(full_spod_mode);
%             spod_mode_w = spod_mode_w/norm(full_spod_mode);
             spod_mode_mag_tke = spod_mode_u.*conj(spod_mode_u) + spod_mode_v.*conj(spod_mode_v) + spod_mode_w.*conj(spod_mode_w);
             spod_mode_mag_rho = spod_mode_rho.*conj(spod_mode_rho);
             alpha_tke = sum(sum((spod_mode_mag_tke.*weight_theta)));
             alpha_rho = sum(sum((spod_mode_mag_rho.*weight_theta)));
             energy_partition(x,fn,Nb) = alpha_rho/(alpha_tke+alpha_rho);
             energy_rho(x,fn,Nb)       = alpha_rho;
        end
    end
end

%% Saving the mat file for future post-processing 

%
%
%
%
%% Plotting the ratio of TPE/TKE 

figure;
hold on;
% h1 = plot(x_sampled, energy_partition(:,1,1), 'ro', 'Markersize', 10);
% h2 = plot(x_sampled, energy_partition(:,2,1), 'ko', 'Markersize', 10);
% h3 = plot(x_sampled, energy_partition(:,3,1), 'bo', 'Markersize', 10);
h4 = plot(x_sampled, energy_partition(:,21,1), 'ko', 'Markersize', 10);
