%% Written by Sheel Nidhan
%  This code is used to compare the spod_fortran code with MATLAB implementation

clear; clc;
% Reading the datafiles of u_x velocity field

var = 'up';
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


%% Performing the azimuthal FFT of the velocity field 

nr_truncated = nr-2; ntheta_truncated = ntheta-2;
vel = vel(2:nr-1, 2:ntheta-1,:); %% Removing the ghost cells

for i = 1:N
    for j = 1:nr_truncated
            vel(j,:,i) = (1/ntheta_truncated)*fft(vel(j,:,i)); %% Azimuthal FFT at each radial location
%         vel(j,:,i) = fft(vel(j,:,i)); %% Azimuthal FFT at each radial location
 
    end
end

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
%weight_thetar(1) = ( rc(1) + (rc(2)-rc(1))/2)^2 - (rc(1))^2;

for i=2:nothetar-1
    weight_thetar(i) = pi*( rc(i) + (rc(i+1)-rc(i))/2 )^2 - pi*( rc(i) - (rc(i)-rc(i-1))/2 )^2;
    %weight_thetar(i) = ( rc(i) + (rc(i+1)-rc(i))/2 )^2 - ( rc(i) - (rc(i)-rc(i-1))/2 )^2;
end

weight_thetar(nothetar) = pi*rc(end)^2 - pi*( rc(end) - (rc(end)-rc(end-1))/2 )^2;
%weight_thetar(nothetar) = rc(end)^2 - ( rc(end) - (rc(end)-rc(end-1))/2 )^2;

% Weights in azimuthal direction
weight_theta = (2*pi/ntheta)*ones(ntheta,1);   % Check once again SHEEL NIDHAN

if numvar == 4
    weight_rtheta = weight_thetar*weight_theta';
    weight_rtheta_column = weight_rtheta(:);
elseif numvar == 3
    weight_rtheta = weight_thetar;
    weight_rtheta_column = weight_rtheta(:);
end

%% Performing the steps of SPOD for each azimuthal wavenumber m

Nfreq = 512;
Novlp = 256;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));

% Creating the hamming window manually

window = zeros(Nfreq,1);
for i = 1:Nfreq
  window(i,1) = 0.54 - 0.46*cos((2*pi*(i-1))/(Nfreq-1));
%   window(i,1) = 0.5*(1-cos(2*pi*(i-1)/(Nfreq-1)));
%  window(i,1) = 1;
end

%winweight = Nfreq/(sum(window));   % Amplitude correction in windowing
winweight = 1.587845072201937; % Hamming window energy correction
%winweight = 1.63; % Hanning window energy correction

qstart = zeros(Nblk,1); qend = zeros(Nblk,1);

for i = 1:Nblk
    qstart(i,1) = (i-1)*(Nfreq-Novlp) + 1;
    qend(i,1)   = qstart(i,1) + Nfreq - 1;
end


vel_blk = zeros(nr_truncated,Nfreq,Nblk);
vel_fft = zeros(nr_truncated,Nfreq,Nblk);
vel_k   = zeros(nr_truncated,Nblk,Nfreq);
spod_eigvalues  = zeros(Nblk, Nfreq, ntheta_truncated);

for j = 1:10
    disp(j);
    vel_spod = squeeze(vel(:,j,:)); % Azimuthal mode over which FFT is performed
    vel_mean = mean(vel_spod,2);
    vel_spod = vel_spod - vel_mean;
    
    for i = 1:Nblk                  % Separating the snapshot matrix into blocks
        vel_blk(:,1:Nfreq,i) = vel_spod(:,qstart(i,1):qend(i,1));
    end
    
    for i = 1:Nblk                  % Windowing the velocity fields 
        for k = 1:Nfreq
            vel_blk(:,k,i) = vel_blk(:,k,i)*window(k,1);
        end
    end
    
    vel_blk = (winweight/Nfreq)*vel_blk;   % Normalizing the velocity fields before taking FFT
    
    for i = 1:Nblk                         % Taking FFT in time 
        for k = 1:nr_truncated
            vel_fft(k,:,i) = fft(vel_blk(k,:,i));
        end
    end
    
    for i = 1:Nblk                        % Arranging the spectral FFTed field
        for k = 1:Nfreq
            vel_k(:,i,k) = vel_fft(:,k,i);
        end
    end
 
    S = zeros(Nblk, Nblk, Nfreq);
   
    for i = 1:Nfreq                    % Calculating the cross spectral density 
        S(:,:,i) =  (squeeze(vel_k(:,:,i))')*(weight_rtheta_column.*squeeze(vel_k(:,:,i)));
    end
    
    S = S/Nblk;
    
    for i = 1:Nfreq                   % Calculating the eigenvalues of CSD matrix 
        [v, e] = eig(S(:,:,i));
        spod_eigvalues(:,i,j) = real(diag(e));
        spod_eigenmodes(:,:,i,j) = v;
    end
    
    for i = 1:Nfreq
        disp(i);
        spod_mode_scaled(:,:,j,i) = squeeze(vel_k(:,:,i))*squeeze(spod_eigenmodes(:,:,i,j));
    end
   
end


%% Plotting the eigenvalues 

time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);

time_blk = zeros(Nblk,Nfreq);
for i = 1:Nblk
    time_blk(i,:) = time_spod(qstart(i,1):qend(i,1),1)';
end

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

%%
% h1 = loglog(f, real(squeeze(spod_eigvalues(1,:,1))),'k','Linewidth',2);
% hold on
% h2 = loglog(f, real(squeeze(spod_eigvalues(1,:,2))),'r','Linewidth',2);
% h3 = loglog(f, real(squeeze(spod_eigvalues(1,:,3))),'b','Linewidth',2);