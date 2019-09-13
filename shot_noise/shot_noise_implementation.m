%% Implementing shot-noise decomposition for characteristic eddy decomposition
%  Section 6.3 in Delville et al. 1999

clc; clear;

%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7200;
mode  = 0;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;

dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_50/eigenmodes/';
disp(dir2);
dir3 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/run_2.0/x_D_50/eigenspectra/';
disp(dir3);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);
ntheta = 256;
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

%% Reading the time file

time = importdata ('/home/sheel/Work/projects/spod_re5e4/post/frinf/time_stamps/time_stamp_1892600_2613200_uniform.txt');
time_spod = time(1:N,1);
dt   = time_spod(2:end,1) - time_spod(1:end-1,1);

%%  Partitioning the time data in different blocks

Nblk = floor((N-Novlp)/(Nfreq-Novlp));

qstart = zeros(Nblk,1);
qend  = zeros(Nblk,1);

for i = 1:Nblk
    qstart(i,1) = (i-1)*(Nfreq-Novlp) + 1;
    qend(i,1) = qstart(i) + Nfreq - 1;
end

%% Separating time data into blocks

time_blk = zeros(Nblk,Nfreq);

for i = 1:Nblk
    time_blk(i,:) = time_spod(qstart(i,1):qend(i,1),1)';
end

%% Fixing the frequency of SPOD spectrum

f = (0:Nfreq-1)/dt(1)/Nfreq;

if mod(Nfreq,2) == 0
    f(Nfreq/2 + 1:end) = f(Nfreq/2 + 1:end)-1/dt(1);
else
    f((Nfreq+1)/2 + 1:end) = f((Nfreq+1)/2 + 1:end) - 1/dt(1);
end

%% Reading first SPOD eigenmodes 

Nf_sampled = 10;
m_sampled  = 19;
mode       = linspace(0, m_sampled, m_sampled + 1)';
eigmodes = zeros(Nrows, Nf_sampled);
spod_mode_1 = zeros(Nrows_permode, m_sampled, Nf_sampled);

for j = 1:size(mode,1)
    for i = 1:Nf_sampled
        disp(i);
        freq = sprintf('%04d',i);   
        m = sprintf('%03d',mode(j,1));
        filename = strcat(dir2, 'eigenmode_freq_',freq,'_',m,'.mod');
        disp(filename);
        A = importdata(filename);
        eigmodes(:,i) = A(:,1) + sqrt(-1).*A(:,2);
        spod_mode_1(:,j,i) = eigmodes(Nrows-Nrows_permode+1:Nrows,i);
    end
end

%% Reshaping the eigenmodes 

for j = 1:size(mode,1)
    for i = 1:Nf_sampled
        spod_modes1_arranged(:,:,j,i) =  reshape(spod_mode_1(:,j,i), [nr, numvar]);
    end
end

clear spod_mode_1

%% Reading the eigenvalues of the modes of interest

eigvalue = zeros(Nfreq,Nblk,size(mode,1));

for k = 1:size(mode,1)
    for i = 1:Nfreq
        m = sprintf('%03d',mode(k,1));
        freq = sprintf('%04d',i);   
        filename = strcat(dir3, 'eigenvalues_freq_',freq,'_',m,'.txt');
        disp(filename);
        A = importdata(filename);
        eigvalue(i,:,k) = A(:,1);
    end
end

for k = 1:size(mode,1)
    for i = 1:Nfreq
        eigvalue(i,:,k) = sort(eigvalue(i,:,k), 'descend');
    end
end

eigvalue_spod_mode1 = squeeze(eigvalue(1:Nf_sampled, 1, 1:m_sampled+1));

clear eigvalue;

%% Normalizing the eigenmodes 
%  Eqn. (3.3) from Moin and Moser 1988
for m = 1:m_sampled
    for fn = 1:Nf_sampled
        for var = 1:numvar
            spod_mode = spod_modes1_arranged(:,var,m,fn);
            spod_mode_mag = spod_mode.*conj(spod_mode);
            alpha = dot(weight_thetar, spod_mode_mag);
            spod_modes1_arranged(:,var,m,fn) = spod_modes1_arranged(:,var,m,fn)/sqrt(alpha);
        end
    end
end
%% Algorithm in section 4.4 of Moin and Moser 1988

for m = 1:m_sampled
    for fn = 1:Nf_sampled
        for var = 1:numvar
            spod_modes1_arranged(:,var,m,fn) = spod_modes1_arranged(:,var,m,fn)*sqrt(eigvalue_spod_mode1(fn,m));
            spod_mode = spod_modes1_arranged(:,var,m,fn);
            % Calculate gamma
            spod_mode_gamma = spod_mode.*rc;
            gamma = trapz(rc,spod_mode_gamma);
            mult_factor = conj(gamma)/abs(gamma);
            dominant_eddy_fft(:,var,m,fn) = spod_modes1_arranged(:,var,m,fn)*mult_factor;
        end
    end
end

%% Inverse fft dominant eddy to the real space 

Ntheta = 256;
theta = linspace(0,2*pi,Ntheta)';


for i = 1:nr
    for fn = 1:Nf_sampled
        for var = 1:numvar
            dominant_eddy_physical(i,:,var,fn) = ifft(dominant_eddy_fft(i,var,1:m_sampled,fn),Ntheta); 
        end
    end
end


% for i = 1:nr
%     for j = 1:Ntheta
%         for var = 1:numvar
%             dominant_eddy(i,j,:,var) = ifft(dominant_eddy_physical(i,j,var,1:Nf_sampled),10);
%         end
%     end
% end
%% Plot 

figure;
[C,h] = polarcont(rc,theta,real(dominant_eddy_physical(:,:,1,6)),5);
set(h,'LineColor','none')
colorbar;
axis equal

