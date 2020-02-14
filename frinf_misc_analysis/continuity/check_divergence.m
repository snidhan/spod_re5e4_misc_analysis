%% Uses the Taylor hypothesis to check the divergence of modes
%  Eqn. 5.11 in Delville et al. 1999

clc; clear;
%% SPOD Parameters

Nfreq = 512;
Novlp = 256;
N     = 7000;
mode  = 1;
stride = 100;
nstart = 1892600;
nend = nstart + (N-1)*stride;
nr = 354;
numvar = 3;
Nblk = floor((N-Novlp)/(Nfreq-Novlp));
Nrows = numvar*nr*Nblk;
Nrows_permode = numvar*nr;

dir2 = '/home/sheel/Work2/projects_data/spod_re5e4/frinf/spod_data/x_D_80/eigenmodes/';
disp(dir2);

%% Loading the grid file in radial direction

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
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

%% Reading SPOD eigenvalues files

eigmodes = zeros(Nrows, Nfreq/2);
m = sprintf('%03d',mode);

for i = 1:10
    
    freq = sprintf('%04d',i);   
    filename = strcat(dir2, 'eigenmode_freq_',freq,'_',m,'.mod');
    disp(filename);
    A = importdata(filename);
    eigmodes(:,i) = A(:,1) + sqrt(-1).*A(:,2);
   
end

%% Rearranging the eigmodes for visualization

for i = 1:10
    for j = 1:26
        eigenmodes_arranged(:,j,i) = eigmodes((j-1)*nr*3+1:j*nr*3,i);
    end
end

for i = 1:10
    for j = 1:Nblk
        for k = 1:3
            eigenmodes_separated(1:nr,k,j,i) = ...
                           eigenmodes_arranged(nr*(k-1)+1:nr*k,j,i);
        end
    end
end

u_eigenmode = squeeze(eigenmodes_separated(:,1,:,:));
v_eigenmode = squeeze(eigenmodes_separated(:,2,:,:));
w_eigenmode = squeeze(eigenmodes_separated(:,3,:,:));

u_eigenmode = flip(u_eigenmode, 2);   % Flip eigenmodes to match the order of eigenvalues
v_eigenmode = flip(v_eigenmode, 2);
w_eigenmode = flip(w_eigenmode, 2);

%% Normalizing the eigenmodes
% for k = 1:Nblk
%     for j = 1:Nfreq/2
%         u_normalize(:,k,j) = u_eigenmode(:,k,j)/norm(u_eigenmode(:,k,j));
%         v_normalize(:,k,j) = v_eigenmode(:,k,j)/norm(v_eigenmode(:,k,j));
%         w_normalize(:,k,j) = w_eigenmode(:,k,j)/norm(w_eigenmode(:,k,j));
%     end
% end

clear u_normalize v_normalize w_normalize
clear eigmodes
%% Calculating the residual of continuity equation :
%  D = 2*pi*i*kx*Phi_x(r) + (im/r)*Phi_theta(r) + (1/r)*d/dr(r*Phi_r(r))

Uc = 1; % Convection velocity
D_continuity = 0;
for k = 1:1
    for j = 6:6
        Dx = 0; Dtheta = 0; D_continuity = 0;
        Dx(1:nr-1,1)     = sqrt(-1).*2*pi.*(-f(j)/Uc).*0.5.*(w_eigenmode(2:nr,k,j) + w_eigenmode(1:nr-1,k,j));
        Dtheta(1:nr-1,1) = 2.*(sqrt(-1)*mode).*0.5.*(v_eigenmode(2:nr,k,j) + v_eigenmode(1:nr-1,k,j))./(rc(2:nr,1) + rc(1:nr-1,1));
        for a = 2:nr
            Dr(a-1,1) = (2/(rc(a,1) + rc(a-1,1)))*(rc(a,1)*u_eigenmode(a,k,j) - rc(a-1,1)*u_eigenmode(a-1,k,j)) ...
                        /(rc(a,1) - rc(a-1,1));
        end
        D_continuity(1:nr-1,1,1) = Dx + Dtheta + Dr;
    end
end


D_continuity_abs = abs(D_continuity)/(0.014/2.6);
plot(rc(1:nr-1), D_continuity_abs,'k--');

%% Taking the maximum of D_continuity

D_continuity = abs(D_continuity);

for k = 1:10
    for j = 1:size(D_continuity,3)
        D_continuity_max(k,j) = max(D_continuity(:,k,j));
    end
end

[n kx] = meshgrid(f(2:81), linspace(1,10,10));

h1 = surf(D_continuity_max', n, kx);
