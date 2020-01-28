function [RHO,U,V,W,P,T,REY_XR,TKE] = get_LES_meanFlow_M04(x_LST,r_LST,nx_smooth,nr_smooth,xStart,pipe_idx_top,pipe_idx_bottom)
%%
root    = '/mnt/Jetyard/M09/LES/10M/CylGrid-stats2000000/JetPlume/';
mach    = 0.9;

%%%%%%%%
% Grid %
%%%%%%%%
file    = [root '10M_CylGrid_656_138_128-stats.xyz'];
nx      = 656;
ny      = 138;
nz      = 128;

nr      = ny;
nth     = nz;

fileID  = fopen(file);
np      = fread(fileID,1,'int');
format  = fread(fileID,1,'int');
xyz     = fread(fileID,[3,np],'float');
fclose(fileID);

x       = reshape(squeeze(xyz(1,:)),nz,ny,nx);
y       = reshape(squeeze(xyz(2,:)),nz,ny,nx);
z       = reshape(squeeze(xyz(3,:)),nz,ny,nx);

x_1D    = squeeze(x(1,1,:));

% azimuthal angle
th      = 2*pi*[0:nth-1]'/nth;
r_1D    = squeeze(y(1,:,1));
[x,r]   = meshgrid(x_1D,r_1D);

%%%%%%%%
% Mean %
%%%%%%%%

file    = [root '10M_CylGrid_656_138_128-stats.2000000.dat'];
fileID  = fopen(file);

np      = fread(fileID,1,'int');
nvar    = fread(fileID,1,'int');
step    = fread(fileID,1,'int');
format  = fread(fileID,1,'int');
time    = fread(fileID,1,'float');
data    = fread(fileID,[nvar,np],'float');
fclose(fileID);

%  1 "cv2no(RHO_AVG)"
%  2 "cv2no(T_AVG)"
%  3 "cv2no(P_AVG-1/1.4)"
%  4 "cv2no(U_AVG)-x"
%  5 "cv2no(U_AVG)-y"
%  6 "cv2no(U_AVG)-z"
%  7 "cv2no(RHO_RMS)"
%  8 "cv2no(T_RMS)"
%  9 "cv2no(P_RMS)"
%  10 "cv2no(U_RMS)-x"
%  11 "cv2no(U_RMS)-y"
%  12 "cv2no(U_RMS)-z"
%  13 "cv2no(U_REY)-x"
%  14 "cv2no(U_REY)-y"
%  15 "cv2no(U_REY)-z"

rho_mean    = reshape(squeeze(data(1,:)),nth,nr,nx);
T_mean      = reshape(squeeze(data(2,:)),nth,nr,nx);
p_mean      = reshape(squeeze(data(3,:)),nth,nr,nx)+1/1.4;
u_x_mean    = reshape(squeeze(data(4,:)),nth,nr,nx);
u_y_mean    = reshape(squeeze(data(5,:)),nth,nr,nx);
u_z_mean    = reshape(squeeze(data(6,:)),nth,nr,nx);
rho_rms     = reshape(squeeze(data(7,:)),nth,nr,nx);
T_rms       = reshape(squeeze(data(8,:)),nth,nr,nx);
p_rms       = reshape(squeeze(data(9,:)),nth,nr,nx);
u_x_rms     = reshape(squeeze(data(10,:)),nth,nr,nx);
u_y_rms     = reshape(squeeze(data(11,:)),nth,nr,nx);
u_z_rms     = reshape(squeeze(data(12,:)),nth,nr,nx);
u_yz_rey    = reshape(squeeze(data(13,:)),nth,nr,nx);
u_xz_rey    = reshape(squeeze(data(14,:)),nth,nr,nx);
u_xy_rey    = reshape(squeeze(data(15,:)),nth,nr,nx);

% cylindrical velocities
th           = 2*pi*[0:nth-1]'/nth;
u_r_mean     = zeros(nth,nr,nx);
u_t_mean     = zeros(nth,nr,nx);
u_r_rms      = zeros(nth,nr,nx);
u_t_rms      = zeros(nth,nr,nx);
for jr = 1:nr
    for jx = 1:nx
        u_r_mean(:,jr,jx) =  cos(th).*u_y_mean(:,jr,jx) + sin(th).*u_z_mean(:,jr,jx);
        u_t_mean(:,jr,jx) = -sin(th).*u_y_mean(:,jr,jx) + cos(th).*u_z_mean(:,jr,jx);
        u_r_rms(:,jr,jx)  =  cos(th).*u_y_rms(:,jr,jx)  + sin(th).*u_z_rms(:,jr,jx);
        u_t_rms(:,jr,jx)  = -sin(th).*u_y_rms(:,jr,jx)  + cos(th).*u_z_rms(:,jr,jx);    
    end
end

% average 3-D mean in azimuth
tke         = squeeze(mean(0.5*rho_mean.*(u_x_rms.^2 + u_y_rms.^2 + u_z_rms.^2),1));
rho_mean    = squeeze(mean(rho_mean,1));
u_x_mean    = squeeze(mean(u_x_mean,1));
u_r_mean    = squeeze(mean(u_r_mean,1));
u_t_mean    = squeeze(mean(u_t_mean,1));
p_mean      = squeeze(mean(p_mean,1));
T_mean      = squeeze(mean(T_mean,1));

% %
% figure
% subplot(1,3,1)
% pcolor(x,r,u_x_mean), shading interp, axis equal tight, hold on
% subplot(1,3,2)
% pcolor(x,r,p_mean), shading interp, axis equal tight, hold on
% subplot(1,3,3)
% pcolor(x,r,rho_mean), shading interp, axis equal tight, hold on

% change non-dimensionalization to jet U_jet
u_x_mean     = u_x_mean/mach;
u_r_mean     = u_r_mean/mach;
u_t_mean     = u_t_mean/mach;
p_mean       = p_mean/mach^2;
tke          = tke/mach^2;

x_mean_1D   = x(1,:);
xStart_i    = findnearest(x_mean_1D, xStart);
if xStart_i>1
    x_mean      = x(:,xStart_i:end);
    r_mean      = r(:,xStart_i:end);
    rho_mean    = rho_mean(:,xStart_i:end);
    u_x_mean    = u_x_mean(:,xStart_i:end);
    u_r_mean    = u_r_mean(:,xStart_i:end);
    p_mean      = p_mean(:,xStart_i:end);
    T_mean      = T_mean(:,xStart_i:end);  
    tke         = tke(:,xStart_i:end);  
    reynolds_xr = zeros(size(x_mean));%reynolds_xr(:,xStart_i:end);
else
    x_mean      = x;
    r_mean      = r;
    rho_mean    = rho_mean;
    u_x_mean    = u_x_mean;
    u_r_mean    = u_r_mean;
    p_mean      = p_mean;
    T_mean      = T_mean;
    tke         = tke;
    reynolds_xr = zeros(size(x_mean));%reynolds_xr;    
end
x   = x_LST;
r   = r_LST;

if ~isempty(pipe_idx_top)
    r_1D_pipe   = r_mean(r_mean(:,1)<0.5,1);
    Nr_1D_pipe  = length(r_1D_pipe);
    u_x_mean(Nr_1D_pipe+1:end,1) = 0;
    u_r_mean(:,1)                = 0;
end

% enlagred domain for interpolation
x_mean_extrap       = zeros(size(x_mean)+[1 2]);
r_mean_extrap       = x_mean_extrap;
rho_mean_extrap     = x_mean_extrap;
u_x_mean_extrap     = x_mean_extrap;
u_r_mean_extrap     = x_mean_extrap;
p_mean_extrap       = x_mean_extrap;
T_mean_extrap       = x_mean_extrap;
rey_xr_extrap       = x_mean_extrap;
tke_extrap          = x_mean_extrap;

x_mean_extrap(1:end-1,2:end-1)  =  x_mean;
x_mean_extrap(end,2:end-1)      =  x_mean(end,:);
x_mean_extrap(:,1)              =  -100; 
x_mean_extrap(:,end)            =  100; 

r_mean_extrap(1:end-1,2:end-1)  =  r_mean;
r_mean_extrap(end,:)            =  100;
r_mean_extrap(1:end-1,1)        =  r_mean(:,1);
r_mean_extrap(1:end-1,end)      =  r_mean(:,end);

rho_mean_extrap(1:end-1,2:end-1)  =  rho_mean;
rho_mean_extrap(end,1)            =  rho_mean(end,1);
rho_mean_extrap(end,end)          =  rho_mean(end,end);
rho_mean_extrap(end,2:end-1)      =  rho_mean(end,:);
rho_mean_extrap(1:end-1,1)        =  rho_mean(:,1);
rho_mean_extrap(1:end-1,end)      =  rho_mean(:,end);

u_x_mean_extrap(1:end-1,2:end-1)  =  u_x_mean;
u_x_mean_extrap(end,1)            =  u_x_mean(end,1);
u_x_mean_extrap(end,end)          =  u_x_mean(end,end);
u_x_mean_extrap(end,2:end-1)      =  u_x_mean(end,:);
u_x_mean_extrap(1:end-1,1)        =  u_x_mean(:,1);
u_x_mean_extrap(1:end-1,end)      =  u_x_mean(:,end);

u_r_mean_extrap(1:end-1,2:end-1)  =  u_r_mean;
u_r_mean_extrap(end,1)            =  u_r_mean(end,1);
u_r_mean_extrap(end,end)          =  u_r_mean(end,end);
u_r_mean_extrap(end,2:end-1)      =  u_r_mean(end,:);
u_r_mean_extrap(1:end-1,1)        =  u_r_mean(:,1);
u_r_mean_extrap(1:end-1,end)      =  u_r_mean(:,end);

p_mean_extrap(1:end-1,2:end-1)  =  p_mean;
p_mean_extrap(end,1)            =  p_mean(end,1);
p_mean_extrap(end,end)          =  p_mean(end,end);
p_mean_extrap(end,2:end-1)      =  p_mean(end,:);
p_mean_extrap(1:end-1,1)        =  p_mean(:,1);
p_mean_extrap(1:end-1,end)      =  p_mean(:,end);

T_mean_extrap(1:end-1,2:end-1)  =  T_mean;
T_mean_extrap(end,1)            =  T_mean(end,1);
T_mean_extrap(end,end)          =  T_mean(end,end);
T_mean_extrap(end,2:end-1)      =  T_mean(end,:);
T_mean_extrap(1:end-1,1)        =  T_mean(:,1);
T_mean_extrap(1:end-1,end)      =  T_mean(:,end);

rey_xr_extrap(1:end-1,2:end-1)  =  reynolds_xr;
rey_xr_extrap(end,1)            =  reynolds_xr(end,1);
rey_xr_extrap(end,end)          =  reynolds_xr(end,end);
rey_xr_extrap(end,2:end-1)      =  reynolds_xr(end,:);
rey_xr_extrap(1:end-1,1)        =  reynolds_xr(:,1);
rey_xr_extrap(1:end-1,end)      =  reynolds_xr(:,end);

tke_extrap(1:end-1,2:end-1)  =  tke;
tke_extrap(end,1)            =  tke(end,1);
tke_extrap(end,end)          =  tke(end,end);
tke_extrap(end,2:end-1)      =  tke(end,:);
tke_extrap(1:end-1,1)        =  tke(:,1);
tke_extrap(1:end-1,end)      =  tke(:,end);

% figure, contourf(x_mean_extrap,r_mean_extrap,rho_mean_extrap,100,'edgecolor','none')
% figure, contourf(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,100,'edgecolor','none')
% figure, contourf(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,100,'edgecolor','none')
% figure, contourf(x_mean_extrap,r_mean_extrap,p_mean_extrap,100,'edgecolor','none')
% figure, contourf(x_mean_extrap,r_mean_extrap,T_mean_extrap,100,'edgecolor','none')

r1      = r_mean(end,1);
x0      = x_mean(1,1);
x1      = x_mean(1,end);
inner_i = x>x0&x<x1&r<r1;
[nr,nx] = size(x);

% linear extrapolation, cubic interpolation
RHO             = interp2(x_mean_extrap,r_mean_extrap,rho_mean_extrap,x,r,'linear');
RHO(inner_i)    = interp2(x_mean_extrap,r_mean_extrap,rho_mean_extrap,x(inner_i),r(inner_i),'spline');
W               = interp2(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,x,r,'linear');
W(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,x(inner_i),r(inner_i),'spline');
U               = interp2(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,x,r,'linear');
U(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,x(inner_i),r(inner_i),'spline');
T               = interp2(x_mean_extrap,r_mean_extrap,T_mean_extrap,x,r,'linear');
T(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,T_mean_extrap,x(inner_i),r(inner_i),'spline');
P               = interp2(x_mean_extrap,r_mean_extrap,p_mean_extrap,x,r,'linear');
P(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,p_mean_extrap,x(inner_i),r(inner_i),'spline');
REY_XR          = interp2(x_mean_extrap,r_mean_extrap,rey_xr_extrap,x,r,'linear');
REY_XR(inner_i) = interp2(x_mean_extrap,r_mean_extrap,rey_xr_extrap,x(inner_i),r(inner_i),'spline');
V               = 0*x;
TKE             = interp2(x_mean_extrap,r_mean_extrap,tke_extrap,x,r,'linear');
TKE(inner_i)    = interp2(x_mean_extrap,r_mean_extrap,tke_extrap,x(inner_i),r(inner_i),'spline');

x_gt_0_i        = 1;%findnearest(0,x(1,:));

% figure, contourf(x(:,x_gt_0_i:end),r(:,x_gt_0_i:end),W(:,x_gt_0_i:end),100,'edgecolor','none')
% figure, contourf(x,r,W,100,'edgecolor','none')

% smooth data in r
if nr_smooth>0
    W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
    RHO     = triSmooth(RHO,1,nr_smooth);
    W(:,x_gt_0_i:end)       = triSmooth(W(:,x_gt_0_i:end),  1,nr_smooth);
    U(:,x_gt_0_i:end)       = triSmooth(U(:,x_gt_0_i:end),  1,nr_smooth);
%     W       = triSmooth(W,  1,nr_smooth);
%     U       = triSmooth(U,  1,nr_smooth);
    T       = triSmooth(T,  1,nr_smooth);
    P       = triSmooth(P,  1,nr_smooth);
    TKE     = triSmooth(TKE,1,nr_smooth);
    REY_XR  = triSmooth(REY_XR,  1,nr_smooth);
    W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
end
% smooth data in x
if nx_smooth>0
%     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
    RHO     = triSmooth(RHO,2,nx_smooth);
    W       = triSmooth(W,  2,nx_smooth);
    U       = triSmooth(U,  2,nx_smooth);
    T       = triSmooth(T,  2,nx_smooth);
    P       = triSmooth(P,  2,nx_smooth);
    TKE     = triSmooth(TKE,  2,nx_smooth);
    REY_XR  = triSmooth(REY_XR,  2,nx_smooth);
%     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
end

% figure, contourf(x,r,RHO,100,'edgecolor','none'), drawnow