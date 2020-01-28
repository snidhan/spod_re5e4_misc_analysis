function [RHO,U,V,W,P,T,REY_XR] = get_LES_meanFlow_B118(x_LST,r_LST,nx_smooth,nr_smooth,xStart,pipe_idx_top,pipe_idx_bottom)

% get grid and long-time mean
file_xyz    = '/home/oschmidt/Codes/matlab/linStab2D_cyl/baseFlows/CylGrid_B118-stats/CylGrid_698_136_128-stats.xyz';
file_stats  = '/home/oschmidt/Codes/matlab/linStab2D_cyl/baseFlows/CylGrid_B118-stats/CylGrid_698_136_128-stats.2000000.dat';
addpath('/home/oschmidt/Codes/matlab/LESdataBaseProcessing');
[x,r,th,x_1D,r_1D,nx,nr,nth]    = read_CylGrid_698_136_128_xyz(file_xyz);
% [rho_mean,u_x_mean,u_r_mean,u_t_mean,p_mean,T_mean] ...
%                                 = read_CylGrid_698_136_128_stats(root,nx,nr,nth);
[rho_mean,u_x_mean,u_r_mean,u_t_mean,p_mean,T_mean,R_rr,R_rt,R_rz,R_tt,R_tz,R_zz] = read_CylGrid_698_136_128_stats(file_stats,nx,nr,nth);

% change non-dimensionalization to jet U_jet
mach         = 1.5;
u_x_mean     = u_x_mean/mach;
u_r_mean     = u_r_mean/mach;
u_t_mean     = u_t_mean/mach;
p_mean       = p_mean/mach^2;
    
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
    reynolds_xr = zeros(size(x_mean));%reynolds_xr(:,xStart_i:end);
else
    x_mean      = x;
    r_mean      = r;
    rho_mean    = rho_mean;
    u_x_mean    = u_x_mean;
    u_r_mean    = u_r_mean;
    p_mean      = p_mean;
    T_mean      = T_mean;    
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
    REY_XR  = triSmooth(REY_XR,  2,nx_smooth);
%     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
end

% figure, contourf(x,r,RHO,100,'edgecolor','none'), drawnow