function [RHO,U,V,W,P,T] = get_LES_meanFlow_69M(x,r,nx_smooth,nr_smooth,xStart)

load('aux/CylGrid_656_138_longTimeMean_69M.mat');
disp(['    NOTE: loading the 69M base-state']);

x_mean_1D   = x_mean(1,:);
xStart_i    = findnearest(x_mean_1D, xStart);
if xStart_i>1
    x_mean      = x_mean(:,xStart_i:end);
    r_mean      = r_mean(:,xStart_i:end);
    rho_mean    = rho_mean(:,xStart_i:end);
    u_x_mean    = u_x_mean(:,xStart_i:end);
    u_r_mean    = u_r_mean(:,xStart_i:end);
    p_mean      = p_mean(:,xStart_i:end);
    T_mean      = T_mean(:,xStart_i:end);    
end

% enlagred domain for interpolation
x_mean_extrap       = zeros(size(x_mean)+[1 2]);
r_mean_extrap       = x_mean_extrap;
rho_mean_extrap     = x_mean_extrap;
u_x_mean_extrap     = x_mean_extrap;
u_r_mean_extrap     = x_mean_extrap;
p_mean_extrap       = x_mean_extrap;
T_mean_extrap       = x_mean_extrap;

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
V               = 0*x;

% smooth data in r
if nr_smooth>0
    RHO = triSmooth(RHO,1,nr_smooth);
    W   = triSmooth(W,  1,nr_smooth);
    U   = triSmooth(U,  1,nr_smooth);
    T   = triSmooth(T,  1,nr_smooth);
    P   = triSmooth(P,  1,nr_smooth);
end
% smooth data in x
if nx_smooth>0
    RHO = triSmooth(RHO,2,nx_smooth);
    W   = triSmooth(W,  2,nx_smooth);
    U   = triSmooth(U,  2,nx_smooth);
    T   = triSmooth(T,  2,nx_smooth);
    P   = triSmooth(P,  2,nx_smooth);
end

% figure, contourf(x,r,RHO,100,'edgecolor','none')