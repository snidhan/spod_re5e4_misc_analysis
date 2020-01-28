function [rho_new,u_r_new,u_th_new,u_x_new,p_new,T_new,tke_new] = get_LES_meanFlow_M04(x_new,r_new,nx_smooth,nr_smooth,xStart,pipe_idx_top,pipe_idx_bottom)

%function [RHO,U,V,W,P,T,REY_XR,TKE] = get_LES_meanFlow_M04(x_LST,r_LST,nx_smooth,nr_smooth,xStart,pipe_idx_top,pipe_idx_bottom)
%%
root    = './baseFlows/CylGrid_M04-stats/';
Ma      = 0.4;
kappa   = 1.4;

%%%%%%%%
% Grid %
%%%%%%%%
file    = [root 'CylGrid_656_138_128-stats.xyz'];
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

file    = [root 'CylGrid_656_138_128-stats.2000000.dat'];
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

%%
% figure
% subplot(1,3,1)
% pcolor(x,r,u_x_mean), shading interp, axis equal tight, hold on
% subplot(1,3,2)
% pcolor(x,r,p_mean), shading interp, axis equal tight, hold on
% subplot(1,3,3)
% pcolor(x,r,rho_mean), shading interp, axis equal tight, hold on

% change non-dimensionalization to jet U_jet
u_x_mean     = u_x_mean/Ma;
u_r_mean     = u_r_mean/Ma;
u_t_mean     = u_t_mean/Ma;
p_mean       = p_mean/Ma^2;
tke          = tke/Ma^2;

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
% x   = x_LST;
% r   = r_LST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COPY AND PASTE FROM 0.9 CASE -> SELF-SIMILAR EXTRAPOLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_inf     = rho_mean(end,1);
u_x_inf     = u_x_mean(end,1);
T_inf       = T_mean(end,1);
rho_mean    = rho_mean  - rho_inf;
u_x_mean    = u_x_mean  - u_x_inf;
T_mean      = T_mean    - T_inf;

x_ss_start  = 27;   % obtain self-similarity parameters from x in [x_ss_start,x_end]

u_x_prof    = u_x_mean(:,x_1D>x_ss_start);
u_r_prof    = u_r_mean(:,x_1D>x_ss_start);
rho_prof    = rho_mean(:,x_1D>x_ss_start);
T_prof      = T_mean(:,x_1D>x_ss_start);
tke_prof    = tke(:,x_1D>x_ss_start);
x_1D        = x_mean(1,:);
x_prof      = x_1D(x_1D>x_ss_start);

r_half      = zeros(length(x_prof),1);
u_x_0       = zeros(length(x_prof),1);
u_r_half    = zeros(length(x_prof),1);
rho_0       = zeros(length(x_prof),1);
T_0         = zeros(length(x_prof),1);
tke_0       = zeros(length(x_prof),1);
for i = 1:length(x_prof)
    r_half(i)   = interp1(u_x_prof(:,i)/u_x_prof(1,i),r_1D,0.5);
    u_x_0(i)    = u_x_prof(1,i);
    u_r_half(i) = u_x_prof(1,i); % u_r scales as u_x_0
    rho_0(i)    = rho_prof(1,i);
    T_0(i)      = T_prof(1,i);
    tke_0(i)    = tke_prof(1,i);
end
r_half      = r_half';
u_x_0       = u_x_0';
u_r_half    = u_r_half';
rho_0       = rho_0';
T_0         = T_0';
tke_0       = tke_0';

% regression for r_half and u_x_0
p_r_half    = polyfit(x_prof,r_half,1);                     % r_1/2 scales with x
fun         = @(b,x)b(1)./x+b(2);
options.Display                 = 'off';
options.OptimalityTolerance     = 1e-10;
options.FunctionTolerance       = 1e-10;
b_guess     = (1./x_prof(:))\rho_0(:);
b_rho       = lsqcurvefit(fun,[b_guess 0],x_prof,rho_0,[],[],options);    % centerline velocity scales with ~1/x
b_guess     = (1./x_prof(:))\u_x_0(:);
b_u         = lsqcurvefit(fun,[b_guess 0],x_prof,u_x_0,[],[],options);
b_guess     = (1./x_prof(:))\T_0(:);
b_T         = lsqcurvefit(fun,[b_guess 0],x_prof,T_0,[],[],options);
b_guess     = (1./x_prof(:))\tke_0(:);
b_tke       = lsqcurvefit(fun,[b_guess 0],x_prof,tke_0,[],[],options);

% mean self-similar profile
r_ss    = zeros(length(r_1D),length(x_prof));
u_x_ss  = zeros(length(r_1D),length(x_prof));
u_r_ss  = zeros(length(r_1D),length(x_prof));
rho_ss  = zeros(length(r_1D),length(x_prof));
T_ss    = zeros(length(r_1D),length(x_prof));
tke_ss  = zeros(length(r_1D),length(x_prof));
for i = 1:length(x_prof)
    r_ss(:,i)   = r_1D/r_half(i);
    u_x_ss(:,i) = u_x_prof(:,i)/u_x_0(i);
    u_r_ss(:,i) = u_r_prof(:,i)/u_r_half(i);
    rho_ss(:,i) = rho_prof(:,i)/rho_0(i);
    T_ss(:,i)   = T_prof(:,i)/T_0(i);
    tke_ss(:,i) = tke_prof(:,i)/tke_0(i);    
end
%%
r_ss    = mean(r_ss,2);
u_x_ss  = mean(u_x_ss,2);
u_r_ss  = mean(u_r_ss,2);
rho_ss  = mean(rho_ss,2);
T_ss    = mean(T_ss,2);
tke_ss  = mean(tke_ss,2);

% figure('name','u_x')
% subplot(1,3,1)
% plot(r_1D, u_x_prof)
% subplot(1,3,2)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), u_x_prof(:,i)), hold on; end
% subplot(1,3,3)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), u_x_prof(:,i)/u_x_0(i)), hold on; end
% plot(r_ss,u_x_ss,'r--','linewidth',2)
% 
% figure('name','rho')
% subplot(1,3,1)
% plot(r_1D, rho_prof)
% subplot(1,3,2)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), rho_prof(:,i)), hold on; end
% subplot(1,3,3)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), rho_prof(:,i)/rho_0(i)), hold on; end
% plot(r_ss,rho_ss,'r--','linewidth',2)
% 
% figure('name','T')
% subplot(1,3,1)
% plot(r_1D, T_prof)
% subplot(1,3,2)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), T_prof(:,i)), hold on; end
% subplot(1,3,3)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), T_prof(:,i)/T_0(i)), hold on; end
% plot(r_ss,T_ss,'r--','linewidth',2)
% 
% figure('name','u_r')
% subplot(1,3,1)
% plot(r_1D, u_r_prof)
% subplot(1,3,2)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), u_r_prof(:,i)), hold on; end
% subplot(1,3,3)
% for i = 1:length(x_prof), plot(r_1D/r_half(i), u_r_prof(:,i)/u_r_half(i)), hold on; end
% plot(r_ss,u_r_ss,'r--','linewidth',2)

[x,r]         = meshgrid(x_1D,r_1D);

u_x_new     = interp2(x,r,u_x_mean,x_new,r_new,'spline',NaN); % don't let spline extrapolate
u_r_new     = interp2(x,r,u_r_mean,x_new,r_new,'spline',NaN);
rho_new     = interp2(x,r,rho_mean,x_new,r_new,'spline',NaN);
T_new       = interp2(x,r,T_mean,x_new,  r_new,'spline',NaN);
tke_new     = interp2(x,r,tke,x_new,     r_new,'spline',0);

rho_inlet   = interp1(r_1D,rho_mean(:,1),r_new(:,1),'spline',NaN);
u_x_inlet   = interp1(r_1D,u_x_mean(:,1),r_new(:,1),'spline',NaN);
u_r_inlet   = interp1(r_1D,u_r_mean(:,1),r_new(:,1),'spline',NaN);
T_inlet     = interp1(r_1D,T_mean(:,1),  r_new(:,1),'spline',NaN);
tke_inlet   = interp1(r_1D,tke(:,1),     r_new(:,1),'spline',NaN);
% constant inlet solution for x<x0_LES

for xi = 1:length(x_new)
    if x_new(1,xi)<x(1)
        rho_new(:,xi)   = rho_inlet;
        u_x_new(:,xi)   = u_x_inlet;
        u_r_new(:,xi)   = u_r_inlet;
        T_new(:,xi)     = T_inlet;
        tke_new(:,xi)   = tke_inlet;
    end
end

% self-similar solution for x>x1_LES
for xi = 1:length(x_new)
    if x_new(1,xi)>x(end)
        x_loc           = x_new(1,xi);
        r_half_loc      = p_r_half(1)*x_loc + p_r_half(2);
        u0_loc          = b_u(1)/x_loc      + b_u(2);
        rho0_loc        = b_rho(1)/x_loc    + b_rho(2);
        T0_loc          = b_T(1)/x_loc      + b_T(2);
        tke0_loc        = b_tke(1)/x_loc    + b_tke(2);
        u_x_loc         = u_x_ss*u0_loc;
        u_r_loc         = u_r_ss*u0_loc;
        rho_loc         = rho_ss*rho0_loc;
        T_loc           = T_ss*T0_loc;
        tke_loc         = tke_ss*tke0_loc;
        r_loc           = r_ss*r_half_loc;
        u_x_new(:,xi)   = interp1(r_loc,u_x_loc, r_new(:,1));
        u_r_new(:,xi)   = interp1(r_loc,u_r_loc, r_new(:,1));
        rho_new(:,xi)   = interp1(r_loc,rho_loc, r_new(:,1));
        T_new(:,xi)     = interp1(r_loc,T_loc,   r_new(:,1));
        tke_new(:,xi)   = interp1(r_loc,tke_loc, r_new(:,1));
    end
end

rho_new     = rho_new + rho_inf;
u_x_new     = u_x_new + u_x_inf;
T_new       = T_new   + T_inf;

for i = 1:length(x_new)
    f_nan   = isnan(u_x_new(:,i));
    noIsnan = sum(f_nan);
    rho_new(f_nan,i) = rho_new(end-noIsnan,i);
    u_x_new(f_nan,i) = u_x_new(end-noIsnan,i);
    u_r_new(f_nan,i) = u_r_new(end-noIsnan,i);
    T_new(f_nan,i)   = T_new(end-noIsnan,i);
    tke_new(f_nan,i) = tke_new(end-noIsnan,i);
end
u_th_new    = zeros(size(x_new));
p_new       = (rho_new.*T_new)/kappa/Ma^2;

% smooth data in r
x_gt_0_i    = findnearest(0,x_new(1,:));
if nr_smooth>0
    u_x_new(pipe_idx_top) = 0; u_x_new(pipe_idx_bottom) = 0; u_r_new(pipe_idx_top) = 0; u_r_new(pipe_idx_bottom) = 0;
    rho_new     = triSmooth(rho_new,1,nr_smooth);
    u_x_new(:,x_gt_0_i:end)       = triSmooth(u_x_new(:,x_gt_0_i:end),  1,nr_smooth);
    u_r_new(:,x_gt_0_i:end)       = triSmooth(u_r_new(:,x_gt_0_i:end),  1,nr_smooth);
    T_new       = triSmooth(T_new,  1,nr_smooth);
    tke_new     = triSmooth(tke_new,  1,nr_smooth);
    p_new       = triSmooth(p_new,  1,nr_smooth);
    u_x_new(pipe_idx_top) = 0; u_x_new(pipe_idx_bottom) = 0; u_r_new(pipe_idx_top) = 0; u_r_new(pipe_idx_bottom) = 0;
end
% smooth data in x
if nx_smooth>0
    rho_new     = triSmooth(rho_new,    2,nx_smooth);
    u_x_new     = triSmooth(u_x_new,    2,nx_smooth);
    u_r_new     = triSmooth(u_r_new,    2,nx_smooth);
    T_new       = triSmooth(T_new,      2,nx_smooth);
    tke_new     = triSmooth(tke_new,    2,nx_smooth);  
end

p_new               = (rho_new.*T_new)/kappa/Ma^2;
tke_new(tke_new<0)  = 0; 
tke_new(x_new<0)    = 0; 
% figure
% subplot(5,1,1)
% pcolor(x_new,r_new,rho_new), shading interp, axis equal tight, ca = caxis;
% hold on, contour(x_new,r_new,u_x_new-u_x_new(end,1),[0.01 0.01],'r--'), caxis(ca);
% subplot(5,1,2)
% pcolor(x_new,r_new,u_x_new), shading interp, axis equal tight, ca = caxis;
% hold on, contour(x_new,r_new,u_x_new-u_x_new(end,1),[0.01 0.01],'r--'), caxis(ca);
% subplot(5,1,3)
% pcolor(x_new,r_new,rho_new), shading interp, axis equal tight, ca = caxis;
% hold on, contour(x_new,r_new,u_x_new-u_x_new(end,1),[0.01 0.01],'r--'), caxis(ca);
% subplot(5,1,4)
% pcolor(x_new,r_new,T_new), shading interp, axis equal tight, ca = caxis;
% hold on, contour(x_new,r_new,u_x_new-u_x_new(end,1),[0.01 0.01],'r--'), caxis(ca);
% subplot(5,1,5)
% pcolor(x_new,r_new,rey_xr_new), shading interp, axis equal tight, ca = caxis;
% hold on, contour(x_new,r_new,u_x_new-u_x_new(end,1),[0.01 0.01],'r--'), caxis(ca);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end: COPY AND PASTE FROM 0.9 CASE -> SELF-SIMILAR EXTRAPOLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x_mean_1D   = x(1,:);
% xStart_i    = findnearest(x_mean_1D, xStart);
% if xStart_i>1
%     x_mean      = x(:,xStart_i:end);
%     r_mean      = r(:,xStart_i:end);
%     rho_mean    = rho_mean(:,xStart_i:end);
%     u_x_mean    = u_x_mean(:,xStart_i:end);
%     u_r_mean    = u_r_mean(:,xStart_i:end);
%     p_mean      = p_mean(:,xStart_i:end);
%     T_mean      = T_mean(:,xStart_i:end);  
%     tke         = tke(:,xStart_i:end);  
%     reynolds_xr = zeros(size(x_mean));%reynolds_xr(:,xStart_i:end);
% else
%     x_mean      = x;
%     r_mean      = r;
%     rho_mean    = rho_mean;
%     u_x_mean    = u_x_mean;
%     u_r_mean    = u_r_mean;
%     p_mean      = p_mean;
%     T_mean      = T_mean;
%     tke         = tke;
%     reynolds_xr = zeros(size(x_mean));%reynolds_xr;    
% end
% x   = x_LST;
% r   = r_LST;
% 
% if ~isempty(pipe_idx_top)
%     r_1D_pipe   = r_mean(r_mean(:,1)<0.5,1);
%     Nr_1D_pipe  = length(r_1D_pipe);
%     u_x_mean(Nr_1D_pipe+1:end,1) = 0;
%     u_r_mean(:,1)                = 0;
% end
% 
% % enlagred domain for interpolation
% x_mean_extrap       = zeros(size(x_mean)+[1 2]);
% r_mean_extrap       = x_mean_extrap;
% rho_mean_extrap     = x_mean_extrap;
% u_x_mean_extrap     = x_mean_extrap;
% u_r_mean_extrap     = x_mean_extrap;
% p_mean_extrap       = x_mean_extrap;
% T_mean_extrap       = x_mean_extrap;
% rey_xr_extrap       = x_mean_extrap;
% tke_extrap          = x_mean_extrap;
% 
% x_mean_extrap(1:end-1,2:end-1)  =  x_mean;
% x_mean_extrap(end,2:end-1)      =  x_mean(end,:);
% x_mean_extrap(:,1)              =  -100; 
% x_mean_extrap(:,end)            =  100; 
% 
% r_mean_extrap(1:end-1,2:end-1)  =  r_mean;
% r_mean_extrap(end,:)            =  100;
% r_mean_extrap(1:end-1,1)        =  r_mean(:,1);
% r_mean_extrap(1:end-1,end)      =  r_mean(:,end);
% 
% rho_mean_extrap(1:end-1,2:end-1)  =  rho_mean;
% rho_mean_extrap(end,1)            =  rho_mean(end,1);
% rho_mean_extrap(end,end)          =  rho_mean(end,end);
% rho_mean_extrap(end,2:end-1)      =  rho_mean(end,:);
% rho_mean_extrap(1:end-1,1)        =  rho_mean(:,1);
% rho_mean_extrap(1:end-1,end)      =  rho_mean(:,end);
% 
% u_x_mean_extrap(1:end-1,2:end-1)  =  u_x_mean;
% u_x_mean_extrap(end,1)            =  u_x_mean(end,1);
% u_x_mean_extrap(end,end)          =  u_x_mean(end,end);
% u_x_mean_extrap(end,2:end-1)      =  u_x_mean(end,:);
% u_x_mean_extrap(1:end-1,1)        =  u_x_mean(:,1);
% u_x_mean_extrap(1:end-1,end)      =  u_x_mean(:,end);
% 
% u_r_mean_extrap(1:end-1,2:end-1)  =  u_r_mean;
% u_r_mean_extrap(end,1)            =  u_r_mean(end,1);
% u_r_mean_extrap(end,end)          =  u_r_mean(end,end);
% u_r_mean_extrap(end,2:end-1)      =  u_r_mean(end,:);
% u_r_mean_extrap(1:end-1,1)        =  u_r_mean(:,1);
% u_r_mean_extrap(1:end-1,end)      =  u_r_mean(:,end);
% 
% p_mean_extrap(1:end-1,2:end-1)  =  p_mean;
% p_mean_extrap(end,1)            =  p_mean(end,1);
% p_mean_extrap(end,end)          =  p_mean(end,end);
% p_mean_extrap(end,2:end-1)      =  p_mean(end,:);
% p_mean_extrap(1:end-1,1)        =  p_mean(:,1);
% p_mean_extrap(1:end-1,end)      =  p_mean(:,end);
% 
% T_mean_extrap(1:end-1,2:end-1)  =  T_mean;
% T_mean_extrap(end,1)            =  T_mean(end,1);
% T_mean_extrap(end,end)          =  T_mean(end,end);
% T_mean_extrap(end,2:end-1)      =  T_mean(end,:);
% T_mean_extrap(1:end-1,1)        =  T_mean(:,1);
% T_mean_extrap(1:end-1,end)      =  T_mean(:,end);
% 
% rey_xr_extrap(1:end-1,2:end-1)  =  reynolds_xr;
% rey_xr_extrap(end,1)            =  reynolds_xr(end,1);
% rey_xr_extrap(end,end)          =  reynolds_xr(end,end);
% rey_xr_extrap(end,2:end-1)      =  reynolds_xr(end,:);
% rey_xr_extrap(1:end-1,1)        =  reynolds_xr(:,1);
% rey_xr_extrap(1:end-1,end)      =  reynolds_xr(:,end);
% 
% tke_extrap(1:end-1,2:end-1)  =  tke;
% tke_extrap(end,1)            =  tke(end,1);
% tke_extrap(end,end)          =  tke(end,end);
% tke_extrap(end,2:end-1)      =  tke(end,:);
% tke_extrap(1:end-1,1)        =  tke(:,1);
% tke_extrap(1:end-1,end)      =  tke(:,end);
% 
% % figure, contourf(x_mean_extrap,r_mean_extrap,rho_mean_extrap,100,'edgecolor','none')
% % figure, contourf(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,100,'edgecolor','none')
% % figure, contourf(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,100,'edgecolor','none')
% % figure, contourf(x_mean_extrap,r_mean_extrap,p_mean_extrap,100,'edgecolor','none')
% % figure, contourf(x_mean_extrap,r_mean_extrap,T_mean_extrap,100,'edgecolor','none')
% 
% r1      = r_mean(end,1);
% x0      = x_mean(1,1);
% x1      = x_mean(1,end);
% inner_i = x>x0&x<x1&r<r1;
% [nr,nx] = size(x);
% 
% % linear extrapolation, cubic interpolation
% RHO             = interp2(x_mean_extrap,r_mean_extrap,rho_mean_extrap,x,r,'linear');
% RHO(inner_i)    = interp2(x_mean_extrap,r_mean_extrap,rho_mean_extrap,x(inner_i),r(inner_i),'spline');
% W               = interp2(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,x,r,'linear');
% W(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,u_x_mean_extrap,x(inner_i),r(inner_i),'spline');
% U               = interp2(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,x,r,'linear');
% U(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,u_r_mean_extrap,x(inner_i),r(inner_i),'spline');
% T               = interp2(x_mean_extrap,r_mean_extrap,T_mean_extrap,x,r,'linear');
% T(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,T_mean_extrap,x(inner_i),r(inner_i),'spline');
% P               = interp2(x_mean_extrap,r_mean_extrap,p_mean_extrap,x,r,'linear');
% P(inner_i)      = interp2(x_mean_extrap,r_mean_extrap,p_mean_extrap,x(inner_i),r(inner_i),'spline');
% REY_XR          = interp2(x_mean_extrap,r_mean_extrap,rey_xr_extrap,x,r,'linear');
% REY_XR(inner_i) = interp2(x_mean_extrap,r_mean_extrap,rey_xr_extrap,x(inner_i),r(inner_i),'spline');
% V               = 0*x;
% TKE             = interp2(x_mean_extrap,r_mean_extrap,tke_extrap,x,r,'linear');
% TKE(inner_i)    = interp2(x_mean_extrap,r_mean_extrap,tke_extrap,x(inner_i),r(inner_i),'spline');
% 
% x_gt_0_i        = 1;%findnearest(0,x(1,:));
% 
% % figure, contourf(x(:,x_gt_0_i:end),r(:,x_gt_0_i:end),W(:,x_gt_0_i:end),100,'edgecolor','none')
% % figure, contourf(x,r,W,100,'edgecolor','none')
% 
% % smooth data in r
% if nr_smooth>0
%     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
%     RHO     = triSmooth(RHO,1,nr_smooth);
%     W(:,x_gt_0_i:end)       = triSmooth(W(:,x_gt_0_i:end),  1,nr_smooth);
%     U(:,x_gt_0_i:end)       = triSmooth(U(:,x_gt_0_i:end),  1,nr_smooth);
% %     W       = triSmooth(W,  1,nr_smooth);
% %     U       = triSmooth(U,  1,nr_smooth);
%     T       = triSmooth(T,  1,nr_smooth);
%     P       = triSmooth(P,  1,nr_smooth);
%     TKE     = triSmooth(TKE,1,nr_smooth);
%     REY_XR  = triSmooth(REY_XR,  1,nr_smooth);
%     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
% end
% % smooth data in x
% if nx_smooth>0
% %     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
%     RHO     = triSmooth(RHO,2,nx_smooth);
%     W       = triSmooth(W,  2,nx_smooth);
%     U       = triSmooth(U,  2,nx_smooth);
%     T       = triSmooth(T,  2,nx_smooth);
%     P       = triSmooth(P,  2,nx_smooth);
%     TKE     = triSmooth(TKE,  2,nx_smooth);
%     REY_XR  = triSmooth(REY_XR,  2,nx_smooth);
% %     W(pipe_idx_top) = 0; W(pipe_idx_bottom) = 0; U(pipe_idx_top) = 0; U(pipe_idx_bottom) = 0;
% end
% 
% % figure, contourf(x,r,RHO,100,'edgecolor','none'), drawnow