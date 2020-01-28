function [rho_new,u_r_new,u_th_new,u_x_new,p_new,T_new,rey_xr_new] = get_LES_meanFlow(x_new,r_new,nx_smooth,nr_smooth,xStart,pipe_idx_top,pipe_idx_bottom)

load('./baseFlows/CylGrid_656_138_longTimeFavreMean.mat'); Ma = 0.9; kappa = 1.4;

r_1D        = r_mean(:,1);
x_1D        = x_mean(1,:);

xStart_i    = findnearest(x_1D, xStart);
if xStart_i>1
    x_1D        = x_1D(xStart_i:end);
    x_mean      = x_mean(:,xStart_i:end);
    r_mean      = r_mean(:,xStart_i:end);
    rho_mean    = rho_mean(:,xStart_i:end);
    u_x_mean    = u_x_favre(:,xStart_i:end);
    u_r_mean    = u_r_favre(:,xStart_i:end);
    p_mean      = p_mean(:,xStart_i:end);
    T_mean      = T_mean(:,xStart_i:end);    
    reynolds_xr = reynolds_xr(:,xStart_i:end);        
end

if ~isempty(pipe_idx_top)
    r_1D_pipe   = r_mean(r_mean(:,1)<0.5,1);
    Nr_1D_pipe  = length(r_1D_pipe);
    u_x_mean(Nr_1D_pipe+1:end,1) = 0;
    u_r_mean(:,1)                = 0;
end

rho_inf     = rho_mean(end,1);
u_x_inf     = u_x_mean(end,1);
T_inf       = T_mean(end,1);
rho_mean    = rho_mean  - rho_inf;
u_x_mean    = u_x_mean  - u_x_inf;
T_mean      = T_mean    - T_inf;

%%
x_ss_start  = 27;   % obtain self-similarity parameters from x in [x_ss_start,x_end]

u_x_prof    = u_x_mean(:,x_1D>x_ss_start);
u_r_prof    = u_r_mean(:,x_1D>x_ss_start);
rho_prof    = rho_mean(:,x_1D>x_ss_start);
T_prof      = T_mean(:,x_1D>x_ss_start);
x_1D        = x_mean(1,:);
x_prof      = x_1D(x_1D>x_ss_start);

r_half      = zeros(length(x_prof),1);
u_x_0       = zeros(length(x_prof),1);
u_r_half    = zeros(length(x_prof),1);
rho_0       = zeros(length(x_prof),1);
T_0         = zeros(length(x_prof),1);
for i = 1:length(x_prof)
    r_half(i)   = interp1(u_x_prof(:,i)/u_x_prof(1,i),r_1D,0.5);
    u_x_0(i)    = u_x_prof(1,i);
    u_r_half(i) = u_x_prof(1,i); % u_r scales as u_x_0
    rho_0(i)    = rho_prof(1,i);
    T_0(i)      = T_prof(1,i);
end
r_half      = r_half';
u_x_0       = u_x_0';
u_r_half    = u_r_half';
rho_0       = rho_0';
T_0         = T_0';

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

% mean self-similar profile
r_ss    = zeros(length(r_1D),length(x_prof));
u_x_ss  = zeros(length(r_1D),length(x_prof));
u_r_ss  = zeros(length(r_1D),length(x_prof));
rho_ss  = zeros(length(r_1D),length(x_prof));
T_ss    = zeros(length(r_1D),length(x_prof));
for i = 1:length(x_prof)
    r_ss(:,i)   = r_1D/r_half(i);
    u_x_ss(:,i) = u_x_prof(:,i)/u_x_0(i);
    u_r_ss(:,i) = u_r_prof(:,i)/u_r_half(i);
    rho_ss(:,i) = rho_prof(:,i)/rho_0(i);
    T_ss(:,i)   = T_prof(:,i)/T_0(i);
end
%%
r_ss    = mean(r_ss,2);
u_x_ss  = mean(u_x_ss,2);
u_r_ss  = mean(u_r_ss,2);
rho_ss  = mean(rho_ss,2);
T_ss    = mean(T_ss,2);

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
T_new       = interp2(x,r,T_mean,x_new,r_new,'spline',NaN);
rey_xr_new  = interp2(x,r,reynolds_xr,x_new,r_new,'spline',0);

rho_inlet   = interp1(r_1D,rho_mean(:,1),r_new(:,1),'spline',NaN);
u_x_inlet   = interp1(r_1D,u_x_mean(:,1),r_new(:,1),'spline',NaN);
u_r_inlet   = interp1(r_1D,u_r_mean(:,1),r_new(:,1),'spline',NaN);
T_inlet     = interp1(r_1D,T_mean(:,1),  r_new(:,1),'spline',NaN);
% rey_xr_inlet= interp1(r_1D,reynolds_xr(:,1),  r_new(:,1),'spline',NaN);
% constant inlet solution for x<x0_LES

for xi = 1:length(x_new)
    if x_new(1,xi)<x(1)
        rho_new(:,xi)   = rho_inlet;
        u_x_new(:,xi)   = u_x_inlet;
        u_r_new(:,xi)   = u_r_inlet;
        T_new(:,xi)     = T_inlet;
%         rey_xr_new(:,xi)= rey_xr_inlet;
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
        u_x_loc         = u_x_ss*u0_loc;
        u_r_loc         = u_r_ss*u0_loc;
        rho_loc         = rho_ss*rho0_loc;
        T_loc           = T_ss*T0_loc;
        r_loc           = r_ss*r_half_loc;
        u_x_new(:,xi)   = interp1(r_loc,u_x_loc, r_new(:,1));
        u_r_new(:,xi)   = interp1(r_loc,u_r_loc, r_new(:,1));
        rho_new(:,xi)   = interp1(r_loc,rho_loc,r_new(:,1));
        T_new(:,xi)     = interp1(r_loc,T_loc,  r_new(:,1));
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
    rey_xr_new  = triSmooth(rey_xr_new,  1,nr_smooth);
    p_new       = triSmooth(p_new,  1,nr_smooth);
    u_x_new(pipe_idx_top) = 0; u_x_new(pipe_idx_bottom) = 0; u_r_new(pipe_idx_top) = 0; u_r_new(pipe_idx_bottom) = 0;
end
% smooth data in x
if nx_smooth>0
    rho_new     = triSmooth(rho_new,    2,nx_smooth);
    u_x_new     = triSmooth(u_x_new,    2,nx_smooth);
    u_r_new     = triSmooth(u_r_new,    2,nx_smooth);
    T_new       = triSmooth(T_new,      2,nx_smooth);
    rey_xr_new  = triSmooth(rey_xr_new, 2,nx_smooth);
    p_new       = triSmooth(p_new,      2,nx_smooth);    
end

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

p_new   = (rho_new.*T_new)/kappa/Ma^2;
% 
% u_x_new(x_new<0) = 0; 
% u_r_new(x_new<0) = 0;