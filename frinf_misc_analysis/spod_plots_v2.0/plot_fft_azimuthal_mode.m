%% Name - Sheel Nidhan
%  Date - 29th February, 2020
%  Plots for the resolvent gains obtained from the 1D resolvent analysis of wake

clear; close all;
dirout = './';
%% Generating the figure object

close all;
x0=0;
y0=0;
width=15;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);
[ha, pos] = tight_subplot(2,2,[0.1, 0.05],[.1,.1],[.1 .05]);

%% Read the grid file

nr  = 354;
ntheta = 256;

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*SAGROW> % Centered the grid faces to grid centers
end

%% m = 1, St = 0.135 plot

load('x_D_10_fft.mat', 'vel_m1', 'vel_m2', 'nr_sampled', 'nt', 'r');

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

dt = 0.07;
t = dt*linspace(0,nt-1,nt)';
[R T] = meshgrid(rc(1:nr_sampled), t);

axes(ha(1));

[C,h] = contourf(T,R,real(vel_m1)', 20);
colormap hot;
set(h, 'edgecolor', 'none')
colorbar;
caxis([-0.05, 0.05]);
hold on;
hline = plot(t,1*ones(size(t,1),1),'k--', 'Linewidth',3);

ylim([0 5]);

ax = gca;
ax.FontSize = 20;

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(1).Position(1) ha(1).Position(2)+0.375 ha(1).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.05:0.025:0.05);
set(hBar, 'TickLabelInterpreter', 'latex');

box on;

dim = [ha(1).Position(1)+0.35*ha(1).Position(3) .6 .3 .3];
str = '$\Re(\tilde{u}_{x}(x/D = 10, r/D, m=1,t))$';
a = annotation('textbox',dim,'interpreter', 'latex', 'String',str,'FitBoxToText','on', 'EdgeColor', 'none');
a.FontSize = 20;

hXLabel = xlabel('$tU_{\infty}/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$r/D$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.15, 0.95, 0]);

%% FFT of time series

dt = 0.0721;
data = (vel_m1)';
Fs   = 1/dt;
L = size(t,1);
W = hamming(L);

% Plotting data at r/D \approx 0.8;
% h_fft_time_series = plot(t,   data(:,60), 'k-', 'Linewidth',2);

FFT_time_series = fft(W.*data(:,75));

P2 = abs(FFT_time_series/L).^2;
P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

axes(ha(3));

loglog(f,P1, 'k-', 'Linewidth',2); 
hold on;
%loglog(f(500:end), 10^-3*f(500:end).^(-5/3), 'k--', 'Linewidth',2);

ax = gca;
ax.FontSize = 20;

xlim([10^-3 10^1])
xticks([10^-3 10^-2 10^-1 10^0 10^1]);
xticklabels({'$10^{-3}$','$10^{-2}$','$10^{-1}$', '$10^{0}$', '$10^{1}$'});
ylim([10^-15 10^-4]);
yticks([10^-15 10^-10 10^-5]);
yticklabels({'$10^{-15}$','$10^{-10}$','$10^{-5}$'});

box on; grid on;

dim = [.3 .41 .03 .03];
a = annotation('ellipse',dim);
a.Color = 'red';
a.LineWidth = 2;

x = [0.40 0.34];
y = [0.42 0.42];
b = annotation('textarrow',x,y,'interpreter', 'latex', 'String','$\mbox{\textit{St}} \approx 0.135$');
b.LineWidth = 2;
b.FontSize = 20;
b.Color = 'red';

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$E_{\tilde{u}_{x}\tilde{u}_{x}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(c)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.15, 0.95, 0]);

%% m = 2, St = 0 plot

load('x_D_80_fft.mat', 'vel_m1', 'vel_m2', 'nr_sampled', 'nt', 'r');

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

dt = 0.07;
t = dt*linspace(0,nt-1,nt)';
[R T] = meshgrid(rc(1:nr_sampled), t);

axes(ha(2));

[C,h] = contourf(T,R,real(vel_m2)', 20);
colormap hot;
set(h, 'edgecolor', 'none')
colorbar;
caxis([-0.01, 0.01]);
hold on;
hline = plot(t,2*ones(size(t,1),1),'k--', 'Linewidth',3);
ylim([0 5]);
yticks([0 1 2 3 4 5]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.

ax = gca;
ax.FontSize = 20;

hBar = colorbar;
set(hBar, 'Location', 'northoutside');
set(hBar,'Position',[ha(2).Position(1) ha(2).Position(2)+0.375 ha(2).Position(3) 0.02])% To change size
set(hBar, 'YTick', -0.01:0.005:0.01);
set(hBar, 'TickLabelInterpreter', 'latex');

box on;

dim = [ha(2).Position(2)+0.35*ha(2).Position(3) .6 .3 .3];
str = '$\Re(\tilde{u}_{x}(x/D = 80, r/D, m=2,t))$';
a = annotation('textbox',dim,'interpreter', 'latex', 'String',str,'FitBoxToText','on', 'EdgeColor', 'none');
a.FontSize = 20;

hXLabel = xlabel('$tU_{\infty}/D$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$r/D$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.07, 0.95, 0]);

%% FFT of time series

dt = 0.0721;
data = (vel_m2)';
Fs   = 1/dt;
L = size(t,1);
W = ones(L);

% Plotting data at r/D \approx 0.8;
% h_fft_time_series = plot(t,   data(:,60), 'k-', 'Linewidth',2);

FFT_time_series = fft(W.*data(:,137));

P2 = abs(FFT_time_series/L).^2;
P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

axes(ha(4));

loglog(f,P1, 'k-', 'Linewidth',2); 
hold on;
%loglog(f(500:end), 10^-3*f(500:end).^(-5/3), 'k--', 'Linewidth',2);

ax = gca;
ax.FontSize = 20;

xlim([10^-3 10^1])
xticks([10^-3 10^-2 10^-1 10^0 10^1]);
xticklabels({'$10^{-3}$','$10^{-2}$','$10^{-1}$', '$10^{0}$', '$10^{1}$'});
ylim([10^-15 10^-4]);
yticks([10^-15 10^-10 10^-5]);
set(gca,'Yticklabel',[]); %to just get rid of the numbers but leave the ticks.

box on; grid on;

dim = [.58 .36 .05 .05];
a = annotation('ellipse',dim);
a.Color = 'red';
a.LineWidth = 2;

x = [0.59 0.59];
y = [0.32 0.36];
b = annotation('textarrow',x,y,'interpreter', 'latex', 'String','Low $\mbox{\textit{St}}$');
b.LineWidth = 2;
b.FontSize = 20;
b.Color = 'red';

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
% hYLabel = ylabel('$E_{\tilde{u}_{x}\tilde{u}_{x}}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.07, 0.95, 0]);


%% Save the plot

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'fft_azimuthal_modes_m1_m2_hot', '.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'fft_azimuthal_modes_m1_m2_hot', '.eps'),'-depsc2','-r600');
