function [] = reynolds_stress_plot(x)


%% Name - Sheel Nidhan
%  Date - 13th February 2020

%% Read file

dir_in = '/home/sheel/Work2/projects_data/spod_re5e4/fr2/reystresses/';
filename = strcat(dir_in, 'reystress_x_D_', int2str(x), '.mat');
load(filename); %#ok<*LOAD>


%% Loading the grid file in radial direction

ntheta = 256;
fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');  %% Reading the radial grid
D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));
r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  %#ok<*AGROW,*SAGROW> % Centered the grid faces to grid centers
end

theta = linspace(0,2*pi,ntheta)';


%% Plotting the RS and fluxes

close all;
%% Plot reystress_ww_av - <ux'ux'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_ww_av, 10); %#ok<*ASGLU>
xlim([-4 4]);
ylim([-2 2]);
% axis equal;
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);

%% Plot reystress_uu_av - <ur'ur'>
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uu_av, 10); %#ok<*ASGLU>
xlim([-4 4]);
ylim([-2 2]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);

%% Plot reystress_vv_av - <utheta'utheta'>

figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_vv_av, 10); %#ok<*ASGLU>
xlim([-4 4]);
ylim([-2 2]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);

%% Plot reystress_densdens_av - <\rho'\rho'>

figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densdens_av, 10); %#ok<*ASGLU>
xlim([-4 4]);
ylim([-2 2]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
