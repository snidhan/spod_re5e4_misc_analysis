function [] = reynolds_stress_plot(x_D)


%% Name - Sheel Nidhan
%  Date - 13th February 2020

%% Setting up the environment

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Read file

dir_in = '/home/sheel/Work2/projects_data/spod_re5e4/fr2/reystresses/';
filename = strcat(dir_in, 'reystress_x_D_', int2str(x_D), '.mat');
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


%% Plotting the RS and fluxes in cylindrical coordinatess

close all;
%% Plot reystress_ww_av - <ux'ux'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_ww_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-2.5 2.5]);
% axis equal;
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);


%% Plot reystress_uu_av - <ur'ur'>
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uu_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-2.5 2.5]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{r}''u_{r}''\rangle$','interpreter','latex','fontsize',20);

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
hTitle = title('$\langle u_{\theta}''u_{\theta}''\rangle$','interpreter','latex','fontsize',20);

%% Plot reystress_densdens_av - <\rho'\rho'>

figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densdens_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);

%% Plot reystress_uv_av - <ur'utheta'>

figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uv_av, 10); %#ok<*ASGLU>
xlim([-4 4]);
ylim([-2 2]);
% axis equal;
set(h,'LevelList',-1:0.01:1);
set(h,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{r}''u_{\theta}''\rangle$','interpreter','latex','fontsize',20);

%% 

close all;

%% Calculating all the reynolds stresses in the cartesian coordinates

ntheta = 256;
theta_center = 0.5*2*pi/(ntheta) + linspace(0,2*pi,ntheta)';

%% NOTE 
%  z - vertical, y - spanwise, x - streamwise when calculating cartesian tensors
%% Calculatinge <uy'uz'> 
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uyuz_av(i,j) = reystress_uu_av(i,j)*sin(theta_center(j,1))*cos(theta_center(j,1)) - ...
                                 reystress_vv_av(i,j)*sin(theta_center(j,1))*cos(theta_center(j,1)) + ...
                                 reystress_uv_av(i,j)*cos(2*theta_center(j,1));
    end 
end

%% Calculating <ux'uy'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uxuy_av(i,j) = reystress_uw_av(i,j)*cos(theta_center(j,1)) - ...
                                 reystress_wv_av(i,j)*sin(theta_center(j,1));
    end 
end

%% Calculating <ux'uz'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uxuz_av(i,j) = reystress_uw_av(i,j)*sin(theta_center(j,1)) + ...
                                 reystress_wv_av(i,j)*cos(theta_center(j,1));
    end 
end

%% Calculating <uy'uy'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uyuy_av(i,j) = reystress_uu_av(i,j)*cos(theta_center(j,1))*cos(theta_center(j,1)) + ...
                                 reystress_vv_av(i,j)*sin(theta_center(j,1))*sin(theta_center(j,1)) - ...
                                 reystress_uv_av(i,j)*sin(2*(theta_center(j,1)));
    end 
end

%% Calculating <uz'uz'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uzuz_av(i,j) = reystress_uu_av(i,j)*sin(theta_center(j,1))*sin(theta_center(j,1)) + ...
                                 reystress_vv_av(i,j)*cos(theta_center(j,1))*cos(theta_center(j,1)) + ...
                                 reystress_uv_av(i,j)*sin(2*(theta_center(j,1)));
    end 
end

%% Calculating <rho'uy'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_densuy_av(i,j) = reystress_udens_av(i,j)*cos(theta_center(j,1)) - ...
                                   reystress_vdens_av(i,j)*sin(theta_center(j,1));
    end 
end

%% Calculating <rho'uz'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_densuz_av(i,j) = reystress_udens_av(i,j)*sin(theta_center(j,1)) + ...
                                   reystress_vdens_av(i,j)*cos(theta_center(j,1));
    end 
end

%% Calculating <rho'ux'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_densux_av(i,j) = reystress_wdens_av(i,j); 
    end 
end

%% Calculating <ux'ux'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_uxux_av(i,j) = reystress_ww_av(i,j); 
    end 
end

%% Calculating <rho'rho'>
%  Checked

for  i = 1:size(rc,1)
    for j = 1:size(theta_center,1)
        reystress_densdens_av(i,j) = reystress_densdens_av(i,j); 
    end 
end

%% Plotting RS in the x-y-z coordinates

%% 1. Plot reystress_uyuy_av - <uy'uy'> 

figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uyuy_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{y}''u_{y}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uyuy','x_D_',int2str(x_D), '.png'),'-dpng','-r600');

%% 2. Plot reystress_uzuz_av - <uz'uz'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uzuz_av, 10); %#ok<*ASGLU>
xlim([-7.5 7.5]);
ylim([-7.5 7.5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{z}''u_{z}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uzuz','x_D_',int2str(x_D), '.png'),'-dpng','-r600');
%% 3. Plot reystress_uxux_av - <ux'ux'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uxux_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{x}''u_{x}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uxux','x_D_',int2str(x_D), '.png'),'-dpng','-r600');
%% 4. Plot reystress_densdens_av - <rho'rho'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densdens_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle \rho''\rho''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('rhorho','x_D_',int2str(x_D), '.png'),'-dpng','-r600');
%% 5. Plot reystress_uxuz_av - <ux'uz'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uxuz_av, 10); %#ok<*ASGLU>
xlim([-7.5 7.5]);
ylim([-7.5 7.5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{x}''u_{z}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uxuz','x_D_',int2str(x_D), '.png'),'-dpng','-r600');
%% 6. Plot reystress_uxuy_av - <ux'uy'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uxuy_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{x}''u_{y}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uxuy','x_D_',int2str(x_D), '.png'),'-dpng','-r600');
%% 7. Plot reystress_uxuy_av - <uy'uz'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_uyuz_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);

% zmin = (min(reystress_uyuz_av(:))); 
% zmax = (max(reystress_uyuz_av(:)));
% zinc = (zmax - zmin) / 40;
% zlevs = zmin:zinc:zmax;

% axis equal;
set(h1,'Linecolor','none');
% set(h1,'levellist',zlevs);

colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle u_{y}''u_{z}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('uyuz','x_D_',int2str(x_D), '.png'),'-dpng','-r600');

%% 8. Plot reystress_rhouz_av - <rho'uz'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densuz_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle \rho''u_{z}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('rhouz','x_D_',int2str(x_D), '.png'),'-dpng','-r600');

%% 9. Plot reystress_rhouz_av - <rho'ux'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densux_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle \rho''u_{x}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('rhoux','x_D_',int2str(x_D), '.png'),'-dpng','-r600');

%% 10. Plot reystress_rhouz_av - <rho'uy'>  
figure;
[C,h1,x,y] = polarcont(rc, theta, reystress_densuy_av, 10); %#ok<*ASGLU>
xlim([-5 5]);
ylim([-5 5]);
% axis equal;
set(h1,'Linecolor','none');
colorbar;
hXLabel = xlabel('$y/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$z/D$','interpreter','latex','fontsize',20);
hTitle = title('$\langle \rho''u_{y}''\rangle$','interpreter','latex','fontsize',20);

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('rhouy','x_D_',int2str(x_D), '.png'),'-dpng','-r600');

close all;
