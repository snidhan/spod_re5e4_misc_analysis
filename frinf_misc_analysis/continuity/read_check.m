
fid = fopen('wp_02000100_100.res');
h = fread(fid,0,'*uint64'); % May need adjusting
a = fread(fid, 91848, '*double');

for j = 1:258
    for i = 1:356
        w(i,j) = a((j-1)*356 + i, 1);
    end
end


%% Reading the grid files

numvar = 3;   % numvar = 3 (only velocity is used for kernel); 
              % numvar = 4 (velocity and density is used for kernel)

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid
%fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/fr2/x1_grid.in');   %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end


%% Creating polar field
cm_viridis = viridis(1000);
theta = linspace(0,2*pi,256)';

figure;
[C,h] = polarcont(rc,theta,w(1:354, 2:257),5);
set(h,'LineColor','none')
colormap('viridis');
colorbar;
axis equal
