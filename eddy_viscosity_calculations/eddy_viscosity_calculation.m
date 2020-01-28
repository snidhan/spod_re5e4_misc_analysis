%% Author - Sheel Nidhan
%  Date   - 28th October 2019

%% Parameters for the statistical files

clear; clc;
close all;

load('./similarity_w.mat');    % Loading the streamwise velocity file
load('./similarity_uxur.mat'); % Loading the uxur file
loc_planes = [5; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120];


%% Importing Karu's datafiles

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_TKE.dat';
LK_TKE   = importdata(filename);
LK_TKE_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_TKE(:,1)));
    LK_TKE_loc_planes(i,2) = LK_TKE(idx,4);
    LK_TKE_loc_planes(i,1) = LK_TKE(idx,1);

end

filename = '/home/sheel/Dropbox/spod_re5e4/statistical_results/Half_length_zwhazi_WMEAN.dat';
LK_mean   = importdata(filename);
LK_mean_loc_planes = zeros(22,2);
for i = 1:size(loc_planes)
    [val idx] = min(abs(loc_planes(i,1)-LK_mean(:,1)));
    LK_mean_loc_planes(i,2) = LK_mean(idx,4);
    LK_mean_loc_planes(i,1) = LK_mean(idx,1);
end

%% Performing the eddy viscosity analysis at different locations

dudr = zeros(size(reystress_uw_1d,1)-2, size(reystress_uw_1d,2));

for i = 1:size(reystress_uw_1d,2)
    U(:,i)    = mean_w_1d(:,i);
    Ud(:,i)   = U(end,i) - U(:,i); %#ok<*SAGROW>
    uxur(:,i) = reystress_uw_1d(:,i);
    
%     h1 = plot(rc, Ud, 'ko', 'Linewidth',2);
%     hold on;
%     h2 = plot(rc, -35*uxur, 'ro', 'Linewidth',2);   
%     
    % calculating the du/dr at each radial location
    count = 1;
    for j = 2:size(Ud,1)-1
        dudr(count,i) = (Ud(j+1,1) - Ud(j-1,1))/(rc(j+1,1) - rc(j-1,1));
        count = count + 1;
    end
    
end

%% Plotting the values of nu_t

smoothed_dudr = zeros(size(reystress_uw_1d,1)-2, size(reystress_uw_1d,2));

for i = 1:size(reystress_uw_1d,2)
    disp(i);
    smoothed_dudr(:,i) = movmean(dudr(:,i),20);
    %h3 = plot(rc(2:end-1), -smoothed_dudr, 'ro', 'Linewidth', 2);
    nu_t(:,i) = uxur(40:end-1,i)./smoothed_dudr(39:end,i);
    rc_nu_t = rc(40:end-1);
end

close all;
x0=0;
y0=0;
width=20;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


lineStyles = maxdistcolor(12,@srgb_to_Jab);

count = 1;
for i = 1:2:20
    plot(rc_nu_t/LK_mean_loc_planes(i,2), nu_t(:,i), 'Color', lineStyles(count,:), 'Linewidth', 2);
    count = count + 1;
    hold on;
end

for i = 21:22
    plot(rc_nu_t/LK_mean_loc_planes(i,2), nu_t(:,i), 'Color', lineStyles(count,:), 'Linewidth', 2);
    count = count + 1;
    hold on;
end

% xlim([0 2]);
% ylim([-2 0.2]);

ax = gca;
ax.FontSize = 20;

Legend = cell(12,1);
Legend{1} = '$x/D = 5$';
Legend{2} = '$x/D = 10$';
Legend{3} = '$x/D = 15$';
Legend{4} = '$x/D = 20$';
Legend{5} = '$x/D = 25$';
Legend{6} = '$x/D = 30$';
Legend{7} = '$x/D = 35$';
Legend{8} = '$x/D = 40$';
Legend{9} = '$x/D = 45$';
Legend{10} = '$x/D = 50$';
Legend{11} = '$x/D = 55$';
Legend{12} = '$x/D = 60$';


hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$\nu_{t}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.10, 0.45, 0]);
%hTitle  = title('(a)','interpreter','latex','fontsize', 20, 'Units', 'normalized', 'Position', [-0.05, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 20;
hLegend.FontWeight = 'bold';
% hLegend.Position = [0 0 1 1];
hLegend.Location = 'northeast';
