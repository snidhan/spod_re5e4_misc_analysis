%% Reading the CSV file of Ud
figure
ud =importdata('ud_dairay.csv');
h2 = loglog(ud(2:end,1), ud(2:end,2), 'ro');
ud_x_loc = interp1(ud(1:end,1), ud(1:end,2), x_loc, 'spline');
hold on;
h3 = loglog(x_loc, ud_x_loc, 'bo');
%% Reading the CSV file of Ko/Ud^2

x_loc = [10; 20; 30; 40; 50; 60; 70; 80; 90; 100];
ko_ud_2 = importdata('tke_centerline_ud2.csv');
ko_ud_2 = flip(ko_ud_2);
ko = ko_ud_2(:,2).*(ud_x_loc.^2); 
% h1 = loglog(x_loc, ko_ud_2(:,2), 'ko');

%% Fit till x/D = 50

[coeffs_tke_x_D_50, S] = polyfit(log(x_loc(2:5)), log(ko(2:5).^0.5), 1)
[coeffs_tke_x_D_100, S] = polyfit(log(x_loc(6:10)), log(ko(6:10).^0.5), 1)
%% Plotting ko as a function of x/D
close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure;
ko = ko_ud_2(:,2).*(ud_x_loc.^2); 
h4 = loglog(x_loc, ko.^0.5, 'ko');
hold on;
y1 = 0.5*x_loc(2:4,1).^(coeffs_tke_x_D_50(1));
h5 = loglog(x_loc(2:4,1),y1, 'k--', 'Linewidth', 2);
h6 = text(25, 0.03,'$x^{-0.7}$','interpreter','latex','FontSize', 20);
y1 = 1*x_loc(6:10,1).^(coeffs_tke_x_D_100(1));
h5 = loglog(x_loc(6:10,1),y1, 'b--', 'Linewidth', 2);
h6 = text(65, 0.015,'$x^{-0.86}$','interpreter','latex','FontSize', 20);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K_{o}^{1/2}$','interpreter','latex','fontsize',20);

ax = gca;
ax.FontSize = 20; 

xlim([0 100]);
xticks([10 20 30 40 60 80 100]);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'ko_sqrt_dairay_etal_x_D.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'ko_sqrt_dairay_etal_x_D.eps'),'-depsc2','-r600');

%% Finding ko from different route

wake_width = importdata('wake_width.csv');
re_based_on_tke = importdata('re_local_basedon_tke.csv');

re_based_on_tke_x_loc = interp1(re_based_on_tke(:,1), re_based_on_tke(:,2), x_loc, 'spline');
wake_width_x_loc      = interp1(wake_width(:,1), wake_width(:,2), x_loc, 'spline');

ko_square_root = re_based_on_tke_x_loc./wake_width_x_loc*(1/5000);

figure
loglog(x_loc, ko_square_root, 'ko');
hold on
loglog(x_loc, ko.^0.5, 'r*');

%% Reading the CSV file of Ko/Ud^2

rxr_uo_2 = importdata('rxr_uo_square.csv');
rx_uo_2 = flip(rxr_uo_2);
rxr_max = rxr_uo_2(:,2).*(ud_x_loc.^2);
h1 = loglog(x_loc, rxr_max, 'ko');

%% Fit till x/D = 50

[coeffs_rxr_max_x_D_50, S] = polyfit(log(x_loc(2:5)), log(rxr_max(2:5)), 1)
[coeffs_rxr_max_x_D_100, S] = polyfit(log(x_loc(6:10)), log(rxr_max(6:10)), 1)


figure; 
h4 = loglog(x_loc, rxr_max, 'ko');
hold on;
y1 = 1.5*x_loc(2:4,1).^(coeffs_rxr_max_x_D_50(1));
h5 = loglog(x_loc(2:4,1),y1, 'k--', 'Linewidth', 2);
h6 = text(25, 0.0003,'$x^{-2}$','interpreter','latex','FontSize', 20);
y1 = 300*x_loc(6:10,1).^(coeffs_rxr_max_x_D_100(1));
h5 = loglog(x_loc(6:10,1),y1, 'b--', 'Linewidth', 2);
h6 = text(65, 0.000025,'$x^{-3.38}$','interpreter','latex','FontSize', 20);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$-\langle u''_{x}u''_{r} \rangle_{max}$','interpreter','latex','fontsize',20);

ax = gca;
ax.FontSize = 20; 

xlim([0 100]);
xticks([10 20 30 40 60 80 100]);

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat('./', 'rxr_dairay_etal_x_D.png'),'-dpng2','-r600');  
% print(gcf,strcat('./', 'rxr_dairay_etal_x_D.eps'),'-depsc2','-r600');


%% Ratio of Rxr/Ko as a function of x/D 

figure;
ratio = rxr_uo_2(:,2)./ko_ud_2(:,2);
h2 = plot(x_loc, ratio, 'ks');

ax = gca;
ax.FontSize = 20; 

xticks([10 20 30 40 50 60 70 80 90 100]);

hXLabel = xlabel('$x/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$R_{xr}/K_{o}$','interpreter','latex','fontsize',20);


set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat('./', 'rxr_ko_ratio_dairay_etal_x_D.png'),'-dpng2','-r600');  
print(gcf,strcat('./', 'rxr_ko_ratio_dairay_etal_x_D.eps'),'-depsc2','-r600');