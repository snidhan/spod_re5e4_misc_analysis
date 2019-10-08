%% Written by Sheel Nidhan
%  Comparing data of Dairay et al. 2015

%% Best fitting data of Dairay et al. 2015

close all;

A = importdata('udefect_dairay_et_al_2015.csv');
B = importdata('k_to_square_ud_centerline.csv');

x_loc = A(:,1);
ud    = A(:,2);

plot(x_loc, ud, 'ko');

x_candidate = [10; 20; 30; 40; 50; 60; 70; 80; 90; 100];

ud_interpolated = spline(x_loc, ud, x_candidate);

hold on;

plot(x_candidate, ud_interpolated, 'ro');

k_centerline_normalized = B(:,2);

k_centerline_normalized = k_centerline_normalized.*(ud_interpolated).^2;

figure;
loglog(x_candidate, k_centerline_normalized.^0.5,'ko');
hold on;

% Scaling
[coeffs, S] = polyfit(log(x_candidate(2:5,1)), log(k_centerline_normalized(2:5,1).^0.5), 1);
best_fit_k_cent = polyval(coeffs, log(x_candidate));
loglog(x_candidate, exp(best_fit_k_cent), 'k-');