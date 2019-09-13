%% Name - Sheel Nidhan
%  Date - 10th September, 2019
%  Plots for the statistical results for prf_spod_re5e4_frinf paper

%% Scaling of the defect velocity as a function of x/D

filename = './files/Defect_centerline.dat';
udefect  = importdata(filename);

% Scaling till x/D = 50
log_spod_m1_st0136_50 = log_spod_m1_st0136(1:10,1);
log_x_sampled_50 = log_x_sampled(1:10,1);
[coeffs_50, S_50] = polyfit(log_x_sampled_50, log_spod_m1_st0136_50, 1);
y_50_fitted = polyval(coeffs_50, log_x_sampled_50);
CI_50 = polyparci(coeffs_50, S_50, 0.99);