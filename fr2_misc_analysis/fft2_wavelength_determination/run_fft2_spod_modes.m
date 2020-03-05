%% Script to run the fft2_spod_modes.m for plotting the spectra

clc; clear;
close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% x = [10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100];
x = [20; 40; 60; 80; 100];
f_idx = [1; 2; 3; 4];
% Write out the files for plotting in paper

dirout = './';

for i = 1:size(x,1)
    disp(i);
    for j = 1:size(f_idx,1)
        disp(j);
        [max_kx, max_ky] = fft2_spod_modes(x(i,1), f_idx(j,1));
        wavelength(j).wavevector_kx = max_kx;
        wavelength(j).wavevector_ky = max_ky;
        wavelength(j).f_idx         = f_idx(j,1);
    end
    fft2_spod_mode(i).wavenumber = wavelength; %#ok<*SAGROW>
    fft2_spod_mode(i).x_D        = x(i,1);
end

save(strcat(dirout, 'fft2_eigmodes_x_D_20_40_60_80_100.mat'), 'fft2_spod_mode');