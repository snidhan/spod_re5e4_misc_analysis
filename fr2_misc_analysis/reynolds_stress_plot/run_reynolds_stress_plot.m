%% Script to run the spod_spectrum_ploot_imag_data.m for plotting the spectra


clc; clear;
close all;

% x = [10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100];
x = [10; 30; 40; 50; 70; 90; 100];

% Write out the files for plotting in paper

dirout = './';

count = 1;
for i = 1:size(x,1)
    reynolds_stress_plot(x(i));
end