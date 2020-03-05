clear; clc; close all;
%%

k1 = 2; %how many times you want wave to oscillate in x-dir
k2 = 10; %how many times you want wave to oscillate in y-dir
[X,Y] = meshgrid(linspace(-pi,pi));
% [X, Y] = meshgrid(linspace(-6,6,600));
Z = real((exp(sqrt(-1)*(k1*X + k2*Y))) + (exp(sqrt(-1)*(-k1*X - k2*Y))));

contourf(X,Y,Z);

%% FFT2 on the uniform data set 

% interpolated_fft_mode_padded = padarray(interpolated_fft_mode, [2 2], 0, 'both');
FFTed_mode = fft2(Z);
abs_FFTed_mode = abs(FFTed_mode);
%% FFT parameters 

L = size(FFTed_mode,1);

sampling_frequency = 1/((X(1,end) - X(end,1))/100);
sampling_period = 1/sampling_frequency;
kx = sampling_frequency*(-L/2+1:L/2)/L;

sampling_frequency = 1/((Y(end,1) - Y(1,end))/100);
sampling_period = 1/sampling_frequency;
ky = sampling_frequency*(-L/2+1:L/2)/L;


[KX KY] = meshgrid(kx, ky);

%% Manipulating the FFT mode calculation

FFTed_mode_rearranged(:,1:L/2-1) = FFTed_mode(:,L/2+2:L);
FFTed_mode_rearranged(:,L/2:L)   = FFTed_mode(:,1:L/2+1);

FFTed_mode_rearranged_2(1:L/2-1,:) = FFTed_mode_rearranged(L/2+2:L,:);
FFTed_mode_rearranged_2(L/2:L,:)   = FFTed_mode_rearranged(1:L/2+1,:);


%% Plot
figure;
[C,h] = contourf(2*pi*KX, 2*pi*KY, abs(FFTed_mode_rearranged_2), 10);
axis equal;
set(h,'Linecolor','none');
% set(h,'LevelList',-1:0.01:1)
colorbar

%% 1D FFT

k1 = 2;
X = linspace(-pi,pi);
Z = cos(k1*X);
plot(X,Z, 'linewidth',2);
FFT_Z = fft(Z);

L = size(X,2);
sampling_frequency = 1/((X(end)-X(1))/100);
sampling_period = 1/sampling_frequency;
kx = sampling_frequency*(-L/2+1:L/2)/L;

FFTed_mode_rearranged(1:L/2-1) = FFT_Z(L/2+2:L);
FFTed_mode_rearranged(L/2:L)   = FFT_Z(1:L/2+1);


% P2 = abs(FFT_Z/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

plot(2*pi*kx, abs(FFTed_mode_rearranged)/L, 'k-', 'Linewidth', 2);