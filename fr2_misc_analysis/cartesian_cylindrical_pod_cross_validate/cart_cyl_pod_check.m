%% Name - Sheel Nidhan
%  Date - 24th February 2020

clear; clc;
%% Define a matrix in cylindrical coordinate with 2 grid points and same format as used for SPOD code

a = linspace(1,6,6)' + 1i;
b = linspace(7,12,6)' - 1i;
% c = linspace(13,18,6)';

A = [a  b];

%%  Transform each row to get the values in cartesian coordinate

theta1 = (pi/180)*30;
theta2 = (pi/180)*60;


B(1,:) = A(1,:)*cos(theta1) - A(3,:)*sin(theta1);
B(3,:) = A(1,:)*sin(theta1) + A(3,:)*cos(theta1);

B(2,:) = A(2,:)*cos(theta2) - A(4,:)*sin(theta2);
B(4,:) = A(2,:)*sin(theta2) + A(4,:)*cos(theta2);

B(5:6,:) = A(5:6,:);

%% Eigenvalues and eigenvector of the cross-spectral density

[eigmodeB, eigvB] = eig(B'*B);
[eigmodeA, eigvA] = eig(A'*A);


