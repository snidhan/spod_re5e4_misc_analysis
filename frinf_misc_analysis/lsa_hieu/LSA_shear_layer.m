clear all; close all; clc;

%% Define velocity profile and parameters
nu=1/24000.0;
kap=nu/1.0;
g=9.81; rho_0 = 1;
z=[10:-0.1:-10];z=z';
V=-0.5*tanh(z/0.5);


%% Specify range of Richardson number Ri and wavenumber k. 
Ri = 0.1; % 1st value of Ri

for Riloop=1:1
    rho = -Ri/9.81/2*tanh(z/0.5); % density profile for 2 layer
    B = -g*(rho)/rho_0; % convert density to buoyancy
    N2 = ddz(z)*B;  % compute squared buoyany frequency
    k = 0.1; % 1st value of wavenumber
    for kloop=1:30
        %Call SSF to compute the eigenvalue spectrum.
        %[sigs,w,b]=SSF(z,V,B,k,0,nu,kap,[1 1],[1 1],0);
        [sigs,w,b]=SSF(z,V,B,k,0,nu,kap,[0 0],[0 0],0);
        k_store(kloop,Riloop)=k;
        Ri_store(kloop,Riloop)=Ri;
        sigma_store(kloop,Riloop)=sigs(1);
        k = k+0.05; % increment to next value of wavenumber
    end
    Ri = Ri + 0.01;  % increment to next value of Ri
end
save('sigma.mat','k_store','Ri_store','sigma_store');
figure(1);
plot(k_store,real(sigma_store),'-*');
xlabel('kx');ylabel('growth rate');

