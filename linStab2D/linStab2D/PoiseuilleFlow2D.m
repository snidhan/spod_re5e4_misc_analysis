clear all
close all
clc
profile on
addpath aux

% Independant coefficients
Ma      = 0.1;                  % Mach number
Re      = 5000;                 % Reynolds number
Pr      = 0.71;                 % Prandtl number
kappa   = 1.4;                  % specific heat ratio
T_0     = 300;                  % reference temperature [K]

beta    = 1;                    % azimuthal wavenumber
alpha   = 1;                    % defines domain length

% Dependant coefficients
cv      = 1/(kappa*(kappa-1)*Ma^2);
c1      = (kappa-1)*Re*Pr*Ma^2;
c2      = kappa*Ma^2;

% Spatial discretization
DiffOrd = 4;

% Arnoldi settings
noEigs      = 5;
SIGMA       = 0.5043 - 0.3344i;
options.tol = 1e-6;

% Grid
Nz          = 150;
lz          = 2*pi/alpha;
Nr_axis     = 40;
Nr_wall     = 30;
r_wall      = 1;
r_cluster   = 0.85;

Nr          = Nr_axis + Nr_wall;
dr_avg      = r_wall/Nr;
r_1D        = [linspace(0,r_cluster-dr_avg/2,Nr_axis) linspace(r_cluster+dr_avg/2,r_wall,Nr_wall)]';
r_1D        = mvg(mvg(mvg(r_1D)));
z_1D        = linspace(0,lz,Nz)';
[z, r]      = meshgrid(z_1D, r_1D);

% Base-flow
RHO     = ones(Nr,Nz);
U       = zeros(Nr,Nz);
V       = zeros(Nr,Nz);
W       = (1 - r.^2);
T       = ones(Nr,Nz);

% Sutherland's Law
S1     = 110.4;                 % Sutherland temperature
T_ref  = 280.0;                 % reference temperature
mu_ref = 1.735e-5;              % reference dynamic viscosity

Tstar  = T*T_0;
% viscosity and thermal conductivity
mu_0            = mu_ref*(T_ref+S1)./(T_0  +S1).*(T_0  /T_ref).^(3/2);
mu              = mu_ref*(T_ref+S1)./(Tstar+S1).*(Tstar/T_ref).^(3/2);
dmudT           = -mu_ref.*(T_ref+S1)./(Tstar+S1).^2.*(Tstar./T_ref).^(3./2)+3./2.*mu_ref.*(T_ref+S1)./(Tstar+S1).*sqrt(Tstar./T_ref)./T_ref;
d2mudT2         = (2.*mu_ref.*(T_ref+S1)./(Tstar+S1).^3.*(Tstar./T_ref).^(3./2))-3.*mu_ref.*(T_ref+S1)./((Tstar+S1).^2).*sqrt((Tstar./T_ref))./T_ref+3./4.*mu_ref.*(T_ref+S1)./(Tstar+S1).*((Tstar./T_ref).^(-1./2))./(T_ref.^2);
MU              = mu          /mu_0;
dMUdT           = dmudT       /mu_0*T_0;
d2MUdT2         = d2mudT2     /mu_0*T_0^2;

tic
% Differentiation matrix
StencilSize = DiffOrd + 1;
if mod(DiffOrd,2) == 0
    % radial derivative
    Dr_1D       = zeros(Nr,Nr);
    i_minus     = ceil(DiffOrd/2);
    i_plus      = ceil(DiffOrd/2);
    for i=1:Nr
        if      i<=i_minus
            Dr_1D(i,1:StencilSize)       = fdweights(r_1D(i),r_1D(1:StencilSize),1);
        elseif  i>=Nr-i_plus+1
            Dr_1D(i,Nr-StencilSize+1:end)= fdweights(r_1D(i),r_1D(Nr-StencilSize+1:end),1);
        else
            Dr_1D(i,i-i_minus:i+i_plus)  = fdweights(r_1D(i),r_1D(i-i_minus:i+i_plus),1);
        end
    end   
    % axial derivative - CFD
    Stencil     = fdweights(z_1D(ceil(StencilSize/2)), z_1D(1:StencilSize), 1);
    col         = zeros(Nz,1);
    col(1:ceil(StencilSize/2))          = Stencil(ceil(StencilSize/2):end);
    col(end-floor(StencilSize/2)+1:end) = Stencil(1:floor(StencilSize/2));    
    Dz_1D        = toeplitz(col);
    for i=1:Nz, for j=1:Nz, if (i>j), Dz_1D(i,j) = - Dz_1D(i,j); end, end, end 
%     % axial derivative - spectral
%     [dummy, Dz_1D]        = fourdif(Nz,1);
%     Dz_1D                 = 1/(z_1D(end)-z_1D(1))*2*pi*Dz_1D;
else
    error('please use an even order differentiation scheme!')
end
time    = toc;
disp(['Elapsed time - Differentiation matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);

D2z_1D  = Dz_1D^2;
D2r_1D  = Dr_1D^2;
Dz      = sparse(kron(Dz_1D,eye(Nr)));
Dr      = sparse(kron(eye(Nz),Dr_1D));
DR      = kron(eye(5,5),Dr);
DZ      = kron(eye(5,5),Dz);
D2z     = sparse(kron(D2z_1D,eye(Nr)));
D2r     = sparse(kron(eye(Nz),D2r_1D));
D2R     = kron(eye(5,5),D2r);
D2Z     = kron(eye(5,5),D2z);

% Reshaping flowfield
NrNz= Nr*Nz;
RHO = reshape(RHO,   NrNz, 1);
U   = reshape(U,     NrNz, 1);
V   = reshape(V,     NrNz, 1);
W   = reshape(W,     NrNz, 1);
T   = reshape(T,     NrNz, 1);
MU  = reshape(MU,    NrNz, 1);
dmudT       = reshape(dMUdT, NrNz, 1);
d2mudT2     = reshape(d2MUdT2, NrNz, 1);

tic
% Base flow derivatives
dUdr    = Dr*U;     d2Udr2    = Dr*dUdr;
dVdr    = Dr*V;     d2Vdr2    = Dr*dVdr;
dWdr    = Dr*W;     d2Wdr2    = Dr*dWdr;
dRHOdr  = Dr*RHO;   d2RHOdr2  = Dr*dRHOdr;                         
dTdr    = Dr*T;     d2Tdr2    = Dr*dTdr;
dUdz    = Dz*U;     d2Udz2    = Dr*dUdz;
dVdz    = Dz*V;     d2Vdz2    = Dr*dVdz;
dWdz    = Dz*W;     d2Wdz2    = Dr*dWdz;
dRHOdz  = Dz*RHO;   d2RHOdz2  = Dr*dRHOdz;                         
dTdz    = Dz*T;     d2Tdz2    = Dr*dTdz;
d2Udrz  = Dz*dUdr;  d2Vdrz    = Dz*dVdr;
d2Wdrz  = Dz*dWdr;

time    = toc;
disp(['Elapsed time - Base flow derivatives: ' datestr(time/24/3600, 'HH:MM:SS')]);

% figure('name','Diff. test');
% subplot(3,1,1)
% contourf(z,r, reshape(d2Wdrz,Nr,Nz) ), xlabel('z'), ylabel('r'), title('W'), axis equal
% subplot(3,1,2)
% contourf(z,r, reshape(dWdr,Nr,Nz) ), xlabel('z'), ylabel('r'), title('dWdr'), axis equal
% subplot(3,1,3)
% contourf(z,r, reshape(dWdz,Nr,Nz) ), xlabel('z'), ylabel('r'), title('dWdz'), axis equal

tic
% % Coefficient matrices
Z = zeros(NrNz, 1);
I = ones(NrNz, 1);
R = reshape(r, NrNz, 1);
% Continuity
A0_11 = -(U + R .* dUdr + 1i .* beta .* V + R .* dWdz) ./ R;
A0_12 = -(RHO + R .* dRHOdr) ./ R;
A0_13 = -1i ./ R .* RHO .* beta;
A0_14 = -dRHOdz;
A0_15 = Z;
Ar_11 = -U;
Ar_12 = -RHO;
Ar_13 = Z;
Ar_14 = Z;
Ar_15 = Z;
Arr_11 = Z;
Arr_12 = Z;
Arr_13 = Z;
Arr_14 = Z;
Arr_15 = Z;
Az_11 = -W;
Az_12 = Z;
Az_13 = Z;
Az_14 = -RHO;
Az_15 = Z;
Azz_11 = Z;
Azz_12 = Z;
Azz_13 = Z;
Azz_14 = Z;
Azz_15 = Z;
Arz_11 = Z;
Arz_12 = Z;
Arz_13 = Z;
Arz_14 = Z;
Arz_15 = Z;
B_11 = 1i*I;
B_12 = Z;
B_13 = Z;
B_14 = Z;
B_15 = Z;
% Radial momentum
A0_21 = -(2 .* dUdr .* U .* R .* c2 + W .* dUdz .* R .* c2 + dWdz .* U .* R .* c2 + 1i .* beta .* U .* V .* c2 + dTdr .* R - V .^ 2 .* c2 + U .^ 2 .* c2) ./ R ./ c2;
A0_22 = -((3 .* MU .* beta .^ 2 + 6 .* RHO .* U .* Re .* R + 6 .* dRHOdr .* U .* Re .* R .^ 2 + 6 .* dUdr .* RHO .* Re .* R .^ 2 + 3 .* W .* dRHOdz .* Re .* R .^ 2 + 3 .* dWdz .* RHO .* Re .* R .^ 2 + 2 .* dmudT .* dTdr .* R + 3 .* 1i .* RHO .* beta .* V .* Re .* R + 4 .* MU) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_23 = -((2 .* 1i .* dmudT .* dTdr .* beta .* R - 6 .* RHO .* V .* Re .* R + 7 .* 1i .* MU .* beta + 3 .* 1i .* RHO .* U .* beta .* Re .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_24 = -RHO .* dUdz - dRHOdz .* U;
A0_25 = -((2 .* R .^ 2 .* c2 .* dWdz .* d2mudT2 .* dTdr + 2 .* R .* c2 .* U .* d2mudT2 .* dTdr - 4 .* R .^ 2 .* c2 .* dUdr .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dWdr .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c2 .* dUdz .* d2mudT2 .* dTdz + 3 .* 1i .* dmudT .* beta .* V .* c2 - 3 .* 1i .* dmudT .* beta .* dVdr .* R .* c2 + 3 .* Re .* R .^ 2 .* dRHOdr) ./ Re ./ R .^ 2 ./ c2) ./ 0.3e1;
Ar_21 = -(T + U .^ 2 .* c2) ./ c2;
Ar_22 = 0.2e1 ./ 0.3e1 .* (2 .* dmudT .* dTdr .* R - 3 .* RHO .* U .* Re .* R + 2 .* MU) ./ Re ./ R;
Ar_23 = 1i ./ Re ./ R .* MU .* beta ./ 0.3e1;
Ar_24 = 0.1e1 ./ Re .* dmudT .* dTdz;
Ar_25 = -((2 .* dWdz .* dmudT .* R .* c2 + 2 .* U .* dmudT .* c2 - 4 .* dUdr .* dmudT .* R .* c2 + 3 .* RHO .* Re .* R) ./ Re ./ R ./ c2) ./ 0.3e1;
Arr_21 = Z;
Arr_22 = 0.4e1 ./ 0.3e1 .* MU ./ Re;
Arr_23 = Z;
Arr_24 = Z;
Arr_25 = Z;
Az_21 = -W .* U;
Az_22 = -(W .* RHO .* Re - dmudT .* dTdz) ./ Re;
Az_23 = Z;
Az_24 = -((2 .* dmudT .* dTdr + 3 .* RHO .* U .* Re) ./ Re) ./ 0.3e1;
Az_25 = dmudT .* (dWdr + dUdz) ./ Re;
Azz_21 = Z;
Azz_22 = MU ./ Re;
Azz_23 = Z;
Azz_24 = Z;
Azz_25 = Z;
Arz_21 = Z;
Arz_22 = Z;
Arz_23 = Z;
Arz_24 = MU ./ Re ./ 0.3e1;
Arz_25 = Z;
B_21 = 1i .* U;
B_22 = 1i .* RHO;
B_23 = Z;
B_24 = Z;
B_25 = Z;
% Azimuthal momentum
A0_31 = -(1i .* beta .* V .^ 2 .* c2 + 2 .* V .* U .* c2 + W .* dVdz .* R .* c2 + dWdz .* V .* R .* c2 + dVdr .* U .* R .* c2 + V .* dUdr .* R .* c2 + 1i .* beta .* T) ./ R ./ c2;
A0_32 = ((-6 .* RHO .* V .* Re .* R - 3 .* V .* dRHOdr .* Re .* R .^ 2 - 3 .* dVdr .* RHO .* Re .* R .^ 2 + 7 .* 1i .* MU .* beta + 3 .* 1i .* dmudT .* dTdr .* beta .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_33 = -((3 .* MU + 3 .* dmudT .* dTdr .* R + 6 .* RHO .* U .* Re .* R + 3 .* W .* dRHOdz .* Re .* R .^ 2 + 3 .* dWdz .* RHO .* Re .* R .^ 2 + 3 .* dRHOdr .* U .* Re .* R .^ 2 + 3 .* RHO .* dUdr .* Re .* R .^ 2 + 4 .* MU .* beta .^ 2 + 6 .* 1i .* RHO .* V .* beta .* Re .* R) ./ Re ./ R .^ 2) ./ 0.3e1;
A0_34 = -(RHO .* dVdz .* Re .* R + dRHOdz .* V .* Re .* R - 1i .* dmudT .* dTdz .* beta) ./ Re ./ R;
A0_35 = -((3 .* 1i .* RHO .* beta .* Re .* R + 3 .* R .* c2 .* V .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dVdr .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c2 .* dVdz .* d2mudT2 .* dTdz - 4 .* 1i .* dmudT .* beta .* U .* c2 + 2 .* 1i .* dmudT .* beta .* dUdr .* R .* c2 + 2 .* 1i .* dmudT .* beta .* dWdz .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.3e1;
Ar_31 = -V .* U;
Ar_32 = -((3 .* RHO .* V .* Re .* R - 1i .* MU .* beta) ./ Re ./ R) ./ 0.3e1;
Ar_33 = (MU + dmudT .* dTdr .* R - RHO .* U .* Re .* R) ./ Re ./ R;
Ar_34 = Z;
Ar_35 = -dmudT .* (V - dVdr .* R) ./ Re ./ R;
Arr_31 = Z;
Arr_32 = Z;
Arr_33 = MU ./ Re;
Arr_34 = Z;
Arr_35 = Z;
Az_31 = -W .* V;
Az_32 = Z;
Az_33 = (dmudT .* dTdz - W .* RHO .* Re) ./ Re;
Az_34 = -((3 .* RHO .* V .* Re .* R - 1i .* MU .* beta) ./ Re ./ R) ./ 0.3e1;
Az_35 = 0.1e1 ./ Re .* dVdz .* dmudT;
Azz_31 = Z;
Azz_32 = Z;
Azz_33 = MU ./ Re;
Azz_34 = Z;
Azz_35 = Z;
Arz_31 = Z;
Arz_32 = Z;
Arz_33 = Z;
Arz_34 = Z;
Arz_35 = Z;
B_31 = 1i .* V;
B_32 = Z;
B_33 = 1i .* RHO;
B_34 = Z;
B_35 = Z;
% Axial momentum
A0_41 = -(dTdz .* R + W .* U .* c2 + W .* dUdr .* R .* c2 + 2 .* dWdz .* W .* R .* c2 + dWdr .* U .* R .* c2 + 1i .* beta .* V .* W .* c2) ./ R ./ c2;
A0_42 = -((2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re + 3 .* dWdr .* RHO .* R .* Re + 3 .* W .* dRHOdr .* R .* Re) ./ R ./ Re) ./ 0.3e1;
A0_43 = -(1i .* beta .* (2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re) ./ R ./ Re) ./ 0.3e1;
A0_44 = -(MU .* beta .^ 2 + RHO .* U .* R .* Re + RHO .* dUdr .* R .^ 2 .* Re + 2 .* dWdz .* RHO .* R .^ 2 .* Re + dRHOdr .* U .* R .^ 2 .* Re + 2 .* dRHOdz .* W .* R .^ 2 .* Re + 1i .* RHO .* V .* beta .* R .* Re) ./ R .^ 2 ./ Re;
A0_45 = -((3 .* dRHOdz .* R .* Re + 2 .* dUdr .* d2mudT2 .* dTdz .* R .* c2 - 4 .* dWdz .* d2mudT2 .* dTdz .* R .* c2 - 3 .* dWdr .* d2mudT2 .* dTdr .* R .* c2 - 3 .* dUdz .* d2mudT2 .* dTdr .* R .* c2 + 2 .* U .* d2mudT2 .* dTdz .* c2 - 3 .* 1i .* dmudT .* beta .* dVdz .* c2) ./ R ./ Re ./ c2) ./ 0.3e1;
Ar_41 = -W .* U;
Ar_42 = -((2 .* dmudT .* dTdz + 3 .* W .* RHO .* Re) ./ Re) ./ 0.3e1;
Ar_43 = Z;
Ar_44 = -(-MU - dmudT .* dTdr .* R + RHO .* U .* R .* Re) ./ R ./ Re;
Ar_45 = dmudT .* (dWdr + dUdz) ./ Re;
Arr_41 = Z;
Arr_42 = Z;
Arr_43 = Z;
Arr_44 = MU ./ Re;
Arr_45 = Z;
Az_41 = -(T + W .^ 2 .* c2) ./ c2;
Az_42 = ((MU + 3 .* dmudT .* dTdr .* R) ./ R ./ Re) ./ 0.3e1;
Az_43 = 1i ./ R ./ Re .* MU .* beta ./ 0.3e1;
Az_44 = 0.2e1 ./ 0.3e1 .* (2 .* dmudT .* dTdz - 3 .* W .* RHO .* Re) ./ Re;
Az_45 = -((3 .* RHO .* R .* Re + 2 .* dUdr .* dmudT .* R .* c2 - 4 .* dWdz .* dmudT .* R .* c2 + 2 .* U .* dmudT .* c2) ./ R ./ Re ./ c2) ./ 0.3e1;
Azz_41 = Z;
Azz_42 = Z;
Azz_43 = Z;
Azz_44 = 0.4e1 ./ 0.3e1 .* MU ./ Re;
Azz_45 = Z;
Arz_41 = Z;
Arz_42 = MU ./ Re ./ 0.3e1;
Arz_43 = Z;
Arz_44 = Z;
Arz_45 = Z;
B_41 = 1i .* W;
B_42 = Z;
B_43 = Z;
B_44 = 1i .* RHO;
B_45 = Z;
% Energy
A0_51 = -((1i .* beta .* V .^ 3 .* c2 + 2 .* cv .* T .* U .* c2 + V .^ 2 .* dWdz .* R .* c2 + 3 .* U .^ 2 .* dUdr .* R .* c2 + 3 .* W .^ 2 .* dWdz .* R .* c2 + V .^ 2 .* dUdr .* R .* c2 + U .^ 2 .* dWdz .* R .* c2 + W .^ 2 .* dUdr .* R .* c2 + 2 .* 1i .* V .* beta .* T + 2 .* T .* U + U .^ 3 .* c2 + 2 .* 1i .* V .* beta .* cv .* T .* c2 + 2 .* cv .* dTdz .* W .* R .* c2 + 2 .* dUdz .* U .* W .* R .* c2 + 2 .* cv .* T .* dWdz .* R .* c2 + 2 .* dVdz .* V .* W .* R .* c2 + 2 .* dVdr .* V .* U .* R .* c2 + 2 .* dWdr .* U .* W .* R .* c2 + 2 .* cv .* T .* dUdr .* R .* c2 + 2 .* cv .* dTdr .* U .* R .* c2 + 1i .* V .* beta .* U .^ 2 .* c2 + 1i .* V .* beta .* W .^ 2 .* c2 + 2 .* dTdz .* W .* R + 2 .* T .* dUdr .* R + 2 .* dTdr .* U .* R + V .^ 2 .* U .* c2 + W .^ 2 .* U .* c2 + 2 .* T .* dWdz .* R) ./ R ./ c2) ./ 0.2e1;
A0_52 = ((-6 .* T .* RHO .* Re .* R - 8 .* MU .* dWdz .* c2 .* R + 6 .* MU .* d2Udz2 .* R .^ 2 .* c2 - 6 .* T .* dRHOdr .* Re .* R .^ 2 - 6 .* dTdr .* RHO .* Re .* R .^ 2 + 8 .* MU .* d2Udr2 .* R .^ 2 .* c2 + 2 .* MU .* d2Wdrz .* R .^ 2 .* c2 - 6 .* U .* MU .* beta .^ 2 .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R - 6 .* W .* dRHOdz .* U .* Re .* R .^ 2 .* c2 - 18 .* dUdr .* RHO .* U .* Re .* R .^ 2 .* c2 - 6 .* dWdz .* RHO .* U .* Re .* R .^ 2 .* c2 - 6 .* dUdz .* W .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dVdr .* V .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dWdr .* W .* RHO .* Re .* R .^ 2 .* c2 - 6 .* cv .* T .* dRHOdr .* Re .* R .^ 2 .* c2 - 6 .* cv .* dTdr .* RHO .* Re .* R .^ 2 .* c2 + 12 .* 1i .* dVdr .* MU .* beta .* R .* c2 - 9 .* U .^ 2 .* RHO .* Re .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 3 .* W .^ 2 .* RHO .* Re .* c2 .* R - 4 .* W .* dmudT .* dTdz .* c2 .* R - 8 .* U .* dmudT .* dTdr .* c2 .* R + 2 .* 1i .* V .* MU .* beta .* c2 - 3 .* W .^ 2 .* dRHOdr .* Re .* R .^ 2 .* c2 - 9 .* dRHOdr .* U .^ 2 .* Re .* R .^ 2 .* c2 - 3 .* V .^ 2 .* dRHOdr .* Re .* R .^ 2 .* c2 + 6 .* dUdz .* dmudT .* dTdz .* R .^ 2 .* c2 + 6 .* dWdr .* dmudT .* dTdz .* R .^ 2 .* c2 + 8 .* dUdr .* dmudT .* dTdr .* R .^ 2 .* c2 - 4 .* dWdz .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* 1i .* V .* dmudT .* dTdr .* beta .* R .* c2 - 6 .* 1i .* U .* RHO .* V .* beta .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_53 = -((-6 .* MU .* d2Vdr2 .* R .^ 2 .* c2 - 6 .* MU .* d2Vdz2 .* R .^ 2 .* c2 + 6 .* dVdr .* MU .* R .* c2 + 8 .* 1i .* dUdr .* MU .* beta .* R .* c2 + 3 .* 1i .* beta .* RHO .* W .^ 2 .* Re .* R .* c2 + 3 .* 1i .* beta .* RHO .* U .^ 2 .* Re .* R .* c2 + 6 .* 1i .* beta .* RHO .* cv .* T .* Re .* R .* c2 + 8 .* 1i .* beta .* MU .* dWdz .* R .* c2 - 2 .* 1i .* U .* MU .* beta .* c2 + 6 .* 1i .* beta .* RHO .* T .* Re .* R + 4 .* 1i .* W .* dmudT .* dTdz .* beta .* R .* c2 + 8 .* V .* MU .* beta .^ 2 .* c2 + 12 .* V .* dmudT .* dTdr .* R .* c2 - 6 .* dVdz .* dmudT .* dTdz .* R .^ 2 .* c2 - 6 .* dVdr .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* U .* RHO .* V .* Re .* R .* c2 + 6 .* dWdz .* RHO .* V .* Re .* R .^ 2 .* c2 + 6 .* U .* dRHOdr .* V .* Re .* R .^ 2 .* c2 + 6 .* W .* dRHOdz .* V .* Re .* R .^ 2 .* c2 + 6 .* dUdr .* RHO .* V .* Re .* R .^ 2 .* c2 + 6 .* dVdz .* W .* RHO .* Re .* R .^ 2 .* c2 + 6 .* dVdr .* U .* RHO .* Re .* R .^ 2 .* c2 + 4 .* 1i .* U .* dmudT .* dTdr .* beta .* R .* c2 + 9 .* 1i .* beta .* RHO .* V .^ 2 .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_54 = ((6 .* MU .* d2Wdr2 .* R .^ 2 .* c2 - 6 .* T .* dRHOdz .* Re .* R .^ 2 - 6 .* dTdz .* RHO .* Re .* R .^ 2 + 2 .* MU .* d2Udrz .* R .^ 2 .* c2 + 8 .* MU .* d2Wdz2 .* R .^ 2 .* c2 + 6 .* MU .* dWdr .* R .* c2 + 2 .* MU .* dUdz .* R .* c2 - 6 .* W .* MU .* beta .^ 2 .* c2 - 18 .* dWdz .* RHO .* W .* Re .* R .^ 2 .* c2 - 6 .* dUdr .* RHO .* W .* Re .* R .^ 2 .* c2 - 6 .* U .* dRHOdr .* W .* Re .* R .^ 2 .* c2 - 6 .* cv .* dTdz .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dUdz .* U .* RHO .* Re .* R .^ 2 .* c2 - 6 .* cv .* T .* dRHOdz .* Re .* R .^ 2 .* c2 - 6 .* dVdz .* V .* RHO .* Re .* R .^ 2 .* c2 - 6 .* dWdr .* U .* RHO .* Re .* R .^ 2 .* c2 - 6 .* U .* RHO .* W .* Re .* R .* c2 + 12 .* 1i .* beta .* MU .* dVdz .* R .* c2 - 9 .* dRHOdz .* W .^ 2 .* Re .* R .^ 2 .* c2 - 3 .* V .^ 2 .* dRHOdz .* Re .* R .^ 2 .* c2 - 3 .* U .^ 2 .* dRHOdz .* Re .* R .^ 2 .* c2 + 6 .* dUdz .* dmudT .* dTdr .* R .^ 2 .* c2 + 6 .* dWdr .* dmudT .* dTdr .* R .^ 2 .* c2 - 4 .* dUdr .* dmudT .* dTdz .* R .^ 2 .* c2 + 8 .* dWdz .* dmudT .* dTdz .* R .^ 2 .* c2 - 4 .* U .* dmudT .* dTdz .* R .* c2 + 6 .* 1i .* V .* dmudT .* dTdz .* beta .* R .* c2 - 6 .* 1i .* W .* RHO .* V .* beta .* Re .* R .* c2) ./ Re ./ R .^ 2 ./ c2) ./ 0.6e1;
A0_55 = -((2 .* R .* c1 .* c2 .* U .^ 2 .* d2mudT2 .* dTdr + 3 .* R .* c1 .* c2 .* V .^ 2 .* d2mudT2 .* dTdr - 3 .* Re .* R .^ 2 .* c2 .* d2mudT2 .* dTdr .^ 2 - 3 .* Re .* R .^ 2 .* c2 .* d2mudT2 .* dTdz .^ 2 + 3 .* Re .* R .* c1 .* U .* RHO + 3 .* Re .* R .^ 2 .* c1 .* dUdr .* RHO + 3 .* Re .* R .^ 2 .* c1 .* U .* dRHOdr + 3 .* Re .* R .^ 2 .* c1 .* dWdz .* RHO + 3 .* Re .* R .^ 2 .* c1 .* W .* dRHOdz + 3 .* Re .* R .* c1 .* c2 .* U .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* W .* dRHOdz .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* dWdz .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* dUdr .* RHO .* cv + 3 .* Re .* R .^ 2 .* c1 .* c2 .* U .* dRHOdr .* cv - 3 .* R .^ 2 .* c1 .* c2 .* dUdz .* W .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c1 .* c2 .* dWdr .* W .* d2mudT2 .* dTdr - 3 .* R .^ 2 .* c1 .* c2 .* dUdz .* U .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dWdr .* U .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dVdz .* V .* d2mudT2 .* dTdz - 3 .* R .^ 2 .* c1 .* c2 .* dVdr .* V .* d2mudT2 .* dTdr + 2 .* R .^ 2 .* c1 .* c2 .* dUdr .* W .* d2mudT2 .* dTdz - 4 .* R .^ 2 .* c1 .* c2 .* dWdz .* W .* d2mudT2 .* dTdz + 2 .* R .* c1 .* c2 .* U .* W .* d2mudT2 .* dTdz - 4 .* R .^ 2 .* c1 .* c2 .* dUdr .* U .* d2mudT2 .* dTdr + 2 .* R .^ 2 .* c1 .* c2 .* dWdz .* U .* d2mudT2 .* dTdr + 3 .* 1i .* V .* RHO .* beta .* Re .* R .* c1 - 1i .* U .* dmudT .* beta .* V .* c1 .* c2 + 3 .* MU .* beta .^ 2 .* Re .* c2 + 3 .* 1i .* V .* RHO .* cv .* beta .* Re .* R .* c1 .* c2 - 3 .* 1i .* U .* dmudT .* beta .* dVdr .* R .* c1 .* c2 - 3 .* 1i .* W .* dmudT .* beta .* dVdz .* R .* c1 .* c2 + 2 .* 1i .* V .* dmudT .* beta .* dUdr .* R .* c1 .* c2 + 2 .* 1i .* V .* dmudT .* beta .* dWdz .* R .* c1 .* c2) ./ Re ./ R .^ 2 ./ c1 ./ c2) ./ 0.3e1;
Ar_51 = -(U .* (V .^ 2 .* c2 + W .^ 2 .* c2 + 2 .* T + U .^ 2 .* c2 + 2 .* cv .* T .* c2) ./ c2) ./ 0.2e1;
Ar_52 = ((16 .* dUdr .* MU .* c2 .* R - 9 .* U .^ 2 .* RHO .* Re .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 3 .* W .^ 2 .* RHO .* Re .* c2 .* R - 6 .* T .* RHO .* Re .* R - 8 .* MU .* dWdz .* c2 .* R + 2 .* 1i .* V .* MU .* beta .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R - 4 .* W .* dmudT .* dTdz .* c2 .* R + 8 .* U .* dmudT .* dTdr .* c2 .* R) ./ Re ./ c2 ./ R) ./ 0.6e1;
Ar_53 = ((6 .* dVdr .* MU .* R - 3 .* V .* MU + 1i .* U .* MU .* beta + 3 .* V .* dmudT .* dTdr .* R - 3 .* U .* RHO .* V .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Ar_54 = (2 .* dWdr .* MU .* R + W .* MU + 2 .* MU .* dUdz .* R + W .* dmudT .* dTdr .* R + U .* dmudT .* dTdz .* R - W .* RHO .* U .* Re .* R) ./ Re ./ R;
Ar_55 = ((6 .* dmudT .* dTdr .* Re .* R .* c2 - 3 .* Re .* R .* c1 .* U .* RHO - 2 .* U .^ 2 .* dmudT .* c1 .* c2 - 3 .* V .^ 2 .* dmudT .* c1 .* c2 + 3 .* MU .* Re .* c2 + 3 .* dUdz .* W .* dmudT .* R .* c1 .* c2 + 3 .* dWdr .* W .* dmudT .* R .* c1 .* c2 - 3 .* Re .* R .* c1 .* c2 .* U .* RHO .* cv + 3 .* dVdr .* V .* dmudT .* R .* c1 .* c2 + 4 .* dUdr .* U .* dmudT .* R .* c1 .* c2 - 2 .* dWdz .* U .* dmudT .* R .* c1 .* c2) ./ Re ./ R ./ c1 ./ c2) ./ 0.3e1;
Arr_51 = Z;
Arr_52 = 0.4e1 ./ 0.3e1 ./ Re .* U .* MU;
Arr_53 = 0.1e1 ./ Re .* V .* MU;
Arr_54 = 0.1e1 ./ Re .* W .* MU;
Arr_55 = 0.1e1 ./ c1 .* MU;
Az_51 = -(W .* (V .^ 2 .* c2 + W .^ 2 .* c2 + 2 .* T + U .^ 2 .* c2 + 2 .* cv .* T .* c2) ./ c2) ./ 0.2e1;
Az_52 = ((6 .* MU .* dUdz .* R + W .* MU + 6 .* dWdr .* MU .* R + 3 .* W .* dmudT .* dTdr .* R + 3 .* U .* dmudT .* dTdz .* R - 3 .* W .* RHO .* U .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Az_53 = ((6 .* dVdz .* MU .* R + 1i .* W .* MU .* beta + 3 .* V .* dmudT .* dTdz .* R - 3 .* W .* RHO .* V .* Re .* R) ./ Re ./ R) ./ 0.3e1;
Az_54 = ((16 .* MU .* dWdz .* c2 .* R - 3 .* V .^ 2 .* RHO .* Re .* c2 .* R - 9 .* W .^ 2 .* RHO .* Re .* c2 .* R - 3 .* U .^ 2 .* RHO .* Re .* c2 .* R - 8 .* U .* MU .* c2 - 8 .* dUdr .* MU .* c2 .* R - 6 .* T .* RHO .* Re .* R + 2 .* 1i .* V .* MU .* beta .* c2 - 6 .* cv .* T .* RHO .* Re .* c2 .* R + 8 .* W .* dmudT .* dTdz .* c2 .* R - 4 .* U .* dmudT .* dTdr .* c2 .* R) ./ Re ./ c2 ./ R) ./ 0.6e1;
Az_55 = -((-6 .* dmudT .* dTdz .* Re .* R .* c2 + 3 .* W .* RHO .* Re .* R .* c1 - 3 .* dUdz .* U .* dmudT .* R .* c1 .* c2 - 3 .* dWdr .* U .* dmudT .* R .* c1 .* c2 - 3 .* dVdz .* V .* dmudT .* R .* c1 .* c2 + 3 .* W .* RHO .* cv .* Re .* R .* c1 .* c2 + 2 .* dUdr .* W .* dmudT .* R .* c1 .* c2 - 4 .* dWdz .* W .* dmudT .* R .* c1 .* c2 + 2 .* U .* W .* dmudT .* c1 .* c2) ./ Re ./ R ./ c1 ./ c2) ./ 0.3e1;
Azz_51 = Z;
Azz_52 = 0.1e1 ./ Re .* U .* MU;
Azz_53 = 0.1e1 ./ Re .* V .* MU;
Azz_54 = 0.4e1 ./ 0.3e1 ./ Re .* W .* MU;
Azz_55 = 0.1e1 ./ c1 .* MU;
Arz_51 = Z;
Arz_52 = 0.1e1 ./ Re .* W .* MU ./ 0.3e1;
Arz_53 = Z;
Arz_54 = 0.1e1 ./ Re .* U .* MU ./ 0.3e1;
Arz_55 = Z;
B_51 = (1i .* (V .^ 2 + W .^ 2 + U .^ 2 + 2 .* cv .* T)) ./ 0.2e1;
B_52 = 1i .* U .* RHO;
B_53 = 1i .* V .* RHO;
B_54 = 1i .* W .* RHO;
B_55 = 1i .* RHO .* cv;

A0 =    ...
    [   ...
    sparse(diag(A0_11)) sparse(diag(A0_12)) sparse(diag(A0_13)) sparse(diag(A0_14)) sparse(diag(A0_15))
    sparse(diag(A0_21)) sparse(diag(A0_22)) sparse(diag(A0_23)) sparse(diag(A0_24)) sparse(diag(A0_25))
    sparse(diag(A0_31)) sparse(diag(A0_32)) sparse(diag(A0_33)) sparse(diag(A0_34)) sparse(diag(A0_35))
    sparse(diag(A0_41)) sparse(diag(A0_42)) sparse(diag(A0_43)) sparse(diag(A0_44)) sparse(diag(A0_45))
    sparse(diag(A0_51)) sparse(diag(A0_52)) sparse(diag(A0_53)) sparse(diag(A0_54)) sparse(diag(A0_55))    
    ];

Ar =    ...
    [   ...
    sparse(diag(Ar_11)) sparse(diag(Ar_12)) sparse(diag(Ar_13)) sparse(diag(Ar_14)) sparse(diag(Ar_15))
    sparse(diag(Ar_21)) sparse(diag(Ar_22)) sparse(diag(Ar_23)) sparse(diag(Ar_24)) sparse(diag(Ar_25))
    sparse(diag(Ar_31)) sparse(diag(Ar_32)) sparse(diag(Ar_33)) sparse(diag(Ar_34)) sparse(diag(Ar_35))
    sparse(diag(Ar_41)) sparse(diag(Ar_42)) sparse(diag(Ar_43)) sparse(diag(Ar_44)) sparse(diag(Ar_45))
    sparse(diag(Ar_51)) sparse(diag(Ar_52)) sparse(diag(Ar_53)) sparse(diag(Ar_54)) sparse(diag(Ar_55))    
    ];

Az =    ...
    [   ...
    sparse(diag(Az_11)) sparse(diag(Az_12)) sparse(diag(Az_13)) sparse(diag(Az_14)) sparse(diag(Az_15))
    sparse(diag(Az_21)) sparse(diag(Az_22)) sparse(diag(Az_23)) sparse(diag(Az_24)) sparse(diag(Az_25))
    sparse(diag(Az_31)) sparse(diag(Az_32)) sparse(diag(Az_33)) sparse(diag(Az_34)) sparse(diag(Az_35))
    sparse(diag(Az_41)) sparse(diag(Az_42)) sparse(diag(Az_43)) sparse(diag(Az_44)) sparse(diag(Az_45))
    sparse(diag(Az_51)) sparse(diag(Az_52)) sparse(diag(Az_53)) sparse(diag(Az_54)) sparse(diag(Az_55))   
    ];

Arz =   ...
    [   ...
    sparse(diag(Arz_11)) sparse(diag(Arz_12)) sparse(diag(Arz_13)) sparse(diag(Arz_14)) sparse(diag(Arz_15))
    sparse(diag(Arz_21)) sparse(diag(Arz_22)) sparse(diag(Arz_23)) sparse(diag(Arz_24)) sparse(diag(Arz_25))
    sparse(diag(Arz_31)) sparse(diag(Arz_32)) sparse(diag(Arz_33)) sparse(diag(Arz_34)) sparse(diag(Arz_35))
    sparse(diag(Arz_41)) sparse(diag(Arz_42)) sparse(diag(Arz_43)) sparse(diag(Arz_44)) sparse(diag(Arz_45))
    sparse(diag(Arz_51)) sparse(diag(Arz_52)) sparse(diag(Arz_53)) sparse(diag(Arz_54)) sparse(diag(Arz_55))    
    ];

Arr =   ...
    [   ...
    sparse(diag(Arr_11)) sparse(diag(Arr_12)) sparse(diag(Arr_13)) sparse(diag(Arr_14)) sparse(diag(Arr_15))
    sparse(diag(Arr_21)) sparse(diag(Arr_22)) sparse(diag(Arr_23)) sparse(diag(Arr_24)) sparse(diag(Arr_25))
    sparse(diag(Arr_31)) sparse(diag(Arr_32)) sparse(diag(Arr_33)) sparse(diag(Arr_34)) sparse(diag(Arr_35))
    sparse(diag(Arr_41)) sparse(diag(Arr_42)) sparse(diag(Arr_43)) sparse(diag(Arr_44)) sparse(diag(Arr_45))
    sparse(diag(Arr_51)) sparse(diag(Arr_52)) sparse(diag(Arr_53)) sparse(diag(Arr_54)) sparse(diag(Arr_55))    
    ];

Azz =   ...
    [   ...
    sparse(diag(Azz_11)) sparse(diag(Azz_12)) sparse(diag(Azz_13)) sparse(diag(Azz_14)) sparse(diag(Azz_15))
    sparse(diag(Azz_21)) sparse(diag(Azz_22)) sparse(diag(Azz_23)) sparse(diag(Azz_24)) sparse(diag(Azz_25))
    sparse(diag(Azz_31)) sparse(diag(Azz_32)) sparse(diag(Azz_33)) sparse(diag(Azz_34)) sparse(diag(Azz_35))
    sparse(diag(Azz_41)) sparse(diag(Azz_42)) sparse(diag(Azz_43)) sparse(diag(Azz_44)) sparse(diag(Azz_45))
    sparse(diag(Azz_51)) sparse(diag(Azz_52)) sparse(diag(Azz_53)) sparse(diag(Azz_54)) sparse(diag(Azz_55))   
    ];

ZZ = sparse(zeros(NrNz,NrNz));

RHS =   ...
    [   ...
    sparse(diag(B_11))  ZZ                  ZZ                  ZZ                  ZZ
    ZZ                  sparse(diag(B_22))  ZZ                  ZZ                  ZZ
    ZZ                  ZZ                  sparse(diag(B_33))  ZZ                  ZZ
    ZZ                  ZZ                  ZZ                  sparse(diag(B_44))  ZZ
    ZZ                  ZZ                  ZZ                  ZZ                  sparse(diag(B_55))    
    ];

% clean up memory
clear A0_* Ar_* Az_* Arr_* Azz_* Arz_* B_* d* I Z ZZ
LHS     = A0 + Ar*DR + Az*DZ + Arr*D2R + Azz*D2Z + Arz*DR*DZ;
clear A0 Ar Az Arr Azz Arz D2R D2Z
time    = toc;
disp(['Elapsed time - Coefficient matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);

tic
%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%
% Row indices for primitive variables
li = 1 : Nr;
ri = li + NrNz - Nr;
bi = 1 : Nr : NrNz;
ti = bi + Nr - 1;

ti_rho  = ti;
bi_rho  = bi;
ri_rho  = ri;
li_rho  = li;

ti_u    = ti   + NrNz;
bi_u    = bi   + NrNz;
ri_u    = ri   + NrNz;
li_u    = li   + NrNz;

ti_v    = ti   + 2*NrNz;
bi_v    = bi   + 2*NrNz;
ri_v    = ri   + 2*NrNz;
li_v    = li   + 2*NrNz;

ti_w    = ti   + 3*NrNz;
bi_w    = bi   + 3*NrNz;
ri_w    = ri   + 3*NrNz;
li_w    = li   + 3*NrNz;

ti_T    = ti   + 4*NrNz;
bi_T    = bi   + 4*NrNz;
ri_T    = ri   + 4*NrNz;
li_T    = li   + 4*NrNz;

% Column indices for conserved variables
rho_j   =        1 :   NrNz;
u_j     =   NrNz+1 : 2*NrNz;
v_j     = 2*NrNz+1 : 3*NrNz;
w_j     = 3*NrNz+1 : 4*NrNz;
T_j     = 4*NrNz+1 : 5*NrNz;

% wall
rho_wall_i = ti_rho;
  u_wall_i = ti_u; 
  v_wall_i = ti_v; 
  w_wall_i = ti_w; 
  T_wall_i = ti_T;
% axis  
rho_axis_i = bi_rho;
  u_axis_i = bi_u; 
  v_axis_i = bi_v; 
  w_axis_i = bi_w; 
  T_axis_i = bi_T;

% force to have a non-zero diagonal  
LHS = LHS + speye(5*NrNz)*eps;
RHS = RHS + speye(5*NrNz)*eps;

% wall
LHS([u_wall_i v_wall_i w_wall_i T_wall_i], :)        = 0;
RHS([u_wall_i v_wall_i w_wall_i T_wall_i], :)        = 0;
for i = [u_wall_i v_wall_i w_wall_i T_wall_i]
    RHS(i, i)        = 1;
    LHS(i, i)        = 1;
end

% axis
if (beta == 0)
    % zero value: u
    LHS([u_axis_i v_axis_i], :)        = 0;
    RHS([u_axis_i v_axis_i], :)        = 0;
    for i = [u_axis_i v_axis_i]
        RHS(i, i)        = 1;
        LHS(i, i)        = 1;
    end
    %zero gradient: rho, w, T
    LHS([rho_axis_i w_axis_i T_axis_i], :) = 0;
    RHS([rho_axis_i w_axis_i T_axis_i], :) = 0;
    RHS(rho_axis_i, rho_j)   = DR(rho_axis_i, rho_j);
    RHS(  w_axis_i,   w_j)   = DR(  w_axis_i,   w_j);
    RHS(  T_axis_i,   T_j)   = DR(  T_axis_i,   T_j);
elseif (abs(beta) == 1)
    % zero value: w, rho, T
    LHS([w_axis_i rho_axis_i T_axis_i], :)        = 0;
    RHS([w_axis_i rho_axis_i T_axis_i], :)        = 0;
    for i = [w_axis_i rho_axis_i T_axis_i]
        RHS(i, i)        = 1;
        LHS(i, i)        = 1;
    end
    %zero gradient: u, v
    LHS([u_axis_i v_axis_i], :) = 0;
    RHS([u_axis_i v_axis_i], :) = 0;
    RHS(u_axis_i, u_j)   = DR(u_axis_i, u_j);
    RHS(v_axis_i, v_j)   = DR(v_axis_i, v_j); 
elseif (abs(beta) >= 2)
    %zero value: all
    LHS([u_axis_i v_axis_i w_axis_i rho_axis_i T_axis_i], :)        = 0;
    RHS([u_axis_i v_axis_i w_axis_i rho_axis_i T_axis_i], :)        = 0;
    for i = [u_axis_i v_axis_i w_axis_i rho_axis_i T_axis_i]
        RHS(i, i)        = 1;
        LHS(i, i)        = 1;
    end    
end
time    = toc;
disp(['Elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

%%%%%%%%%%%%%%%%%
% Solve EVP
%%%%%%%%%%%%%%%%%

tic
L0               = -RHS\LHS;
clear LHS RHS
whosL0           = whos('L0');
disp(['Matrix memory requirenment: ' num2str(whosL0.bytes/1048576), ' MB']);
opptions.disp    = 3;
[eigVecs, OMEGA] = eigs(L0, noEigs, SIGMA, options);
clear L0
OMEGA            = diag(OMEGA);
time = toc;
disp(['Elapsed time - Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);
clear RHS LHS

%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%
figure('name','spectrum');
hold on
for i=1:noEigs
    plot(real(OMEGA(i)),imag(OMEGA(i)),'go')
    text(real(OMEGA(i)),imag(OMEGA(i)),[' ' num2str(i)])
end
plot(real(SIGMA),imag(SIGMA),'rx');

nrows = floor(sqrt(noEigs));
ncols = ceil(noEigs/nrows);
figure('name', 'w''')
for i=1:noEigs
   subplot(nrows,ncols,i)
   w = reshape(eigVecs(w_j,i), Nr, Nz);
   contourf(z, r, real(w), 'edgecolor','none'), xlabel('z'), ylabel('r'), axis equal, axis tight, title(['mode ' num2str(i)]) 
end

% while 1
%     i = input('Eigenvector #?');
%     Nr_interp          = 200;
%     Nz_interp          = 50;
%     Nphi_interp        = 90;
%     
%     w = reshape(eigVecs(w_j,i), Nr, Nz);
%     
%     phi_1D                  = linspace(0,2*pi,Nphi_interp);
%     [rmesh,phimesh,zmesh]   = meshgrid(r_1D,phi_1D,z_1D);
%     [xmesh,ymesh,zmesh]     = pol2cart(phimesh,rmesh,zmesh);
%     
%     w = w./(norm(w(:),Inf));
%     w = permute(repmat(w,[1 1 Nphi_interp]),[3 1 2]);
%     w = w.*exp(1).^(1i*beta.*phimesh);
%     w = real(w);
%     
%     [xg, yg, zg]    = meshgrid(linspace(-r_1D(end),r_1D(end),Nr_interp),linspace(-r_1D(end),r_1D(end),Nr_interp),linspace(z_1D(1),z_1D(end),Nz_interp));
%     gdata           = griddata(xmesh,ymesh,zmesh,w,xg,yg,zg);
%     gdata(isnan(gdata)) = 0;
%     
%     figure('name',['+/- 0.2*max(real(w'')) for mode ' num2str(i)])
%     % positive contour
%     p = patch(isosurface(xg,yg,zg,gdata,0.2));
%     isonormals(xg,yg,zg,gdata,p)
%     set(p,'FaceColor','yellow','EdgeColor','none');
%     hold on
%     % negative contour
%     p = patch(isosurface(xg,yg,zg,gdata,-0.2));
%     isonormals(xg,yg,zg,gdata,p)
%     set(p,'FaceColor','blue','EdgeColor','none');
%     
%     daspect([1,1,1])
%     view(3); axis tight
%     camlight
%     lighting gouraud
%     
%     drawnow
% end
% 
% profile viewer
% p = profile('info');
















