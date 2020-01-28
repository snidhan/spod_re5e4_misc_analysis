clear all
% close all
clc
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');  set(groot, 'defaultTextInterpreter','latex'); 
set(groot, 'defaultFigureRenderer','painters')
set(groot, 'defaultFigureColor',[1 1 1])
addpath('../aux');
keep_L0     = false;     % keep_L0 must exist for pseudo-spectrum calculation, value is not really checked
global Nfilter Nz Nr

calcMethod  = 'frequency_response'; %'input_output'%, 'timestepper', 'iterate_optimal', 'adjoint', else direct LU (Arnoldi)
baseState   = 'M09';

% Independant coefficients
T_0             = 288.15;
rho_0           = 1.225;
mu_0            = 1.789e-5;
p_0             = 1.013e5;
nu_0            = mu_0/rho_0;
Rgas            = 287;
kappa           = 1.4;
Pr              = 0.7;
c               = sqrt(kappa*Rgas*T_0);
Ma              = 0.9;
u_0             = Ma*c;

% discrete LES frequencies
nfft            = 256;
ovlap           = 0.5;
dt              = 0.2;
St              = [0:nfft/2]/Ma/nfft/dt;
St              = 0.39;%St(findnearest(St,0.2):findnearest(St,1.5));

% azimuthal wave number
m               = 0;

% eddy viscosity model
Re              = 1e4;
useEddyVisc     = false;
eddyVisc        = mu_0*100;                     % number or 'LES' for Boussinesq approx.
PrT             = 0.9;
Nsmooth_nuT     = 100;

% Base-flow treatment
x_LES_inlet     = 0.0;                          % start shortly after LES inlet since u_x(r,x<0.07) is dicont. in d/dr
Nsmooth_r       = 5;
Nsmooth_z       = 5;

% Domain & grid
% r
Nr              = 195;
Nr_jet          = 50;
Nr_shearl       = 70;
r_farf          = 6.0923+1;
r_sponge        = r_farf-1;
r_shearl1       = 0.3;                          % set to zero for equidistant mesh in r
r_shearl2       = 0.7;                          % set to zero for equidistant mesh in r
Nr_smooth       = 15;
% z
Nz              = 950;
z1              = -1;                           % Aaron's St~0.4 marginal instability @ z1 = 3                 
z2              = 30.0087+1;
Nz_core         = 400;
Nz_smooth       = 250;
z_core          = 5;                            % set to zero for equidistant mesh in z
z_sponge1       = z1+1;
z_sponge2       = z2-1;

% BCs
pipeBC          = false;
outletBC        = 'zero_gradient';              % 'zero_gradient', else zero value
farfieldBC      = 'zero_value';                 % 'zero_gradient', else zero value
inletBC         = 'zero_value';                 % 'zero_gradient_pipe', else zero value

% Sponge & damping
spongeType      = 'poly5'; sponge_pow = 2;      % 'poly5' for fifth order smooth polynomial or 'pow' + sponge_pow for ((x-x0)/(x-x_max))^sponge_pow
aSponge_left    = 0;                            % exponential grid stretching factor
aSponge_right   = 2.25;
aSponge_top     = 2;
spongeEps_left  = 50;
spongeEps_right = 2.5;
spongeEps_top   = 0.5;

% Arnoldi settings
noEigs          = 3;
SIGMA_vec       = St*2*pi;

% Save & visualization options
saveFolder      = 'results';
% saveSoln        = true;
% saveEigVecs     = true;
% vizBaseAndRes   = false;
saveSoln        = false;
saveEigVecs     = false;
vizBaseAndRes   = true;

% time stepper parameters
Dt              = 1;
Nfilter         = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%
% Input-output specifics
%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*u
% y = C*q (only for iput/output, NOT for frequency response)
% u: input, y:output
%%%%%%%%%%%%%%%%%%%%%%%%
filter_f      = true % filter forcing

% forcingType   = 'global'
forcingType   = 'global_momentum'
% forcingType   = 'diag_TKE'
% forcingType   = 'farfield'
% forcingType   = 'core'
% forcingType   = 'jet'
% forcingType   = 'shearlayer'
% responseType    = 'global'
responseType    = 'pressure'

% window for input
INPUT_z_min   = z_sponge1;
INPUT_z_max   = z_sponge2;
INPUT_r_min   = 0.0;
INPUT_r_max   = r_sponge;

OUTPUT_z_min  = -Inf;
OUTPUT_z_max  = Inf;
OUTPUT_r_min  = 0;
OUTPUT_r_max  = Inf;

% run LST
Ncalcs          = length(SIGMA_vec);
calcCount       = 1;
custom_comment  = [baseState '_' forcingType]% '_noFilter']
LESmeanflowJet2D