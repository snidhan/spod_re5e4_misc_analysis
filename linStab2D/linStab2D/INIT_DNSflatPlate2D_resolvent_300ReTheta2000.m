clear variables
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');  set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultFigureRenderer','painters')
set(groot, 'defaultFigureColor',[1 1 1])
addpath('aux');
keep_L0     = false;     % keep_L0 must exist for pseudo-spectrum calculation, value is not really checked
global Nfilter Nz Nr

% 300 < Re_\theta < 2000
y1              = 4.4354;
x0              = 55.4149-1;
x1              = 150.4335;
z1              = 4*pi;
caseName        = 'RSEM';
calcMethod      = 'frequency_response'; %'input_output'%, 'timestepper', 'iterate_optimal', 'adjoint', else direct LU (Arnoldi)
root_data       = 'baseflows';
eas_file_mean   = 'mean_flow_all_mittel-z_ReTheta300-2000';
eas_file_digits = 5;
dt              = 0.0067*10*6;
dz              = 0.0491;
nz              = 128;
ATTRLEN         = 10;
UDEFLEN         = 20;

% save folder
saveFolder      = 'results';

% Independant coefficients
T_0             = 288.15;

% run frequency and wavenumber grid for frequency-wavenumber diagramm
freq_vec        = 0.1%linspace(0.01,1.0,30);
kz_vec          = 5%0:0.5:15;
checkLock       = false;
saveSoln        = false;
saveEigVecs     = false;
vizBaseAndRes   = true;
% checkLock       = true;
% saveSoln        = true;
% saveEigVecs     = true;
% vizBaseAndRes   = false;

for kz = kz_vec
    
    % Base-flow treatment
    x_LES_inlet     = 0.0;                          % start shortly after LES inlet since u_x(r,x<0.07) is dicont. in d/dr
    Nsmooth_r       = 5;
    Nsmooth_z       = 5;
    
    % Domain & grid
    FDorder         = 4;
    % r
    Nr              = 95;
    Nr_jet          = 35;
    Nr_shearl       = 35;
    r_farf          = y1+4;
    r_sponge        = r_farf-2;
    r_shearl1       = 0.5;                        % set to zero for equidistant mesh in r
    r_shearl2       = 2;                          % set to zero for equidistant mesh in r
    Nr_smooth       = 50;
    % z
    z_segm          = [x0-1 x0 x0+3 x1-0.5 x1 x1+2];
    Nz_segm         = [30 75 350 40 60];
    Nz_smooth       = 25;
    
    % BCs
    outletBC        = 'zero_gradient';                 % 'zero_gradient', else zero value
    farfieldBC      = 'zero_gradient';                 % 'zero_gradient', else zero value
    inletBC         = 'zero_gradient';                 % 'zero_gradient_pipe', else zero value
    bottomBC        = 'wall';
    
    % Sponge & damping
    spongeType      = 'poly5'; sponge_pow = 2;      % 'poly5' for fifth order smooth polynomial or 'pow' + sponge_pow for ((x-x0)/(x-x_max))^sponge_pow
    aSponge_left    = 3;                            % exponential grid stretching factor
    aSponge_right   = 3;
    aSponge_top     = 3;
    spongeEps_left  = 1e-5;
    spongeEps_right = 1e-4;
    spongeEps_top   = 1e-5;
    
    % Arnoldi settings
    noEigs          = 5;
    SIGMA_vec       = freq_vec*2*pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Input-output specifics
    %%%%%%%%%%%%%%%%%%%%%%%%
    % dq/dt = L*q + B*u
    % y = C*q (only for iput/output, NOT for frequency response)
    % u: input, y:output
    %%%%%%%%%%%%%%%%%%%%%%%%
    filter_f      = true; % filter forcing
    
    % window for input
    forcingType   = 'strip';
    INPUT_z_min   = z_segm(2)+1;
    INPUT_z_max   = INPUT_z_min+3;
    INPUT_r_min   = 1e-10;
    INPUT_r_max   = 3;
    OUTPUT_z_min  = z_segm(2);
    OUTPUT_z_max  = z_segm(end-1);
    OUTPUT_r_min  = 1e-10;
    OUTPUT_r_max  = INPUT_r_max;%r_sponge;
    nSmooth_B_x   = 0;
    nSmooth_B_y   = 5;
    nSmooth_C_x   = 5;
    nSmooth_C_y   = 5;
    
    % run LST
    Ncalcs          = length(SIGMA_vec);
    calcCount       = 1;
    custom_comment  = [caseName '_' forcingType];
    DNSflatPlate2D
    
end