% OTS, 2018
if ~exist('Nsmooth_r','var'),   Nsmooth_r   = 0;        end
if ~exist('Nsmooth_z','var'),   Nsmooth_z   = 0;        end
if ~exist('x_LES_inlet','var'), x_LES_inlet = 0;        end
if ~exist('spongeType','var'),  spongeType  = 'poly5';  end
if ~exist('checkLock','var'),   checkLock   = false;    end
if ~exist('useEddyVisc','var'), useEddyVisc = false;    end; if useEddyVisc, useEddyVisc_string=' with turbulent eddy viscosity'; else useEddyVisc_string=' using the molecular viscosity'; end
if ~exist('Ncalcs','var'),      Ncalcs      = 1;        end
if ~exist('calcCount','var'),   calcCount   = 1;        end
if ~exist('pipeBC','var'),      pipeBC      = 0;        end
if ~exist('forcingType','var'), forcingType = 'global'; end
if ~exist('responseType','var'),responseType= 'global'; end
if ~exist('baseState','var'),   baseState   = 'M09';    end
if ~exist('filter_f','var'),    filter_f    = false;    end
if ~exist('FDorder','var'),     FDorder     = 4;        end

NrNz= Nr*Nz;

disp(['--> calculating for Re=' num2str(Re) ', m=' num2str(m) ', Nr=' num2str(Nr) ', Nz=' num2str(Nz) useEddyVisc_string]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiation Matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if Nz~=0&&Nr~=0
    % temporary equidistant mesh
    dr              = r_farf/(Nr-0.5);
    r_1D_equidist   = linspace(dr/2,r_farf,Nr)';
    z_1D_equidist   = linspace(z1,z2,Nz)';
    dz              = z_1D_equidist(2)-z_1D_equidist(1);
    % fourth order summation by parts after Mattsson & Nordstroem (2004, p. 524)
    DiffOrd = 'SBP42';
    tic
    
    % construct transformed coordinates
    if (r_shearl1~=0&&r_shearl2~=0)
        dr_jet          = r_shearl1/(Nr_jet-0.5);
        dr_shearl       = (r_shearl2-r_shearl1)/(Nr_shearl-1);
        dr_farf         = (r_farf-r_shearl2)/(Nr-Nr_jet-Nr_shearl-1);
        dr_avg1         = (dr_jet+dr_shearl)/2;
        dr_avg2         = (dr_shearl+dr_farf)/2;
        r_1D            = [linspace(dr_jet/2,r_shearl1-dr_avg1/2,Nr_jet) linspace(r_shearl1+dr_avg1/2,r_shearl2-dr_avg2/2,Nr_shearl) linspace(r_shearl2+dr_avg2/2,r_farf,Nr-Nr_jet-Nr_shearl)];
        for i=1:Nr_smooth, r_1D = mvg(r_1D); end
    else
        r_1D        = r_1D_equidist;
        dr_farf     = dr;
    end
    
    if (z_core~=0)
        dz_core         = (z_core-z1)/(Nz_core-1);
        dz_rest         = (z2-z_core)/(Nz-Nz_core-1);
        dz_avg          = (dz_core+dz_rest)/2;
        z_1D            = [linspace(z1,z_core-dz_avg/2,Nz_core) linspace(z_core+dz_avg/2,z2,Nz-Nz_core)];
        for i=1:Nz_smooth, z_1D = mvg(z_1D); end
    else
        z_1D    = z_1D_equidist;
        dz_rest = dz;
        dz_core = dz;
    end
else % use M=0.9 grid if Nz=Nr=0
    disp('--> USING M=0.9 LES GRID!');
    load('/home/oschmidt/Codes/matlab/linStab2D_cyl/baseFlows/CylGrid_656_138_longTimeFavreMean.mat','x_mean','r_mean');
    
    % add sponge regions
    dz_inlet        = x_mean(1,2)   - x_mean(1,1);
    dz_outlet       = x_mean(1,end) - x_mean(1,end-1);
    dr_farf         = r_mean(end,1) - r_mean(end-1,1);
    dr_axis         = r_mean(2,1)   - r_mean(1,1);
    r_mean(1,:)     = dr_axis/2;
    
    % add inlet sponge
    z_1D_sponge_inlet   = x_mean(1,1):-dz_inlet:z1;
    z_1D_sponge_inlet   = z_1D_sponge_inlet(2:end);
    z_1D_sponge_inlet   = fliplr(z_1D_sponge_inlet);
    Nz_sponge_inlet     = length(z_1D_sponge_inlet);
    [z_sponge_inlet,r_sponge_inlet]     = meshgrid(z_1D_sponge_inlet,r_mean(:,1));
    x_new = [z_sponge_inlet x_mean];
    r_new = [r_sponge_inlet r_mean];
    
    % add outlet sponge
    z_1D_sponge_outlet  = x_mean(1,end):dz_outlet:z2;
    z_1D_sponge_outlet  = z_1D_sponge_outlet(2:end);
    Nz_sponge_outlet    = length(z_1D_sponge_outlet);
    [z_sponge_outlet,r_sponge_outlet]   = meshgrid(z_1D_sponge_outlet,r_mean(:,1));
    x_new = [x_new z_sponge_outlet];
    r_new = [r_new r_sponge_outlet];
    
    % add far-field sponge
    r_1D_sponge_farf  = r_mean(end,1):dr_farf:r_farf;
    r_1D_sponge_farf  = r_1D_sponge_farf(2:end);
    Nr_sponge_farf    = length(r_1D_sponge_farf);
    [z_sponge_farf,r_sponge_farf]   = meshgrid(x_new(1,:),r_1D_sponge_farf);
    x_new = [x_new; z_sponge_farf];
    r_new = [r_new; r_sponge_farf];
    
    % retrieve 1D grid profiles
    [Nr,Nz]         = size(x_new);
    r_1D            = r_new(:,1);
    z_1D            = x_new(1,:);
    r_farf          = r_1D(end);    % this has to be overwritten
    dr              = r_farf/(Nr-0.5);
    r_1D_equidist   = linspace(dr/2,r_farf,Nr)';
    z_1D_equidist   = linspace(x_new(1),x_new(end),Nz)';
    dz              = z_1D_equidist(2)-z_1D_equidist(1);
    dr_farf         = r_1D(end)-r_1D(end-1);
    dz_rest         = z_1D(end)-z_1D(end-1);
    dz_core         = z_1D(2)-z_1D(1);
    
    NrNz            = Nr*Nz;
end

% apply exponential stretching in damping region
Nsponge_top      = sum(r_1D>=r_sponge);
Nsponge_left     = sum(z_1D<=z_sponge1);
Nsponge_right    = sum(z_1D>=z_sponge2);
% far-field
if aSponge_top>0
    stretch_fact  = exp(aSponge_top*linspace(0,1,Nsponge_top+1)'); stretch_fact  = stretch_fact(2:end);
    count         = 1;
    for i=Nr-Nsponge_top+1:Nr
        r_1D(i)   = r_1D(i-1)+stretch_fact(count)*dr_farf;
        count     = count+1;
    end
end
% outlet
if aSponge_right>0
    stretch_fact  = exp(aSponge_right*linspace(0,1,Nsponge_right+1)'); stretch_fact  = stretch_fact(2:end);
    count         = 1;
    for i=Nz-Nsponge_right+1:Nz
        z_1D(i)   = z_1D(i-1)+stretch_fact(count)*dz_rest;
        count     = count+1;
    end
end
% inlet
if aSponge_left>0
    stretch_fact  = exp(aSponge_left*linspace(0,1,Nsponge_left+1)'); stretch_fact  = stretch_fact(2:end);
    count         = 1;
    for i=Nsponge_left:-1:1
        z_1D(i)   = z_1D(i+1)-stretch_fact(count)*dz_core;
        count     = count+1;
    end
end

% circumvent problem with derivatives on axis
if strcmpi(calcMethod,'calculate_forcings')
    r_1D            = r_1D - 2*(r_1D(2)-r_1D(1));
    r_1D_equidist   = r_1D_equidist - 2*(r_1D(2)-r_1D(1));
end

% radial derivative
[Dr_1D,D2r_1D,Dr_1D_symm,Dr_1D_asymm,D2r_1D_symm,D2r_1D_asymm] = Dmats_SBP(Nr,dr,FDorder);
[Dr_1D,D2r_1D]              = gridtrans1D(r_1D_equidist,Dr_1D,D2r_1D,r_1D);
[Dr_1D_symm,D2r_1D_symm]    = gridtrans1D(r_1D_equidist,Dr_1D_symm,D2r_1D_symm,r_1D);
[Dr_1D_asymm,D2r_1D_asymm]  = gridtrans1D(r_1D_equidist,Dr_1D_asymm,D2r_1D_asymm,r_1D);
% axial derivative
[Dz_1D,D2z_1D,~,~,~,~]      = Dmats_SBP(Nz,dz,FDorder); % fixed to 4th order for now
[Dz_1D,D2z_1D]              = gridtrans1D(z_1D_equidist,Dz_1D,D2z_1D,z_1D);
if pipeBC
    %%%%%%%%%%%%%%%%
    % INCLUDE PIPE %
    %%%%%%%%%%%%%%%%
    % pipe index in r
    r_idx           = findnearest(0.5,r_1D);
    % pipe index in z
    z_idx           = findnearest(0,z_1D);
    Nz_pipe         = z_idx;
    Nz_jet          = Nz-z_idx;
    
    r_1D_pipe_bottom    = r_1D(1:r_idx);
    r_1D_pipe_top       = r_1D(r_idx+1:end);
    
    % region inside pipe
    Nr_pipe_bottom      = r_idx;
    dr_pipe_bottom      = r_1D_pipe_bottom(end)/(Nr_pipe_bottom-0.5);
    [Dr_1D_pipe_bottom,D2r_1D_pipe_bottom,Dr_1D_pipe_bottom_symm,Dr_1D_pipe_bottom_asymm,D2r_1D_pipe_bottom_symm,D2r_1D_pipe_bottom_asymm] = Dmats_SBP(Nr_pipe_bottom,dr_pipe_bottom,FDorder);
    r_1D_equidist_pipe_bottom               = linspace(r_1D_pipe_bottom(1),r_1D_pipe_bottom(end),Nr_pipe_bottom);
    r_1D_pipe_bottom                        = r_1D(1:r_idx);
    [Dr_1D_pipe_bottom,D2r_1D_pipe_bottom]  = gridtrans1D(r_1D_equidist_pipe_bottom,Dr_1D_pipe_bottom,D2r_1D_pipe_bottom,r_1D_pipe_bottom);
    % Mohseni's pole treatment
    [Dr_1D_pipe_bottom_symm,D2r_1D_pipe_bottom_symm]    = gridtrans1D(r_1D_equidist,Dr_1D_pipe_bottom_symm,D2r_1D_pipe_bottom_symm,r_1D_pipe_bottom);
    [Dr_1D_pipe_bottom_asymm,D2r_1D_pipe_bottom_asymm]  = gridtrans1D(r_1D_equidist,Dr_1D_pipe_bottom_asymm,D2r_1D_pipe_bottom_asymm,r_1D_pipe_bottom);
    
    % region above pipe
    Nr_pipe_top         = Nr-r_idx;
    dr_pipe_top         = (r_1D_pipe_top(end)-r_1D_pipe_top(1))/Nr_pipe_top;
    [Dr_1D_pipe_top,D2r_1D_pipe_top,Dr_1D_pipe_top_symm,Dr_1D_pipe_top_asymm,D2r_1D_pipe_top_symm,D2r_1D_pipe_top_asymm] = Dmats_SBP(Nr_pipe_top,dr_pipe_top,FDorder);
    r_1D_equidist_pipe_top              = linspace(r_1D_pipe_top(1),r_1D_pipe_top(end),Nr_pipe_top);
    [Dr_1D_pipe_top,D2r_1D_pipe_top]    = gridtrans1D(r_1D_equidist_pipe_top,Dr_1D_pipe_top,D2r_1D_pipe_top,r_1D_pipe_top);
    % Mohseni's pole treatment
    [Dr_1D_pipe_top_symm,D2r_1D_pipe_top_symm]      = gridtrans1D(r_1D_equidist,Dr_1D_pipe_top_symm,D2r_1D_pipe_top_symm,r_1D_pipe_top);
    [Dr_1D_pipe_top_asymm,D2r_1D_pipe_top_asymm]    = gridtrans1D(r_1D_equidist,Dr_1D_pipe_top_asymm,D2r_1D_pipe_top_asymm,r_1D_pipe_top);
    
    % build derivative matrix (including row of zeros for pipe index)
    Dr_1D_pipe          = blkdiag(Dr_1D_pipe_bottom,Dr_1D_pipe_top);
    D2r_1D_pipe         = blkdiag(D2r_1D_pipe_bottom,D2r_1D_pipe_top);
    % Mohseni's pole treatment
    Dr_1D_pipe_symm          = blkdiag(Dr_1D_pipe_bottom_symm,Dr_1D_pipe_top_symm);
    D2r_1D_pipe_symm         = blkdiag(D2r_1D_pipe_bottom_symm,D2r_1D_pipe_top_symm);
    Dr_1D_pipe_asymm         = blkdiag(Dr_1D_pipe_bottom_asymm,Dr_1D_pipe_top_asymm);
    D2r_1D_pipe_asymm        = blkdiag(D2r_1D_pipe_bottom_asymm,D2r_1D_pipe_top_asymm);
    
    % combine full differentiation matrix in r
    Dr      = blkdiag(kron(speye(Nz_pipe),Dr_1D_pipe),  kron(speye(Nz_jet),Dr_1D));
    D2r     = blkdiag(kron(speye(Nz_pipe),D2r_1D_pipe), kron(speye(Nz_jet),D2r_1D));
    
    % Mohseni's pole treatment
    Dr_symm      = blkdiag(kron(speye(Nz_pipe),Dr_1D_pipe_symm),  kron(speye(Nz_jet),Dr_1D_symm));
    D2r_symm     = blkdiag(kron(speye(Nz_pipe),D2r_1D_pipe_symm), kron(speye(Nz_jet),D2r_1D_symm));
    Dr_asymm     = blkdiag(kron(speye(Nz_pipe),Dr_1D_pipe_asymm), kron(speye(Nz_jet),Dr_1D_asymm));
    D2r_asymm    = blkdiag(kron(speye(Nz_pipe),D2r_1D_pipe_asymm),kron(speye(Nz_jet),D2r_1D_asymm));
    
    % axial derivative
    z_1D_pipe      = z_1D(1:z_idx-1);
    z_1D_jet       = z_1D(z_idx:end);
    
    % region right of pipe
    dz_jet         = (z_1D_jet(end)-z_1D_jet(1))/(Nz_jet+1);
    [Dz_1D_jet,D2z_1D_jet,~,~,~,~] = Dmats_SBP(Nz_jet+1,dz_jet,FDorder);
    z_1D_equidist_jet       = linspace(z_1D_jet(1),z_1D_jet(end),Nz_jet+1);
    [Dz_1D_jet,D2z_1D_jet]  = gridtrans1D(z_1D_equidist_jet,Dz_1D_jet,D2z_1D_jet,z_1D_jet);
    
    % region left of pipe
    dz_pipe         = (z_1D_pipe(end)-z_1D_pipe(1))/(Nz_pipe-1);
    [Dz_1D_pipe,D2z_1D_pipe,~,~,~,~] = Dmats_SBP(Nz_pipe-1,dz_pipe,FDorder);
    z_1D_equidist_pipe          = linspace(z_1D_pipe(1),z_1D_pipe(end),Nz_pipe-1);
    [Dz_1D_pipe,D2z_1D_pipe]    = gridtrans1D(z_1D_equidist_pipe,Dz_1D_pipe,D2z_1D_pipe,z_1D_pipe);
    
    % build derivative matrix
    Dz_1D_jet           = blkdiag(Dz_1D_pipe,Dz_1D_jet);
    D2z_1D_jet          = blkdiag(D2z_1D_pipe,D2z_1D_jet);
    
    % combine full differentiation matrix in r
    Dz                          = kron(Dz_1D,speye(Nr));
    D2z                         = kron(D2z_1D,speye(Nr));
    % and replace z derivative along pipe coordinate
    pipe_row_idx_bottom         = (1:Nr:Nr*Nz) + (r_idx-1);
    pipe_row_idx_top            = (1:Nr:Nr*Nz) + (r_idx);
    Dz(pipe_row_idx_bottom,pipe_row_idx_bottom)     = Dz_1D_jet;
    D2z(pipe_row_idx_bottom,pipe_row_idx_bottom)    = D2z_1D_jet;
    Dz(pipe_row_idx_top,pipe_row_idx_top)           = Dz_1D_jet;
    D2z(pipe_row_idx_top,pipe_row_idx_top)          = D2z_1D_jet;
    
    % indices for BCs
    pipe_idx_bottom     = pipe_row_idx_bottom(1:z_idx);
    pipe_idx_top        = pipe_row_idx_top(1:z_idx);
else
    Dz          = sparse(kron(Dz_1D,speye((Nr))));
    Dr          = sparse(kron(speye((Nz)),Dr_1D));
    Dr_symm     = sparse(kron(speye((Nz)),Dr_1D_symm));
    Dr_asymm    = sparse(kron(speye((Nz)),Dr_1D_asymm));
    D2z         = sparse(kron(D2z_1D,speye(Nr)));
    D2r         = sparse(kron(speye(Nz),D2r_1D));
    D2r_symm    = sparse(kron(speye((Nz)),D2r_1D_symm));
    D2r_asymm   = sparse(kron(speye((Nz)),D2r_1D_asymm));
    
    pipe_idx_bottom     = [];
    pipe_idx_top        = [];
end

% build mesh
[z, r]  = meshgrid(z_1D, r_1D);

Z = 0*speye(NrNz,NrNz); % need a zero diagonal sparse matrix... is there a better way?
if strcmpi(calcMethod,'calculate_forcings')
    DR  = sparse([ ...
        Dr        Z         Z           Z        Z
        Z         Dr        Z           Z        Z
        Z         Z         Dr          Z        Z
        Z         Z         Z           Dr       Z
        Z         Z         Z           Z        Dr
        ]);
elseif m==0
    DR  = sparse([ ...
        Dr_symm   Z         Z           Z        Z
        Z         Dr_asymm  Z           Z        Z
        Z         Z         Dr_asymm    Z        Z
        Z         Z         Z           Dr_symm  Z
        Z         Z         Z           Z        Dr_symm
        ]);
elseif m==1
    DR  = sparse([ ...
        Dr_asymm  Z         Z           Z        Z
        Z         Dr_symm   Z           Z        Z
        Z         Z         Dr_symm     Z        Z
        Z         Z         Z           Dr_asymm Z
        Z         Z         Z           Z        Dr_asymm
        ]);
elseif m>=2
    DR  = sparse([ ...
        Dr_asymm  Z         Z           Z        Z
        Z         Dr_asymm  Z           Z        Z
        Z         Z         Dr_asymm    Z        Z
        Z         Z         Z           Dr_asymm Z
        Z         Z         Z           Z        Dr_asymm
        ]);
end
DZ      = kron(speye(5,5),Dz);
if strcmpi(calcMethod,'calculate_forcings')
    D2R = sparse([ ...
        D2r       Z         Z           Z        Z
        Z         D2r       Z           Z        Z
        Z         Z         D2r         Z        Z
        Z         Z         Z           D2r      Z
        Z         Z         Z           Z        D2r
        ]);
elseif m==0
    D2R = sparse([ ...
        D2r_symm  Z         Z           Z        Z
        Z         D2r_asymm Z           Z        Z
        Z         Z         D2r_asymm   Z        Z
        Z         Z         Z           D2r_symm Z
        Z         Z         Z           Z        D2r_symm
        ]);
elseif m==1
    D2R = sparse([ ...
        D2r_asymm Z         Z           Z         Z
        Z         D2r_symm  Z           Z         Z
        Z         Z         D2r_symm    Z         Z
        Z         Z         Z           D2r_asymm Z
        Z         Z         Z           Z         D2r_asymm
        ]);
elseif m>=2
    D2R = sparse([ ...
        D2r_asymm Z         Z           Z         Z
        Z         D2r_asymm Z           Z         Z
        Z         Z         D2r_asymm   Z         Z
        Z         Z         Z           D2r_asymm Z
        Z         Z         Z           Z         D2r_asymm
        ]);
end
clear Z
D2Z     = kron(speye(5,5),D2z);

% BC derivative matricies (WRONG!?)
BSr_1D                  = zeros(Nr,Nr);
BSr_1D(1,1:4)           = -[11/6 -3 3/2 -1/3];
BSr_1D(end,  end-3:end) = fliplr([11/6 -3 3/2 -1/3]);
BSr_1D                  = 1/dr* BSr_1D;
BSz_1D          = zeros(Nz,Nz);
BSz_1D(1,1:4)   = -[11/6 -3 3/2 -1/3];
BSz_1D(end,  end-3:end) = fliplr([11/6 -3 3/2 -1/3]);
BSz_1D          = 1/dz * BSz_1D;

if sum((r_1D(:)~=r_1D_equidist(:)))~=0
    [BSr_1D,~]                  = gridtrans1D(r_1D_equidist,BSr_1D,BSr_1D,r_1D);
end
if sum(z_1D(:)~=z_1D_equidist(:))~=0
    [BSz_1D,~]                  = gridtrans1D(z_1D_equidist,BSz_1D,BSz_1D,z_1D);
end

BSz     = sparse(kron(BSz_1D,speye((Nr))));
BSr     = sparse(kron(speye((Nz)),BSr_1D));
BSR     = kron(speye(5,5),BSr);
BSZ     = kron(speye(5,5),BSz);
time = toc;
disp(['    elapsed time - Differentiation matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);

tic
%% jet data
if strcmp(baseState,'69M')
    [RHO,U,V,W,P,T] = get_LES_meanFlow_69M(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet);
elseif strcmp(baseState,'M04')
%     [RHO,U,V,W,P,T,REY_XR,TKE] = get_LES_meanFlow_M04(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet,pipe_idx_top,pipe_idx_bottom);
    [RHO,U,V,W,P,T,TKE] = get_LES_meanFlow_M04(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet,pipe_idx_top,pipe_idx_bottom);
elseif strcmp(baseState,'B118')
    [RHO,U,V,W,P,T] = get_LES_meanFlow_B118(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet,pipe_idx_top,pipe_idx_bottom);
elseif strcmp(baseState,'M09')
    [RHO,U,V,W,P,T,REY_XR] = get_LES_meanFlow(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet,pipe_idx_top,pipe_idx_bottom);
elseif strcmp(baseState,'M09laminar')
    [RHO,U,V,W,P,T,REY_XR] = get_LES_meanFlow_M09laminar(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet,pipe_idx_top,pipe_idx_bottom);
elseif strcmp(baseState,'Cambridge')
    [RHO,U,V,W,P,T] = get_LES_meanFlow_Cambridge(z,r,Nsmooth_z,Nsmooth_r,x_LES_inlet);
else
    error('Unknown base-state!');
end
time = toc;
disp(['    elapsed time - Base-flow interpolation: ' datestr(time/24/3600, 'HH:MM:SS')]);

if (vizBaseAndRes)&&(calcCount==1)
    fig_grid = figure('name','Grid');
    pcolor_h=pcolor(z, r, W); set(pcolor_h,'EdgeColor','none'), hold on
    mesh(z,r,zeros(size(z)),'FaceColor','none','edgecolor','k'), text(z(1), r(end), ' $$w$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
    xlim([z(1) z(1)+1]), ylim([0.25 0.75]), daspect([1 1 1])
    
    fig_base = figure('name','Base-state');
    subplot(3,3,1)
    pcolor_h=pcolor(z, r, RHO); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\rho$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,2)
    pcolor_h=pcolor(z, r, U); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$u$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,3)
    pcolor_h=pcolor(z, r, W); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$w$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    contour(z, r, W,[0.99 0.99],'r-');
    subplot(3,3,4)
    pcolor_h=pcolor(z, r, P); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$p$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    drawnow
end

% Dependant coefficients
cv      = 1/(kappa*(kappa-1)*Ma^2);
c1      = (kappa-1)*Re*Pr*Ma^2;
c2      = kappa*Ma^2;

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

% Reshaping flowfield
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
dUdr    = Dr*U;     d2Udr2    = D2r*U;
dVdr    = Dr*V;     d2Vdr2    = D2r*V;
dWdr    = Dr*W;     d2Wdr2    = D2r*W;
dRHOdr  = Dr*RHO;   d2RHOdr2  = D2r*RHO;
dTdr    = Dr*T;     d2Tdr2    = D2r*T;
dUdz    = Dz*U;     d2Udz2    = D2z*U;
dVdz    = Dz*V;     d2Vdz2    = D2z*V;
dWdz    = Dz*W;     d2Wdz2    = D2z*W;
dRHOdz  = Dz*RHO;   d2RHOdz2  = D2z*RHO;
dTdz    = Dz*T;     d2Tdz2    = D2z*T;
d2Udrz  = Dz*dUdr;
d2Vdrz  = Dz*dVdr;
d2Wdrz  = Dz*dWdr;

% set up sponge region
spongeFun       = zeros(Nr,Nz);
for i = 1:Nr
    for j = 1:Nz
        % left
        if j<=Nsponge_left
            loc     = 1-(j-1)/(Nsponge_left);
            if strcmp(spongeType,'poly5')
                poly    = -(-6*loc.^5+15*loc.^4-10*loc^3);
            elseif strcmp(spongeType,'pow')
                poly    = loc.^sponge_pow;
            end
            if spongeFun(i,j)<spongeEps_left*poly; spongeFun(i,j) = spongeEps_left*poly; end
        end
        % right
        if j>Nz-Nsponge_right
            loc     = 1+(j-Nz)/Nsponge_right;
            if strcmp(spongeType,'poly5')
                poly    = -(-6*loc.^5+15*loc.^4-10*loc^3);
            elseif strcmp(spongeType,'pow')
                poly    = loc.^sponge_pow;
            end
            if spongeFun(i,j)<spongeEps_right*poly; spongeFun(i,j) = spongeEps_right*poly; end
        end
        % top
        if i>Nr-Nsponge_top
            loc     = 1+(i-Nr)/Nsponge_top;
            if strcmp(spongeType,'poly5')
                poly    = -(-6*loc.^5+15*loc.^4-10*loc^3);
            elseif strcmp(spongeType,'pow')
                poly    = loc.^sponge_pow;
            end
            if spongeFun(i,j)<spongeEps_top*poly; spongeFun(i,j) = spongeEps_top*poly; end
        end
    end
end
if exist('spongePipeFactor','var'), spongeFun(z<=0&r<=0.5) = spongeFun(z<=0&r<=0.5)*(spongePipeFactor); end

% if pipeBC, spongeFun(r>=r_1D_pipe_top(1)&z<=z(pipe_idx_top(end))) = spongeFun(pipe_idx_top(1)); end
% spongeFun(r>=0.5&z<=0) = spongeEps_left;

spongeFun   = spongeFun(:);
Asponge     = kron(speye(5),diag(sparse(spongeFun)));

% calculate, mask and smooth eddy viscosity
if useEddyVisc
    if isa(eddyVisc,'numeric')
        muT         = eddyVisc*ones(Nr*Nz,1);
        S_zr        = zeros(Nr*Nz,1);
        REY_XR      = zeros(Nr,Nz);
    else
        S_zr        = 1/2*(dUdz + dWdr);
        muT_mask    = REY_XR>=5e-4;
        muT         = -REY_XR(:)./(2*S_zr).*muT_mask(:);
        muT         = reshape(muT,Nr,Nz);
        muT         = triSmooth(muT,  1, Nsmooth_nuT);
        muT         = triSmooth(muT,  2, Nsmooth_nuT);
        muT         = muT(:);
    end
end
time    = toc;
disp(['    elapsed time - Base flow derivatives: ' datestr(time/24/3600, 'HH:MM:SS')]);


%%
if (vizBaseAndRes)&&(calcCount==1)
    figure(fig_base)
    subplot(3,3,5)
    pcolor_h=pcolor(z, r, reshape(spongeFun,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' sponge', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,6)
    pcolor_h=pcolor(z, r, reshape(W,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$W$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,7)
    pcolor_h=pcolor(z, r, reshape(dWdr,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$w$$/d$$r$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,8)
    pcolor_h=pcolor(z, r, reshape(dWdz,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$w$$/d$$z$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,9)
    pcolor_h=pcolor(z, r, reshape(d2Wdrz,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$^2w$$/d$$r$$d$$z$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    set(gcf,'Position',[171         103        1324         666])
    
    if useEddyVisc
        figure('name','Reynolds stress')
        subplot(3,1,1)
        pcolor_h=pcolor(z, r, REY_XR); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\overline{\rho u_x'' u_r''}$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(3,1,2)
        pcolor_h=pcolor(z, r, reshape(S_zr,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\overline{S_{xr}}=1/2(du_r/dz + du_z/dr)$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(3,1,3)
        pcolor_h=pcolor(z, r, reshape(muT/mu_0,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\mu_T$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    end
    
    drawnow
end
%%
if strcmpi(calcMethod,'les_forcingPOD_response')
    POD_file            = [outsubdir '/POD_m' num2str(m,'%.2i') '_freq' num2str(St_save_i(St_cmp_i),'%.2i') POD_comment '.mat'];
    load(POD_file,'psi_T_hat_POD','psi_rho_hat_POD','psi_u_x_hat_POD','psi_u_r_hat_POD','psi_u_th_hat_POD','n_save_POD');
    V_in                = zeros(5*NrNz,n_save_POD);
    for lambda_i = 1:n_save_POD
        V_in(0*NrNz+1:1*NrNz,lambda_i) = reshape(interp2(x_les,r_les,squeeze(psi_rho_hat_POD(:,:,lambda_i)),z,r,'spline',0),NrNz,1);
        V_in(1*NrNz+1:2*NrNz,lambda_i) = reshape(interp2(x_les,r_les,squeeze(psi_u_r_hat_POD(:,:,lambda_i)),z,r,'spline',0),NrNz,1);
        V_in(2*NrNz+1:3*NrNz,lambda_i) = reshape(interp2(x_les,r_les,squeeze(psi_u_th_hat_POD(:,:,lambda_i)),z,r,'spline',0),NrNz,1);
        V_in(3*NrNz+1:4*NrNz,lambda_i) = reshape(interp2(x_les,r_les,squeeze(psi_u_x_hat_POD(:,:,lambda_i)),z,r,'spline',0),NrNz,1);
        V_in(4*NrNz+1:5*NrNz,lambda_i) = reshape(interp2(x_les,r_les,squeeze(psi_T_hat_POD(:,:,lambda_i)),z,r,'spline',0),NrNz,1);
    end
    
    if (vizBaseAndRes)&&(calcCount==1)
        %%
        figure('name','forcing POD modes','position',[100 500 1600 400])
        subplot(2,3,1)
        pcolor_h=pcolor(z, r, real(reshape(V_in(0*NrNz+1:1*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{\rho}''}^{(' num2str(1) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,2)
        pcolor_h=pcolor(z, r, real(reshape(V_in(1*NrNz+1:2*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{r}''}^{(' num2str(1) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,3)
        pcolor_h=pcolor(z, r, real(reshape(V_in(2*NrNz+1:3*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{\theta}''}^{(' num2str(1) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,4)
        pcolor_h=pcolor(z, r, real(reshape(V_in(3*NrNz+1:4*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{x}''}^{(' num2str(1) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,5)
        pcolor_h=pcolor(z, r, real(reshape(V_in(4*NrNz+1:5*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{T}''}^{(' num2str(1) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        drawnow
    end
end

%%
if strcmpi(calcMethod,'les_forcingFFT_response')
    V_in                = zeros(5*NrNz,nblks);
    for blk_i = 1:nblks
        disp(['    --> loading FFT block ' num2str(blk_i) '/' num2str(nblks)]);
        FFT_file    = [outsubdir '/m' num2str(m,'%.2i') '_block' num2str(blk_i,'%.4i') '.mat'];
        load(FFT_file, 'p_hat','rho_hat','u_x_hat','u_r_hat','u_th_hat');
        
        rho_hat                     = reshape(interp2(x_les,r_les,squeeze(rho_hat(St_cmp_i,:,:)),z,r,'spline',0),NrNz,1);
        V_in(0*NrNz+1:1*NrNz,blk_i) = rho_hat;
        V_in(1*NrNz+1:2*NrNz,blk_i) = reshape(interp2(x_les,r_les,squeeze(u_r_hat(St_cmp_i,:,:)),z,r,'spline',0),NrNz,1);
        V_in(2*NrNz+1:3*NrNz,blk_i) = reshape(interp2(x_les,r_les,squeeze(u_th_hat(St_cmp_i,:,:)),z,r,'spline',0),NrNz,1);
        V_in(3*NrNz+1:4*NrNz,blk_i) = reshape(interp2(x_les,r_les,squeeze(u_x_hat(St_cmp_i,:,:)),z,r,'spline',0),NrNz,1);
        p_hat                       = reshape(interp2(x_les,r_les,squeeze(p_hat(St_cmp_i,:,:)),z,r,'spline',0),NrNz,1);
        V_in(4*NrNz+1:5*NrNz,blk_i) = ((P(:)+p_hat)*kappa*Ma^2 - RHO.*T - rho_hat.*T)./(RHO+rho_hat); % get perturbation temperature from equation of state
    end
    clear p_hat rho_hat u_x_hat u_r_hat u_th_hat
    
    if (vizBaseAndRes)&&(calcCount==1)
        %%
        figure('name','forcing FFT modes','position',[100 500 1600 400])
        subplot(2,3,1)
        pcolor_h=pcolor(z, r, real(reshape(V_in(0*NrNz+1:1*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{\rho}''}^{(' num2str(blk_i) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,2)
        pcolor_h=pcolor(z, r, real(reshape(V_in(1*NrNz+1:2*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{r}''}^{(' num2str(blk_i) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,3)
        pcolor_h=pcolor(z, r, real(reshape(V_in(2*NrNz+1:3*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{\theta}''}^{(' num2str(blk_i) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,4)
        pcolor_h=pcolor(z, r, real(reshape(V_in(3*NrNz+1:4*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{u}_{x}''}^{(' num2str(blk_i) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        subplot(2,3,5)
        pcolor_h=pcolor(z, r, real(reshape(V_in(4*NrNz+1:5*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['$\Psi_{\hat{T}''}^{(' num2str(blk_i) ')} (m='  num2str(m) ', St=' num2str(St_save(St_cmp_i),'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
        drawnow
    end
end
%%
tic
Z = zeros(NrNz, 1);
I = ones(NrNz, 1);
R = reshape(r, NrNz, 1);

% Coefficient matrices
if useEddyVisc
    % add turbulent viscosity
    MUT = muT/mu_0;
    getCoeffsTurbEddyVisc
else
    getCoeffsLaminar
end


A0 =    ...
    [   ...
    diag(sparse(A0_11)) diag(sparse(A0_12)) diag(sparse(A0_13)) diag(sparse(A0_14)) diag(sparse(A0_15))
    diag(sparse(A0_21)) diag(sparse(A0_22)) diag(sparse(A0_23)) diag(sparse(A0_24)) diag(sparse(A0_25))
    diag(sparse(A0_31)) diag(sparse(A0_32)) diag(sparse(A0_33)) diag(sparse(A0_34)) diag(sparse(A0_35))
    diag(sparse(A0_41)) diag(sparse(A0_42)) diag(sparse(A0_43)) diag(sparse(A0_44)) diag(sparse(A0_45))
    diag(sparse(A0_51)) diag(sparse(A0_52)) diag(sparse(A0_53)) diag(sparse(A0_54)) diag(sparse(A0_55))
    ];

Ar =    ...
    [   ...
    diag(sparse(Ar_11)) diag(sparse(Ar_12)) diag(sparse(Ar_13)) diag(sparse(Ar_14)) diag(sparse(Ar_15))
    diag(sparse(Ar_21)) diag(sparse(Ar_22)) diag(sparse(Ar_23)) diag(sparse(Ar_24)) diag(sparse(Ar_25))
    diag(sparse(Ar_31)) diag(sparse(Ar_32)) diag(sparse(Ar_33)) diag(sparse(Ar_34)) diag(sparse(Ar_35))
    diag(sparse(Ar_41)) diag(sparse(Ar_42)) diag(sparse(Ar_43)) diag(sparse(Ar_44)) diag(sparse(Ar_45))
    diag(sparse(Ar_51)) diag(sparse(Ar_52)) diag(sparse(Ar_53)) diag(sparse(Ar_54)) diag(sparse(Ar_55))
    ];

Az =    ...
    [   ...
    diag(sparse(Az_11)) diag(sparse(Az_12)) diag(sparse(Az_13)) diag(sparse(Az_14)) diag(sparse(Az_15))
    diag(sparse(Az_21)) diag(sparse(Az_22)) diag(sparse(Az_23)) diag(sparse(Az_24)) diag(sparse(Az_25))
    diag(sparse(Az_31)) diag(sparse(Az_32)) diag(sparse(Az_33)) diag(sparse(Az_34)) diag(sparse(Az_35))
    diag(sparse(Az_41)) diag(sparse(Az_42)) diag(sparse(Az_43)) diag(sparse(Az_44)) diag(sparse(Az_45))
    diag(sparse(Az_51)) diag(sparse(Az_52)) diag(sparse(Az_53)) diag(sparse(Az_54)) diag(sparse(Az_55))
    ];

Arz =   ...
    [   ...
    diag(sparse(Arz_11)) diag(sparse(Arz_12)) diag(sparse(Arz_13)) diag(sparse(Arz_14)) diag(sparse(Arz_15))
    diag(sparse(Arz_21)) diag(sparse(Arz_22)) diag(sparse(Arz_23)) diag(sparse(Arz_24)) diag(sparse(Arz_25))
    diag(sparse(Arz_31)) diag(sparse(Arz_32)) diag(sparse(Arz_33)) diag(sparse(Arz_34)) diag(sparse(Arz_35))
    diag(sparse(Arz_41)) diag(sparse(Arz_42)) diag(sparse(Arz_43)) diag(sparse(Arz_44)) diag(sparse(Arz_45))
    diag(sparse(Arz_51)) diag(sparse(Arz_52)) diag(sparse(Arz_53)) diag(sparse(Arz_54)) diag(sparse(Arz_55))
    ];

Arr =   ...
    [   ...
    diag(sparse(Arr_11)) diag(sparse(Arr_12)) diag(sparse(Arr_13)) diag(sparse(Arr_14)) diag(sparse(Arr_15))
    diag(sparse(Arr_21)) diag(sparse(Arr_22)) diag(sparse(Arr_23)) diag(sparse(Arr_24)) diag(sparse(Arr_25))
    diag(sparse(Arr_31)) diag(sparse(Arr_32)) diag(sparse(Arr_33)) diag(sparse(Arr_34)) diag(sparse(Arr_35))
    diag(sparse(Arr_41)) diag(sparse(Arr_42)) diag(sparse(Arr_43)) diag(sparse(Arr_44)) diag(sparse(Arr_45))
    diag(sparse(Arr_51)) diag(sparse(Arr_52)) diag(sparse(Arr_53)) diag(sparse(Arr_54)) diag(sparse(Arr_55))
    ];

Azz =   ...
    [   ...
    diag(sparse(Azz_11)) diag(sparse(Azz_12)) diag(sparse(Azz_13)) diag(sparse(Azz_14)) diag(sparse(Azz_15))
    diag(sparse(Azz_21)) diag(sparse(Azz_22)) diag(sparse(Azz_23)) diag(sparse(Azz_24)) diag(sparse(Azz_25))
    diag(sparse(Azz_31)) diag(sparse(Azz_32)) diag(sparse(Azz_33)) diag(sparse(Azz_34)) diag(sparse(Azz_35))
    diag(sparse(Azz_41)) diag(sparse(Azz_42)) diag(sparse(Azz_43)) diag(sparse(Azz_44)) diag(sparse(Azz_45))
    diag(sparse(Azz_51)) diag(sparse(Azz_52)) diag(sparse(Azz_53)) diag(sparse(Azz_54)) diag(sparse(Azz_55))
    ];

ZZ = 0*speye(NrNz,NrNz); % need a zero diagonal sparse matrix... is there a better way?

RHS =     ...
    [   ...
    diag(sparse(B_11)) diag(sparse(B_12)) diag(sparse(B_13)) diag(sparse(B_14)) diag(sparse(B_15))
    diag(sparse(B_21)) diag(sparse(B_22)) diag(sparse(B_23)) diag(sparse(B_24)) diag(sparse(B_25))
    diag(sparse(B_31)) diag(sparse(B_32)) diag(sparse(B_33)) diag(sparse(B_34)) diag(sparse(B_35))
    diag(sparse(B_41)) diag(sparse(B_42)) diag(sparse(B_43)) diag(sparse(B_44)) diag(sparse(B_45))
    diag(sparse(B_51)) diag(sparse(B_52)) diag(sparse(B_53)) diag(sparse(B_54)) diag(sparse(B_55))
    ];

% clean up memory
clear A0_* Ar_* Az_* Arr_* Azz_* Arz_* B_* d2* dU* dV* dW* dRHO* dT* dmu* I Z ZZ
LHS     = A0 + Ar*DR + Az*DZ + Arr*D2R + Azz*D2Z + Arz*DR*DZ - Asponge;
clear A0 Ar Az Arr Azz Arz D2R D2Z Asponge

time    = toc;
disp(['    elapsed time - Coefficient matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);

tic
%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%
% Row indices for primitive variables
li = 1 : Nr;
ri = li + NrNz - Nr;
bi = 1 : Nr : NrNz;
ti = bi + Nr - 1;

% inlet (pipe section)
li_pipe = li(r_1D<=0.5);
% inlet (farfield section)
li_farf = li(r_1D>0.5);

ti_rho  = ti;
bi_rho  = bi;
ri_rho  = ri;
li_rho  = li;
li_rho_pipe  = li_pipe;
li_rho_farf  = li_farf;

ti_u    = ti   + NrNz;
bi_u    = bi   + NrNz;
ri_u    = ri   + NrNz;
li_u    = li   + NrNz;
li_u_pipe  = li_pipe + NrNz;
li_u_farf  = li_farf + NrNz;

ti_v    = ti   + 2*NrNz;
bi_v    = bi   + 2*NrNz;
ri_v    = ri   + 2*NrNz;
li_v    = li   + 2*NrNz;
li_v_pipe  = li_pipe + 2*NrNz;
li_v_farf  = li_farf + 2*NrNz;

ti_w    = ti   + 3*NrNz;
bi_w    = bi   + 3*NrNz;
ri_w    = ri   + 3*NrNz;
li_w    = li   + 3*NrNz;
li_w_pipe  = li_pipe + 3*NrNz;
li_w_farf  = li_farf + 3*NrNz;

ti_T    = ti   + 4*NrNz;
bi_T    = bi   + 4*NrNz;
ri_T    = ri   + 4*NrNz;
li_T    = li   + 4*NrNz;
li_T_pipe  = li_pipe + 4*NrNz;
li_T_farf  = li_farf + 4*NrNz;

% Column indices for conserved variables
rho_j   =        1 :   NrNz;
u_j     =   NrNz+1 : 2*NrNz;
v_j     = 2*NrNz+1 : 3*NrNz;
w_j     = 3*NrNz+1 : 4*NrNz;
T_j     = 4*NrNz+1 : 5*NrNz;

% farfield
rho_farf_i  = ti_rho;
u_farf_i  = ti_u;
v_farf_i  = ti_v;
w_farf_i  = ti_w;
T_farf_i  = ti_T;
% axis
rho_axis_i  = bi_rho;
u_axis_i  = bi_u;
v_axis_i  = bi_v;
w_axis_i  = bi_w;
T_axis_i  = bi_T;
% inlet
rho_in_i    = li_rho;
u_in_i    = li_u;
v_in_i    = li_v;
w_in_i    = li_w;
T_in_i    = li_T;
% inlet (pipe section)
rho_pipe_in_i  = li_rho_pipe;
u_pipe_in_i    = li_u_pipe;
v_pipe_in_i    = li_v_pipe;
w_pipe_in_i    = li_w_pipe;
T_pipe_in_i    = li_T_pipe;
% inlet (farfield section)
rho_farf_in_i  = li_rho_farf;
u_farf_in_i    = li_u_farf;
v_farf_in_i    = li_v_farf;
w_farf_in_i    = li_w_farf;
T_farf_in_i    = li_T_farf;
% outlet
rho_out_i = ri_rho;
u_out_i   = ri_u;
v_out_i   = ri_v;
w_out_i   = ri_w;
T_out_i   = ri_T;

%%
% nozzle
rho_pipe_i    = [pipe_idx_bottom pipe_idx_top];
u_pipe_i      = rho_pipe_i + 1*NrNz;
v_pipe_i      = rho_pipe_i + 2*NrNz;
w_pipe_i      = rho_pipe_i + 3*NrNz;
T_pipe_i      = rho_pipe_i + 4*NrNz;

% figure('name','Grid')
% pcolor_h = pcolor(z, r, reshape(W,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal tight, hold on
% plot(z(pipe_idx_bottom), r(pipe_idx_bottom),'rx')
% plot(z(pipe_idx_top), r(pipe_idx_top),'bx')
% drawnow
%%

size_LHS = [5*NrNz 5*NrNz];

% use 3rd order boudary derivative for SBP scheme
DR = BSR;
DZ = BSZ;
clear BS*

if ~strcmpi(calcMethod,'calculate_forcings')
    % far-field
    if strcmpi(farfieldBC,'zero_gradient')
        index_set = [u_farf_i v_farf_i w_farf_i rho_farf_i T_farf_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        
        [RHS_i,RHS_j,RHS_val] = find(RHS);
        DR_i    = zeros(4*length(index_set),1); % 5-point stencil with 4 entries
        DR_j    = DR_i;
        DR_val  = DR_i;
        for row_i=1:length(index_set);
            DR_i((row_i-1)*4+1:row_i*4)    = [index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
            [~,jj,val] = find(DR(index_set(row_i),:));
            DR_j((row_i-1)*4+1:row_i*4)    = jj;
            DR_val((row_i-1)*4+1:row_i*4)  = val;
        end
        % append derivative matrix elements
        RHS_i       = [RHS_i; DR_i];
        RHS_j       = [RHS_j; DR_j];
        RHS_val     = [RHS_val; DR_val];
        RHS = sparse(RHS_i,RHS_j,RHS_val);
        %     RHS(  u_farf_i,   u_j)   = DR(  u_farf_i,   u_j);
        %     RHS(  v_farf_i,   v_j)   = DR(  v_farf_i,   v_j);
        %     RHS(  w_farf_i,   w_j)   = DR(  w_farf_i,   w_j);
        %     RHS(rho_farf_i, rho_j)   = DR(rho_farf_i, rho_j);
        %     RHS(  T_farf_i,   T_j)   = DR(  T_farf_i,   T_j);
    else
        index_set = [u_farf_i v_farf_i w_farf_i T_farf_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        RHS = RHS + Diag;
    end
    
    % inlet
    if strcmpi(inletBC,'zero_gradient_pipe')
        % inlet (pipe section): zero gradient: all but rho
        index_set = [u_pipe_in_i v_pipe_in_i w_pipe_in_i rho_pipe_in_i T_pipe_in_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        RHS(  u_pipe_in_i,   u_j)   = DZ(  u_pipe_in_i,   u_j);
        RHS(  v_pipe_in_i,   v_j)   = DZ(  v_pipe_in_i,   v_j);
        RHS(  w_pipe_in_i,   w_j)   = DZ(  w_pipe_in_i,   w_j);
        RHS(rho_pipe_in_i, rho_j)   = DZ(rho_pipe_in_i, rho_j);
        RHS(  T_pipe_in_i,   T_j)   = DZ(  T_pipe_in_i,   T_j);
        
        % inlet (farfield section)
        index_set = [u_farf_in_i v_farf_in_i w_farf_in_i T_farf_in_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        RHS = RHS + Diag;
    elseif strcmpi(inletBC,'zero_gradient')
        index_set = [u_in_i v_in_i w_in_i rho_in_i T_in_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        
        [RHS_i,RHS_j,RHS_val] = find(RHS);
        DZ_i    = zeros(4*length(index_set),1); % 5-point stencil with 4 entries
        DZ_j    = DZ_i;
        DZ_val  = DZ_i;
        for row_i=1:length(index_set);
            DZ_i((row_i-1)*4+1:row_i*4)    = [index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
            [~,jj,val] = find(DZ(index_set(row_i),:));
            DZ_j((row_i-1)*4+1:row_i*4)    = jj;
            DZ_val((row_i-1)*4+1:row_i*4)  = val;
        end
        % append derivative matrix elements
        RHS_i       = [RHS_i; DZ_i];
        RHS_j       = [RHS_j; DZ_j];
        RHS_val     = [RHS_val; DZ_val];
        RHS = sparse(RHS_i,RHS_j,RHS_val);
        %     RHS(  u_in_i,   u_j)   = DZ(  u_in_i,   u_j);
        %     RHS(  v_in_i,   v_j)   = DZ(  v_in_i,   v_j);
        %     RHS(  w_in_i,   w_j)   = DZ(  w_in_i,   w_j);
        %     RHS(  rho_in_i, rho_j) = DZ(rho_in_i, rho_j);
        %     RHS(  T_in_i,   T_j)   = DZ(  T_in_i,   T_j);
    else
        % inlet: zero value
        index_set = [u_in_i v_in_i w_in_i T_in_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        RHS = RHS + Diag;
    end
    
    % outlet
    if strcmpi(outletBC,'zero_gradient')
        index_set = [u_out_i v_out_i w_out_i rho_out_i T_out_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        
        [RHS_i,RHS_j,RHS_val] = find(RHS);
        DZ_i    = zeros(4*length(index_set),1); % 5-point stencil with 4 entries
        DZ_j    = DZ_i;
        DZ_val  = DZ_i;
        for row_i=1:length(index_set);
            DZ_i((row_i-1)*4+1:row_i*4)    = [index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
            [~,jj,val] = find(DZ(index_set(row_i),:));
            DZ_j((row_i-1)*4+1:row_i*4)    = jj;
            DZ_val((row_i-1)*4+1:row_i*4)  = val;
        end
        % append derivative matrix elements
        RHS_i       = [RHS_i; DZ_i];
        RHS_j       = [RHS_j; DZ_j];
        RHS_val     = [RHS_val; DZ_val];
        RHS = sparse(RHS_i,RHS_j,RHS_val);
        %     RHS(  u_out_i,   u_j)   = DZ(  u_out_i,   u_j);
        %     RHS(  v_out_i,   v_j)   = DZ(  v_out_i,   v_j);
        %     RHS(  w_out_i,   w_j)   = DZ(  w_out_i,   w_j);
        %     RHS(  rho_out_i, rho_j) = DZ(rho_out_i, rho_j);
        %     RHS(  T_out_i,   T_j)   = DZ(  T_out_i,   T_j);
    else
        index_set = [u_out_i v_out_i w_out_i T_out_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)=1;
        LHS = LHS + Diag;
        RHS = RHS + Diag;
    end
    if pipeBC
        % Nozzle
        index_set = [u_pipe_i v_pipe_i w_pipe_i];% T_pipe_i];
        LHS(index_set, :)   = 0;
        RHS(index_set, :)   = 0;
        diag_index_set      = sub2ind(size_LHS, index_set, index_set);
        Diag                = speye(5*NrNz)*0;
        Diag(diag_index_set)= 1;
        LHS = LHS + Diag;
        RHS = RHS + Diag;
    end
end

time    = toc;
disp(['    elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

clear D2* Dr* Dz* DR DZ Diag DB* RHS_* DZ_* DR_*

%%%%%%%%%%%%%
% Solve EVP %
%%%%%%%%%%%%%

% undo Fourier ansatz in time if dq/dt = L0*q is needed, else keep (L0-omega*I)*q = 0 form as EVP
if strcmpi(calcMethod,'timestepper')||strcmpi(calcMethod,'iterate_optimal')||strcmpi(calcMethod,'frequency_response_PSD')||strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'les_forcingFFT_response')||strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'pseudospectrum')||strcmpi(calcMethod,'calculate_forcings')
    RHS = RHS/1i;
end

% reduce to generalized to standard EVP
L0               = -RHS\LHS;
clear LHS RHS
%%
for SIGMA =  reshape(SIGMA_vec,1,size(SIGMA_vec,1)*size(SIGMA_vec,2)) % makes sure that SIGMA_vec can be a matrix (pseudospectrum)
    % set up save file name
    [saveFile] = getSaveFileName(saveFolder,custom_comment,calcMethod,Re,m,Nr,Nz,r_farf,z1,z2,SIGMA);
    % check if case is already calculating or already exists, and skip if so
    if (~exist(saveFile,'file')&&~exist([saveFile '.lock'],'file'))||(~checkLock)||strcmpi(calcMethod,'pseudospectrum')
        if (checkLock)&&~strcmpi(calcMethod,'pseudospectrum'), save([saveFile '.lock'],'SIGMA'); end
        
        if strcmpi(calcMethod,'calculate_forcings')
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % forcings from LES data %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            f_tilde = zeros(5*Nr*Nz,nblks);
            q_tilde = zeros(5*Nr*Nz,nblks);
            I       = speye(5*NrNz);
            St_i    = findnearest(St_LES,SIGMA/2/pi);
                       
            for iblk = 1:nblks
                disp(['    approximating forces from LES data for block ' num2str(iblk) '/' num2str(nblks)]);
                file_FFT    = matfile([tfft_root '/FFT_m' num2str(m) '_block' num2str(iblk,'%.4i') '.mat']);
                rho_tilde   = interp2(z_LES,r_LES,squeeze(file_FFT.rho_tilde(St_i,:,:)),z,r,'linear',0);
                u_x_tilde   = interp2(z_LES,r_LES,squeeze(file_FFT.u_x_tilde(St_i,:,:)),z,r,'linear',0);
                u_r_tilde   = interp2(z_LES,r_LES,squeeze(file_FFT.u_r_tilde(St_i,:,:)),z,r,'linear',0);
                u_t_tilde   = interp2(z_LES,r_LES,squeeze(file_FFT.u_t_tilde(St_i,:,:)),z,r,'linear',0);
                T_tilde     = interp2(z_LES,r_LES,squeeze(file_FFT.T_tilde(St_i,:,:)),z,r,'linear',0);
                              
                q    = [reshape(rho_tilde,1,NrNz) reshape(u_x_tilde,1,NrNz) reshape(u_r_tilde,1,NrNz) reshape(u_t_tilde,1,NrNz) reshape(T_tilde,1,NrNz)];                
                f    = (L0-1i*SIGMA*I)*transpose(q);
                
                q_tilde(:,iblk)   = q;
                f_tilde(:,iblk)   = f;
            end

            f_tilde(isnan(f_tilde)) = 0;
            %%           
            fig_modes = figure('name',['Forcings for St=' num2str(St_LES(St_i))]);
            subplot(5,2,1)
            pcolor(z,r,reshape(mean(abs(f_tilde(w_j,:)),2),Nr,Nz)); shading interp, axis equal tight, colorbar, colormap(jetwhite)
            
            % forcing SPOD
            S_ff        = f_tilde'*f_tilde; %K'*(M.*K);
            S_ff        = S_ff./nblks;
            S_ff        = 0.5*(S_ff+S_ff');                   % unneseccary but cleans up roundoff issue
            [a,lambda]  = eig(S_ff);
            lambda      = diag(lambda);
            [lambda, sort_i] = sort(lambda,'descend');
            a           = a(:,sort_i);
            f_POD       = f_tilde*a;
            for i=1:nblks
                a(:,i)      = a(:,i)  *sqrt(lambda(i))*sqrt(nblks); % normalize such that mean(a(i)*a(j)) = lambda(i)*delta(i,j);
                f_POD(:,i)  = f_POD(:,i)/sqrt(lambda(i))/sqrt(nblks); % normalize such that q(t) = sum_i=1:nblks(a(i)*psi(i)), or equiv psi(:,i)'*(M(:,1).*psi(:,i)) = 1 (orthonormal in M-induced norm)
            end
            
            for modei = 1:4
                subplot(5,2,modei*2+1)
                pcolor(z,r,reshape(real(f_POD(w_j,modei)),Nr,Nz)); shading interp, axis equal tight, colorbar, colormap(jetblack), caxis([-1 1]*min(abs(caxis)))
            end
            
            fig_energy = figure;
            semilogy(cumsum(lambda/sum(lambda)),'r+-'), hold on
            
            %%           
            figure(fig_modes);
            subplot(5,2,2)
            pcolor(z,r,reshape(mean(abs(q_tilde(w_j,:)),2),Nr,Nz)); shading interp, axis equal tight, colorbar, colormap(jetwhite)
            
            % forcing SPOD
            S_qq        = q_tilde'*q_tilde; %K'*(M.*K);
            S_qq        = S_qq./nblks;
            S_qq        = 0.5*(S_qq+S_qq');                   % unneseccary but cleans up roundoff issue
            [a,lambda]  = eig(S_qq);
            lambda      = diag(lambda);
            [lambda, sort_i] = sort(lambda,'descend');
            a           = a(:,sort_i);
            q_POD       = q_tilde*a;
            for i=1:nblks
                a(:,i)      = a(:,i)  *sqrt(lambda(i))*sqrt(nblks); % normalize such that mean(a(i)*a(j)) = lambda(i)*delta(i,j);
                q_POD(:,i)  = q_POD(:,i)/sqrt(lambda(i))/sqrt(nblks); % normalize such that q(t) = sum_i=1:nblks(a(i)*psi(i)), or equiv psi(:,i)'*(M(:,1).*psi(:,i)) = 1 (orthonormal in M-induced norm)
            end
            
            for modei = 1:4
                subplot(5,2,modei*2+2)
                pcolor(z,r,reshape(real(q_POD(w_j,modei)),Nr,Nz)); shading interp, axis equal tight, colorbar, colormap(jetblack), caxis([-1 1]*min(abs(caxis)))
            end
            
            figure(fig_energy)
            semilogy(cumsum(lambda/sum(lambda)),'b+-')            
        %%    
        elseif strcmpi(calcMethod,'adjoint')
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjoint               %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['    calling Arnoldi (direct LU) of ADJOINT with SIGMA=' num2str(SIGMA)]);
            tic
            
            % weights induced by compressible energy norm
            intWeight   = reshape(trapzWeightsPolar(r_1D,z_1D),NrNz,1);
            vol         = pi*r_farf^2*(z2-z1);
            F  = ...
                1/vol*0.5*[
                (T./(RHO*kappa*Ma^2)).*intWeight
                RHO.*intWeight
                RHO.*intWeight
                RHO.*intWeight
                (RHO./(kappa*(kappa-1)*T*Ma.^2)).*intWeight
                ];
            invF    = 1./F;
            L0adj   = diag(sparse(invF))*ctranspose(L0)*diag(sparse(F));
            
            [eigVecs, OMEGA] = eigs(L0adj, noEigs, SIGMA);
            OMEGA            = conj(diag(OMEGA));
            time = toc;
            disp(['    elapsed time - Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);
            clear L0adj weight F invF
        elseif strcmpi(calcMethod,'timestepper')
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Timestepper           %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['    calling timestepper-based matrix-free Arnoldi with SIGMA=' num2str(SIGMA)]);
            tic
            fig_timestepper = figure('name','Timestepper','position',[392 470 1378 387],'PaperSize',[1378 387],'PaperPositionMode','auto','InvertHardcopy', 'off','Renderer','painters');
            % construct initial condition
            rho_t0  = zeros(Nr*Nz,1);
            u_t0    = zeros(Nr*Nz,1);
            v_t0    = zeros(Nr*Nz,1);
            w_t0    = zeros(Nr*Nz,1);
            T_t0    = zeros(Nr*Nz,1);
            [zz, rr]  = meshgrid(linspace(z_1D(1),z_1D(end),Nz_equi), linspace(r_1D(1),r_1D(end),round(Nz_equi*r_farf/z2)));
            pert_t0   = rand(size(zz))-0.5;
            pert_t0(rr>r_sponge|zz>z_sponge2|zz<z_sponge1) = 0;
            for smooth_i = 1:Nsmooth_pert
                pert_t0(rr==r_1D(1)|rr>=r_sponge|zz>=z_sponge2|zz<=z_sponge1) = 0;   % no perturbation on centerline or in sponge region
                pert_t0  = triSmooth(pert_t0,1,1);
                pert_t0  = triSmooth(pert_t0,2,1);
                pert_t0(rr==r_1D(1)|rr>=r_sponge|zz>=z_sponge2|zz<=z_sponge1) = 0;
            end
            pert_t0  = interp2(zz,rr,pert_t0,z,r,'spline');
            pert_t0  = reshape(pert_t0,Nr*Nz,1);
            pert_t0  = pert_t0/max(abs(pert_t0));
            q_init   = [pert_t0; u_t0; v_t0; w_t0; T_t0];
            clear zz rr pert_t0 smooth_i
            
            % prepare probing
            Nsamples = length(0:dt_probe:t_end);
            Nprobes  = length(probe_z);
            w_probe  = NaN(Nprobes,Nsamples);
            time_probe  = NaN(Nsamples);
            St       = (0:nfft_probe/2-1)/Ma/nfft_probe/dt_probe;
            sample_i = 1;
            
            time_vec   = [0 t_end 0:dt_filter:t_end 0:dt_vis:t_end 0:dt_probe:t_end]; time_vec = sort(unique(round(time_vec,10)));
            wp_max     = zeros(length(0:dt_vis:t_end),1);
            time_vis   = zeros(length(0:dt_vis:t_end),1);
            vis_i      = 1;
            
            global calls
            calls   = 0;
            
            figure(fig_timestepper), drawnow
            LNSEfun = @(t,q1)LNSE(t,q1,L0);
            for t_i = 2:length(time_vec)   % clumsy because time integration gets interrupted but needed for filtering
                time    = time_vec(t_i);
                options = odeset('Stats','off');
                [~,q2]  = ode45(LNSEfun,linspace(time_vec(t_i-1),time,3),q_init,options);
                q2      = squeeze(q2(3,:));
                
                % filter primitive perturbation variables
                if mod(time,dt_filter)==0
                    alpha_f    = 0.4;
                    q2(rho_j)  = reshape(filter_x(reshape(q2(rho_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(u_j)    = reshape(filter_x(reshape(q2(u_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(v_j)    = reshape(filter_x(reshape(q2(v_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(w_j)    = reshape(filter_x(reshape(q2(w_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(T_j)    = reshape(filter_x(reshape(q2(T_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(rho_j)  = reshape(filter_y(reshape(q2(rho_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(u_j)    = reshape(filter_y(reshape(q2(u_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(v_j)    = reshape(filter_y(reshape(q2(v_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(w_j)    = reshape(filter_y(reshape(q2(w_j),Nr,Nz),alpha_f),Nr*Nz,1);
                    q2(T_j)    = reshape(filter_y(reshape(q2(T_j),Nr,Nz),alpha_f),Nr*Nz,1);
                end
                
                % probe solution
                if mod(time,dt_probe)==0
                    for probe_i = 1:Nprobes
                        w_probe(probe_i,sample_i)   = interp2(z, r, reshape(q2(w_j),Nr,Nz), probe_z(probe_i), probe_r(probe_i));
                        time_probe(sample_i)        = time;
                        if mod(sample_i,nfft_probe)==0
                            w_hat   = fft(w_probe(probe_i, sample_i-nfft_probe+1:sample_i));
                            PSD     = fftshift(abs(w_hat).^2/nfft_probe);
                            PSDdB   = 10*log10(PSD);
                        end
                        sample_i    = sample_i + 1;
                    end
                end
                
                % visualize solution
                if mod(time,dt_vis)==0
                    figure(fig_timestepper)
                    subplot(2,2,1:2)
                    wp_max(vis_i)   = max(abs(q2(w_j(:)))); time_vis(vis_i) = time;
                    wp              = reshape(q2(w_j),Nr,Nz)/wp_max(vis_i);
                    %pcolor_h        = pcolor(z,r,real(wp)); set(pcolor_h,'EdgeColor','none'), caxis([-1 1])
                    contourf(z,r,real(wp),linspace(-1,1,30),'EdgeColor','none');
                    hold on
                    for probe_i = 1:Nprobes, plot(probe_z(probe_i), probe_r(probe_i),'+'), text(probe_z(probe_i), probe_r(probe_i), num2str(probe_i),'VerticalAlignment','Bottom','HorizontalAlignment','left','BackgroundColor','none','FontSize',6), end
                    axis equal, axis tight, hold off
                    colormap(gray)
                    xlabel('$$x$$'); ylabel('$$r$$'); title(['$$t=$$' num2str(time) ' (call No. ' num2str(calls) ')'], 'Interpreter','latex'); caxis([-1 1])
                    c = colorbar('east'); c.Label.String  = '$$\frac{u''_x}{|u''_x|_\infty}$$'; c.TickLabelInterpreter = 'latex';
                    subplot(2,2,3)
                    semilogy(time_vis,wp_max,'+-'); xlabel('$$t$$'); ylabel('$$|u''_x|_\infty$$'); hold off
                    subplot(2,2,4)
                    if time<=nfft_probe*dt_probe
                        plot(time_probe,real(w_probe(probe_i,:))); xlabel('$$t$$'); ylabel('$$u''_x|_{probe 1}$$'); hold off
                    else
                        stem(St, PSD(nfft_probe/2+1:end)); ylabel('PSD$$_{u_x,probe 1}$$'); xlabel('St');
                    end
                    drawnow, pause(0.1)
                    print(['animation/' custom_comment '_frame' num2str(t_i-1,'%.5d')],'-dpng','-r200')
                    drawnow, pause(0.1)
                    vis_i      = vis_i + 1;
                end
                
                q_init  = q2;
            end
            
            % just for visulaisation routine
            OMEGA   = 0;
            eigVecs = permute(q2,[2 1]); clear q2 wp
            
            time = toc;
            disp(['    elapsed time - timestepper based Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);
        elseif strcmpi(calcMethod,'iterate_optimal')
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Optimal perturbation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['    calling timestepper with SIGMA=' num2str(SIGMA)]);
            tic
            
            %         % weights induced by compressible energy norm
            %         intWeight   = reshape(trapzWeightsPolar(r_1D,z_1D),NrNz,1);
            %         vol         = pi*r_farf^2*(z2-z1);
            %         F  = ...
            %             1/vol*0.5*[
            %             (T./(RHO*kappa*Ma^2)).*intWeight
            %             RHO.*intWeight
            %             RHO.*intWeight
            %             RHO.*intWeight
            %             (RHO./(kappa*(kappa-1)*T*Ma.^2)).*intWeight
            %             ];
            %         invF    = 1./F;
            %         L0adj = diag(sparse(invF))*ctranspose(L0)*diag(sparse(F));
            % CHECK IF CHOLESKY OR NOT!!!
            
            rho_t0  = zeros(Nr*Nz,1);
            u_t0    = zeros(Nr*Nz,1);
            v_t0    = zeros(Nr*Nz,1);
            w_t0    = zeros(Nr*Nz,1);
            T_t0    = zeros(Nr*Nz,1);
            
            NTriSmth = 50;
            pert_t0  = rand(Nr,Nz)-0.5;
            pert_t0  = triSmooth(pert_t0,1,NTriSmth);
            pert_t0  = triSmooth(pert_t0,2,1);
            pert_t0  = pert_t0.*(1-reshape(spongeFun,Nr,Nz));
            pert_t0  = reshape(pert_t0,Nr*Nz,1);
            
            q_init  = [pert_t0; u_t0; v_t0; w_t0; T_t0];
            
            calls = 0;
            
            for no_loops = 1:10
                % direct solution
                tic
                disp(' ')
                disp('time integration')
                disp('----------------')
                tspan   = linspace(0,t_end,50);
                LNSEfun = @(t,q)LNSE(t,q,L0);
                options = odeset('Stats','on');
                [~,q_soln] = ode45(LNSEfun,tspan,q_init,options);
                time = toc;
                disp(['    elapsed time - time stepping: ' datestr(time/24/3600, 'HH:MM:SS')]);
                
                % adjoint solution
                q_init = q_soln(end,:); % set direct calculation output as initial condition for adjoint
                LNSEfun = @(t,q)LNSE(t,q,L0adj);
                tic
                disp(' ')
                disp('adjoint time integration')
                disp('------------------------')
                options = odeset('Stats','on');
                [~,q_soln_adj] = ode45(LNSEfun,tspan,q_init,options);
                time = toc;
                disp(['    elapsed time - adjoint time stepping: ' datestr(time/24/3600, 'HH:MM:SS')]);
                q_init = q_soln_adj(end,:); % set dajoint calculation output as initial condition for direct
            end
            
            % combine solutions for animation
            q_soln  = [q_soln; q_soln_adj];
            tspan   = [tspan fliplr(tspan)];
            
            while 1
                cmax = max(abs(pert_t0(:))); cmin = -cmax;
                figure, contourf(z, r, reshape(real(q_soln(1,w_j)),Nr,Nz), 20, 'edgecolor','none'); axis equal, axis tight
                caxis([cmin cmax])
                drawnow
                input('start animation?')
                for ti = 1:length(tspan)
                    f = real(reshape(q_soln(ti,w_j), Nr, Nz));
                    cmax = max(abs(f(:))); cmin = -cmax;
                    contourf(z, r, f, 20, 'edgecolor','none'); axis equal, axis tight
                    caxis([cmin cmax]);
                    drawnow
                    pause(0.01)
                end
            end
            
            vizBaseAndRes = false;
            
        elseif (strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'frequency_response_PSD')||strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response'))
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Input-output analysis %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % see Nichols, CTR report 2014 "Input-output analysis of high-speed jet noise"
            disp(['    calling intput-output analysis with SIGMA=' num2str(SIGMA)]);
            tic
            
            % weights induced by compressible energy norm
            intWeight   = reshape(trapzWeightsPolar(r_1D,z_1D),NrNz,1);
            vol         = pi*r_farf^2*(z2-z1);
            F  = ...
                [
                (T./(RHO*kappa*Ma^2)).*intWeight
                RHO.*intWeight
                RHO.*intWeight
                RHO.*intWeight
                (RHO./(kappa*(kappa-1)*T*Ma.^2)).*intWeight
                ];
            
            invF    = 1./F;
            F       = diag(sparse(F));
            invF    = diag(sparse(invF));
            
            % get input and output matrices from user-defined routine
            [B,C] = setUpBandC(r_farf,Nr,z1,z2,Nz,r,z,Nsponge_left,Nsponge_right,Nsponge_top,spongeFun,RHO,U,V,W,T,MU,T_0,rho_0,mu_0,p_0,nu_0,Rgas,kappa,Pr,c,Re,Ma,u_0,OUTPUT_z_min,OUTPUT_z_max,OUTPUT_r_min,OUTPUT_r_max,INPUT_z_min,INPUT_z_max,INPUT_r_min,INPUT_r_max,forcingType,responseType,baseState,vizBaseAndRes);
            
            % LU-decomposition of L-sigma*I
            tic
            LsI                 = L0-1i*SIGMA*speye(5*Nr*Nz);
            [LL,UU,pp,qq,rr]    = lu(LsI);
            time = toc;
            disp(['    elapsed time - LU-decomposition of L-sigma*I: ' datestr(time/24/3600, 'HH:MM:SS')]);
            
            % complex conjugate transpose
            LL_ct   = LL';
            UU_ct   = UU';
            pp_ct   = pp';
            qq_ct   = qq';
            rr_ct   = rr';
            C_ct    = C';
            B_ct    = B';
            
            %%
            
            if filter_f
                %%
                disp('--> filtering forcing')
                Z       = 0*speye(NrNz,NrNz);
                alpha_f = 0.4;
                %                 [Mz_1D,Fz_1D] = filter_O10_diffMat(Nz,alpha_f);
                %                 [Mr_1D,Fr_1D] = filter_O10_diffMat(Nr,alpha_f);
                %                 Mz      = sparse(kron(Mz_1D,speye((Nr))));
                %                 Mr      = sparse(kron(speye((Nz)),Mr_1D));
                %                 Fz      = sparse(kron(Fz_1D,speye((Nr))));
                %                 Fr      = sparse(kron(speye((Nz)),Fr_1D));
                
                filter_z_i_start    = findnearest(z_1D,z_sponge1,2);
                filter_z_i_end      = findnearest(z_1D,z_sponge2,1);
                filter_r_i_start    = 1;
                filter_r_i_end      = findnearest(r_1D,r_sponge,1);
                
                Nz_filter           = filter_z_i_end - filter_z_i_start + 1;
                Nr_filter           = filter_r_i_end - filter_r_i_start + 1;
                
                alpha_f = 0.4;
                
                [Mz_1D,Fz_1D] = filter_O10_diffMat(Nz_filter,alpha_f);
                tmp       = eye(Nz);
                tmp(filter_z_i_start:filter_z_i_start+Nz_filter-1,filter_z_i_start:filter_z_i_start+Nz_filter-1) = Mz_1D;
                Mz_1D   = tmp;
                tmp(filter_z_i_start:filter_z_i_start+Nz_filter-1,filter_z_i_start:filter_z_i_start+Nz_filter-1) = Fz_1D;
                Fz_1D   = tmp;
                [Mr_1D,Fr_1D] = filter_O10_diffMat(Nr_filter,alpha_f);
                tmp       = eye(Nr);
                tmp(filter_r_i_start:filter_r_i_start+Nr_filter-1,filter_r_i_start:filter_r_i_start+Nr_filter-1) = Mr_1D;
                Mr_1D   = tmp;
                tmp(filter_r_i_start:filter_r_i_start+Nr_filter-1,filter_r_i_start:filter_r_i_start+Nr_filter-1) = Fr_1D;
                Fr_1D   = tmp;
                
                Mz      = speye(Nr*Nz);
                Mr      = speye(Nr*Nz);
                Fz      = speye(Nr*Nz);
                Fr      = speye(Nr*Nz);
                
                % filter in z (at one r-index at a time)
                for r_i = filter_r_i_start:filter_r_i_start+Nr_filter-1
                    idx             = r_i:Nr:r_i+(Nz-1)*Nr;
                    Mz(idx,idx)     = Mz_1D;
                    Fz(idx,idx)     = Fz_1D;
                end
                
                % filter in r (at one z-index at a time)
                for z_i = filter_z_i_start:filter_z_i_start+Nz_filter-1
                    idx             = (z_i-1)*Nr+1:1:(z_i-1)*Nr+Nr;
                    Mr(idx,idx)     = Mr_1D;
                    Fr(idx,idx)     = Fr_1D;
                end
                
                MZ  = sparse([ ...
                    Mz        Z         Z           Z        Z
                    Z         Mz        Z           Z        Z
                    Z         Z         Mz          Z        Z
                    Z         Z         Z           Mz       Z
                    Z         Z         Z           Z        Mz
                    ]);
                MR  = sparse([ ...
                    Mr        Z         Z           Z        Z
                    Z         Mr        Z           Z        Z
                    Z         Z         Mr          Z        Z
                    Z         Z         Z           Mr       Z
                    Z         Z         Z           Z        Mr
                    ]);
                FZ  = sparse([ ...
                    Fz        Z         Z           Z        Z
                    Z         Fz        Z           Z        Z
                    Z         Z         Fz          Z        Z
                    Z         Z         Z           Fz       Z
                    Z         Z         Z           Z        Fz
                    ]);
                FR  = sparse([ ...
                    Fr        Z         Z           Z        Z
                    Z         Fr        Z           Z        Z
                    Z         Z         Fr          Z        Z
                    Z         Z         Z           Fr       Z
                    Z         Z         Z           Z        Fr
                    ]);
                clear Z Mz Mr Fz Fr Mz_1D Mr_1D Fz_1D Fr_1D
                
                MR_ct   = MR';
                FR_ct   = FR';
                MZ_ct   = MZ';
                FZ_ct   = FZ';
            end
            
            %%
            % anonymous function that applies H*ctranspose(H) on vector w
            % NOTE: ctranspose(H)*H gives right-singular vectors    -> input basis V
            %       H*ctranspose(H) gives left-singular vectors     -> output basis U
            if filter_f
                HtrH    = @(v) invF*(   FZ_ct*(MZ_ct\(FR_ct*(MR_ct\(   B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(F*(C*(qq*(UU\(LL\(pp*(rr\(B*    (MR\(FR*(MZ\(FZ*    v)))))))))))))))))))))));% CM^-1BS (f only)
            else
                HtrH    = @(v) invF*(B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(F*(C*(qq*(UU\(LL\(pp*(rr\(B*v)))))))))))))));
            end
            
            %%
            opts.isreal         = 0;
            opts.tol            = eps;
            opts.disp           = 3;
            if exist('V_in','var')
                opts.v0         = V_in(:,1);
            end
            
            % input
            tic
            if ~(strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response'))
                [V_in, OMEGA]       = eigs(HtrH, 5*Nz*Nr, noEigs, 'lm', opts);
                gain                = sqrt(real(diag(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))
            else
                OMEGA               = eye(noEigs);
            end
            
            time = toc;
            disp(['    elapsed time - SVD via ''eigs'' using matrix-free Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);            %%
            % output
            U_out = zeros(size(V_in));
            for eig_i = 1:noEigs
                % normalization in inner product norm
                V_in(:,eig_i)  = V_in(:,eig_i)/sqrt(V_in(:,eig_i)'*F*V_in(:,eig_i));
                if strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response')
                    gain(eig_i)     = real(sqrt(V_in(:,eig_i)'*(B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(F*(C*(qq*(UU\(LL\(pp*(rr\(B*V_in(:,eig_i)))))))))))))))))); % real() because of O(eps) imaginary error
                end
                U_out(:,eig_i) = C*(qq*(UU\(LL\(pp*(rr\(B*(V_in(:,eig_i)/gain(eig_i))))))));
            end
            
            if strcmpi(calcMethod,'frequency_response_PSD')
                p       = (U_out(T_j,:).*RHO + U_out(rho_j,:).*T).*gain'/kappa/Ma^2;
                p_sum   = sum(p,2); 
                psd     = sum(abs(p),2); 
                % psd     = reshape(sum(abs(p),2),Nr,Nz);                
                U_out   = [psd p_sum];%psd; % only save low-rank PSD approximation
                V_in    = V_in(:,1);
            end
                       
            % % check orthogonality
            % abs(V_in'*V_in)         % I/O should be orthonormal in 2-norm
            % abs(V_in'*F*V_in)       % frequency response should be orthonormal in (induced) energy norm
            % abs(U_out'*U_out)
            % abs(U_out'*F*U_out)
            %%
            clear  HtrH LsI LL  UU pp qq rr C B LL_ct UU_ct pp_ct qq_ct rr_ct C_ct B_ct MR_ct MR FR_ct FR MZ_ct MZ FZ_ct FZ F invF        
            OMEGA               = diag(OMEGA);
            %%
            if (vizBaseAndRes)
                %%
                figure('name','gain factors');
                stem(gain);
                xlabel('mode number N'); ylabel('singular value $$\sigma$$');
                
                fig_f_p   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - forcing pressure'],'fileName',['resolvent_f_pressure_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                fig_f_w   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - forcing u_x'],'fileName',['resolvent_f_u_x_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                fig_q_p   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - response pressure'],'fileName',['resolvent_q_pressure_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                fig_q_w   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - response u_x'],'fileName',['resolvent_q_u_x_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                
                if noEigs<5, n_row = noEigs; else n_row = 5; end
                n_col = ceil(noEigs/5);
                for i = 1:noEigs
                    %i=i+5;
                    rho_f   = reshape(V_in(rho_j,i), Nr, Nz);
                    w_f     = reshape(V_in(w_j,i), Nr, Nz);
                    t_f     = reshape(V_in(T_j,i), Nr, Nz);
                    p_f     = (reshape(RHO, Nr, Nz).*t_f + rho_f.*reshape(T, Nr, Nz))/kappa/Ma^2;
                    rho_q   = reshape(U_out(rho_j,i), Nr, Nz);
                    w_q     = reshape(U_out(w_j,i), Nr, Nz);
                    t_q     = reshape(U_out(T_j,i), Nr, Nz);
                    p_q     = (reshape(RHO, Nr, Nz).*t_q + rho_q.*reshape(T, Nr, Nz))/kappa/Ma^2;
                    %i=i-5;
                    figure(fig_f_p)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(p_f)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1])
                    hold on
                    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
                    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
                    %                     rectangle('position',[INPUT_z_min INPUT_r_min (INPUT_z_max-INPUT_z_min) (INPUT_r_max-INPUT_r_min)],'edgecolor','b','LineStyle','--')
                    %                     rectangle('position',[z_sponge1 -0.01 (z_sponge2-z_sponge1) r_sponge],'edgecolor','m','LineStyle','-')
                    ylabel('$r$');
                    % xlim([z_sponge1 z_sponge2]);
                    % ylim([r(1) r_sponge]);
                    
                    figure(fig_f_w)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(w_f)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1])
                    hold on
                    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
                    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
                    %                     rectangle('position',[INPUT_z_min INPUT_r_min (INPUT_z_max-INPUT_z_min) (INPUT_r_max-INPUT_r_min)],'edgecolor','b','LineStyle','--')
                    %                     rectangle('position',[z_sponge1 -0.01 (z_sponge2-z_sponge1) r_sponge],'edgecolor','m','LineStyle','-')
                    ylabel('$r$');
                    % xlim([z_sponge1 z_sponge2]);
                    % ylim([r(1) r_sponge]);
                    
                    figure(fig_q_p)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(p_q)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1])
                    hold on
                    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
                    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
                    %                     rectangle('position',[OUTPUT_z_min OUTPUT_r_min (OUTPUT_z_max-OUTPUT_z_min) (OUTPUT_r_max-OUTPUT_r_min)],'edgecolor','b','LineStyle','--')
                    %                     rectangle('position',[z_sponge1 -0.01 (z_sponge2-z_sponge1) r_sponge],'edgecolor','m','LineStyle','-')
                    ylabel('$r$');
                    % xlim([z_sponge1 z_sponge2]);
                    % ylim([r(1) r_sponge]);
                    
                    figure(fig_q_w)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(w_q)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1])
                    hold on
                    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
                    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
                    %                     rectangle('position',[OUTPUT_z_min OUTPUT_r_min (OUTPUT_z_max-OUTPUT_z_min) (OUTPUT_r_max-OUTPUT_r_min)],'edgecolor','b','LineStyle','--')
                    %                     rectangle('position',[z_sponge1 -0.01 (z_sponge2-z_sponge1) r_sponge],'edgecolor','m','LineStyle','-')
                    ylabel('$r$');
                    % xlim([z_sponge1 z_sponge2]);
                    % ylim([r(1) r_sponge]);
                end
                
                this_fig = {fig_f_p,fig_f_w,fig_q_p,fig_q_w};
                for figi = 1:4
                    figure(this_fig{figi})
                    xlabel('$x$');
                    colormap(gray)
                    fig_pos = get(gcf,'position'); set(gcf,'Position',[fig_pos(1) fig_pos(2) n_col*fig_pos(3) fig_pos(3)*n_row*((1.3*r_sponge)/z_sponge2)]);
                    fig_w   = fig_pos(3); fig_h = fig_pos(4);
                end
                clear fig_f_p fig_f_w fig_q_p fig_q_w rho_f w_f t_f p_f rho_q w_q t_q p_q
                
                vizBaseAndRes   = false;
            end
            time = toc;
            disp(['    elapsed time - intput-output analysis: ' datestr(time/24/3600, 'HH:MM:SS')]);
            
        elseif strcmpi(calcMethod,'pseudospectrum')
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Pseudospectrum        %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            opts.isreal         = 0;
            opts.tol            = eps;
            opts.disp           = 0;
            
            % weights induced by compressible energy norm
            intWeight = reshape(trapzWeightsPolar(r_1D,z_1D),NrNz,1);
            weight  = ...
                0.5*[
                (T./(RHO*kappa*Ma^2)).*intWeight
                RHO.*intWeight
                RHO.*intWeight
                RHO.*intWeight
                (RHO./(kappa*(kappa-1)*T*Ma.^2)).*intWeight
                ];
            % Cholesky decomposition of weight matrix, see e.g. Trefethen "Computation of Pseudospectra", p. 6 (books.bib)
            F    = sqrt(weight);
            invF = 1./F;
            LW   = diag(sparse(F))*L0*diag(sparse(invF));
            clear L0 weight F invF
            
            disp(['    calculating ' num2str(nSIGMA_real) 'x' num2str(nSIGMA_imag) ' pseudo-spectrum']);
            tic
            count = 0;
            for i = 1:length(SIGMA_vec(:))
                SIGMA = SIGMA_vec(i);
                % set up save file name
                [saveFile] = getSaveFileName(saveFolder,custom_comment,calcMethod,Re,m,Nr,Nz,r_farf,z1,z2,SIGMA);
                if (~exist(saveFile,'file')&&~exist([saveFile '.lock'],'file'))
                    save([saveFile '.lock'],'SIGMA');
                    
                    disp(['    -> calculating resolvent norm for SIGMA = ' num2str(SIGMA)]);
                    SIGMA = SIGMA_vec(i);
                    
                    % LU-decomposition of i*sigma*I - LW
                    LsI              = LW - 1i*SIGMA*speye(5*Nr*Nz);
                    [L,U,pp,qq,rr]   = lu(LsI);
                    
                    % complex conjugate transpose
                    L_ct    = L';
                    U_ct    = U';
                    pp_ct   = pp';
                    qq_ct   = qq';
                    rr_ct   = rr';
                    
                    % anonymous function that applies ctranspose(R)*R (R=(LsI)^-1)
                    % NOTE: ctranspose(R)*R gives right-singular vectors    -> input basis V
                    RtrR    = @(v) rr_ct\(pp_ct*(L_ct\(U_ct\(qq_ct*(qq*(U\(L\(pp*(rr\v))))))))); % returns (RtrR)^-1*v
                    
                    [V_in, OMEGA_max]    = eigs(RtrR, 5*Nz*Nr, 1, 'lm', opts);  % largest singular value given by sqrt of largest eigenvalue of R^H*R
                    OMEGA                = 1/sqrt(real(OMEGA_max));             % sqrt of eigenvalue of RtrR gives singular value of R
                    
                    count = count+1;
                    
                    if saveEigVecs
                        U_out                = qq*(U\(L\(pp*(rr\(V_in/OMEGA_max)))));
                    else
                        clear V_in
                    end
                    clear LsI L U pp qq rr L_ct U_ct pp_ct qq_ct rr_ct RtrR
                    if saveSoln
                        disp(['    writing output to ' saveFile]);
                        save(saveFile, '-regexp', '^(?!(LW)$).');
                    end
                    delete([saveFile '.lock']);
                else
                    disp(['--> skipping SIGMA=' num2str(SIGMA) ' bacause case is running or finished!']);
                end
            end
            
            time = toc;
            disp(['    elapsed time - pseudo-spectrum calculation (' num2str(count) ' points): ' datestr(time/24/3600, 'HH:MM:SS')]);
            return
        else
            disp(['    calling Arnoldi (direct LU) with SIGMA=' num2str(SIGMA)]);
            tic
            
            [eigVecs, OMEGA] = eigs(L0, noEigs, SIGMA);
            
            OMEGA            = diag(OMEGA);
            time = toc;
            disp(['    elapsed time - Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);
        end
        
        % saving
        if (saveSoln)
            clear count
            disp(['    writing output to ' saveFile]);
            if (saveEigVecs)
                save(saveFile, '-regexp', '^(?!(L0)$).');
            else
                save(saveFile, '-regexp', '^(?!(L0|eigVecs|V_in|U_out)$).');
            end
        end
        
        disp(['    done with  ' num2str(calcCount) '/' num2str(Ncalcs) '!']);
        calcCount = calcCount + 1;
        if exist([saveFile '.lock'],'file'), delete([saveFile '.lock']); end % unlock case
    elseif exist([saveFile '.lock'],'file')
        disp(['--> NOTE!: case ' saveFile ' is already running, proceeding to next SIGMA...']);
    else
        disp(['--> NOTE!: solution ' saveFile ' already exists, proceeding to next SIGMA...']);
    end
end
clear L0
%%
if (vizBaseAndRes)
    figure('name','spectrum');
    plot(real(OMEGA),imag(OMEGA),'go')
    hold on
    for i=1:noEigs
        text(real(OMEGA(i)),imag(OMEGA(i)),[' ' num2str(i)])
    end
    plot(real(SIGMA),imag(SIGMA),'rx');
    
    nrows = floor(sqrt(noEigs));
    ncols = ceil(noEigs/nrows);
    fig_overview_w = figure('name', 'w''');
    fig_overview_p = figure('name', 'p''');
    for i=1:noEigs
        figure(fig_overview_w)
        subplot(nrows,ncols,i)
        w = real(reshape(eigVecs(w_j,i), Nr, Nz));
        pcolor_h=pcolor(z(1,:), r(:,1), w); set(pcolor_h,'EdgeColor','none'), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[]), axis equal, axis tight
        colormap(gray)
        text(z(1), r(end), ['mode ' num2str(i)], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        [val idx] = max(abs(w(:)));
        caxis([-max(abs(w(idx))) max(abs(w(idx)))]);
        hold on
        figure(fig_overview_p)
        subplot(nrows,ncols,i)
        rho = reshape(eigVecs(rho_j,i), Nr, Nz);
        t   = reshape(eigVecs(T_j,i), Nr, Nz);
        p = (reshape(RHO, Nr, Nz).*t + rho.*reshape(T, Nr, Nz))/kappa/Ma^2;
        pcolor_h=pcolor(z(1,:), r(:,1), real(p)); set(pcolor_h,'EdgeColor','none'), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[]), axis equal, axis tight
        colormap(gray)
        text(z(1), r(end), ['mode ' num2str(i)], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        [val idx] = max(abs(p(:)));
        caxis([-max(abs(p(idx))) max(abs(p(idx)))]);
        hold on
    end
    
    while 1
        i = input('Eigenvector # (or 0 to quit)?');
        if i==0, break, end
        
        % detailed individual mode plot
        figure('name',['mode#' num2str(i)],'position',[2         452        1346         508])
        rho = reshape(eigVecs(rho_j,i), Nr, Nz);
        w = reshape(eigVecs(w_j,i), Nr, Nz);
        u = reshape(eigVecs(u_j,i), Nr, Nz);
        v = reshape(eigVecs(v_j,i), Nr, Nz);
        t = reshape(eigVecs(T_j,i), Nr, Nz);
        p = (reshape(RHO, Nr, Nz).*t + rho.*reshape(T, Nr, Nz))/kappa/Ma^2;
        sponge = reshape(spongeFun, Nr, Nz);
        
        [~,idx]       = max(abs(w(:)));
        [z_idx,r_idx]   = ind2sub(size(w),idx);
        % contourplot streamwise plane
        subplot(2,4,1:3)
        plotvar = real(p);
        pcolor_h=pcolor(z, r, plotvar); set(pcolor_h,'EdgeColor','none'), axis tight, hold on
        plot([z(1,r_idx) z(1,r_idx)],[r(1) r(end)],'r:')
        plot([z(1) z(end)],[r(z_idx,1) r(z_idx,1)],'r:')
        contour(z, r, sponge, [10*eps 10*eps], 'k:');
        caxis([-0.1*(norm(plotvar(:),Inf)) 0.1*(norm(plotvar(:),Inf))]);
        colormap(gray)
        c = colorbar('east');
        c.Label.String = 'real(p'')';
        xlabel('x'); ylabel('r');
        hold on
        [c1,h1]=contour(z, r, real(w), linspace(-norm(w(:),Inf),norm(w(:),Inf),10), 'edgecolor','r');
        contour_negDashedposSolid(c1,h1)
        % radial profiles
        subplot(2,4,[4 8])
        plot(abs(  w(:,r_idx)), r_1D, 'r-'), hold on
        plot(abs(  u(:,r_idx)), r_1D, 'k:')
        plot(abs(rho(:,r_idx)), r_1D, 'k--')
        plot(abs(  t(:,r_idx)), r_1D, 'k-.')
        plot(abs(  p(:,r_idx)), r_1D, 'b-')
        % plot(abs(v(:,r_idx)), r_1D, 'k.')
        legend('|u''|','|v''|','|\rho''|','|T''|','|p''|')
        xlabel('r'), ylabel('amplitude'), axis tight
        % streamwise profiles
        subplot(2,4,5:7)
        plot(z_1D, abs(  w(z_idx,:)), 'r-'), hold on
        plot(z_1D, abs(  real(w(z_idx,:))), 'k-', 'color', [.5 .5 .5])
        plot(z_1D, abs(  u(z_idx,:)), 'k:')
        plot(z_1D, abs(rho(z_idx,:)), 'k--')
        plot(z_1D, abs(  t(z_idx,:)), 'k-.')
        plot(z_1D, abs(  p(z_idx,:)), 'b-')
        legend('|u''|','real(u'')','|v''|','|\rho''|','|T''|','|p''|','location','northeast')
        xlabel('x'), ylabel('amplitude'), axis tight
        
        % pressure and velocity contour plots
        figure('name',['mode#' num2str(i)],'position',[575           2        1345         453])
        subplot(2,1,1)
        plotvar = real(p(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right)); varname = 'p''';
        plotvar = plotvar./norm(plotvar(:),Inf);
        pcolor_h=pcolor(z(1,Nsponge_left+1:end-Nsponge_right), r(1:end-Nsponge_top,1), plotvar); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on
        caxis([-0.5 0.5]);
        colormap(gray)
        c = colorbar('eastoutside');
        c.Label.String = varname;
        xlabel('x'); ylabel('r'); title('pressure');
        subplot(2,1,2)
        plotvar = real(w(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right)); varname = 'u_x''';
        plotvar = plotvar./norm(plotvar(:),Inf);
        pcolor_h=pcolor(z(1,Nsponge_left+1:end-Nsponge_right), r(1:end-Nsponge_top,1), plotvar); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on
        caxis([-0.5 0.5]);
        colormap(gray)
        c = colorbar('eastoutside');
        c.Label.String = varname;
        xlabel('x'); ylabel('r'); title('streamwise velocity');
        drawnow
    end
end

% profile viewer
% p = profile('info');
















