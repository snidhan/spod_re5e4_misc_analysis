if ~exist('Nsmooth_r','var'),   Nsmooth_r   = 0;        end
if ~exist('Nsmooth_z','var'),   Nsmooth_z   = 0;        end
if ~exist('spongeType','var'),  spongeType  = 'poly5';  end
if ~exist('checkLock','var'),   checkLock   = false;    end
if ~exist('Ncalcs','var'),      Ncalcs      = 1;        end
if ~exist('calcCount','var'),   calcCount   = 1;        end
if ~exist('forcingType','var'), forcingType = 'global'; end
if ~exist('filter_f','var'),    filter_f    = false;    end
if ~exist('filter_q','var'),    filter_q    = false;    end
if ~exist('FDorder','var'),     FDorder     = 4;        end
if ~exist('parallelFlow','var'),parallelFlow= false;    end
% list of variables to be saved in .info.mat file (if they exist)
varList ={'r','z','Nr','Nz','RHO','U','V','W','T','P','MU','rho_j','u_j','v_j','w_j','T_j','z_sponge1','z_sponge2','r_sponge', ...
'Nsmooth_r','Nsmooth_z','FDorder','Nr_jet','Nr_shearl','r_farf','r_sponge','r_shearl1','r_shearl2','Nr_smooth','z_segm',...
'Nz_segm','Nz_smooth','outletBC','farfieldBC','inletBC','bottomBC','spongeType','aSponge_left','aSponge_right','aSponge_top',...
'spongeEps_left','spongeEps_right','spongeEps_top','noEigs','SIGMA_vec','filter_f','filter_q','forcingType','INPUT_z_min','INPUT_z_max',...
'INPUT_r_min','INPUT_r_max','OUTPUT_r_min','OUTPUT_r_max','OUTPUT_z_min','OUTPUT_z_max','nSmooth_B_x','nSmooth_B_y','nSmooth_C_x', ...
'nSmooth_C_y','custom_comment','calcMethod','Re','Pr','kappa','Ma','S1','T_ref','mu_ref','z1','z2'};

Nz      = sum(Nz_segm);
NrNz    = Nr*Nz;

disp(['--> calculating for kz=' num2str(kz) ', Nr=' num2str(Nr) ', Nz=' num2str(Nz)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiation Matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% temporary equidistant mesh
dr              = r_farf/(Nr-0.5);
r_1D_equidist   = linspace(dr/2,r_farf,Nr)';
z_1D_equidist   = linspace(z_segm(1),z_segm(end),Nz)';
dz              = z_1D_equidist(2)-z_1D_equidist(1);
tic

% construct transformed coordinates
if (r_shearl1~=0&&r_shearl2~=0)
    dr_jet          = r_shearl1/(Nr_jet-0.5);
    dr_shearl       = (r_shearl2-r_shearl1)/(Nr_shearl-1);
    dr_farf         = (r_farf-r_shearl2)/(Nr-Nr_jet-Nr_shearl-1);
    dr_avg1         = (dr_jet+dr_shearl)/2;
    dr_avg2         = (dr_shearl+dr_farf)/2;
    r_1D            = [linspace(0,r_shearl1-dr_avg1/2,Nr_jet) linspace(r_shearl1+dr_avg1/2,r_shearl2-dr_avg2/2,Nr_shearl) linspace(r_shearl2+dr_avg2/2,r_farf,Nr-Nr_jet-Nr_shearl)];
    for i=1:Nr_smooth, r_1D = mvg(r_1D); end
else
    r_1D        = r_1D_equidist;
    dr_farf     = dr;
end

% if (z_core~=0)
%     dz_core         = (z_core-z1)/(Nz_core-1);
%     dz_rest         = (z2-z_core)/(Nz-Nz_core-1);
%     dz_avg          = (dz_core+dz_rest)/2;
%     z_1D            = [linspace(z1,z_core-dz_avg/2,Nz_core) linspace(z_core+dz_avg/2,z2,Nz-Nz_core)];
%     for i=1:Nz_smooth, z_1D = mvg(z_1D); end
% else
%     z_1D    = z_1D_equidist;
%     dz_rest = dz;
%     dz_core = dz;
% end

% z_segm      = [x0-1 x0 x0+0.5 x1-0.5 x1 x1+2];
% Nz_segm     = [15 30 200 30 30];

N_segm          = length(z_segm)-1;
dz_segm         = diff(z_segm)./Nz_segm;
Nz_cumsum       = cumsum(Nz_segm);

% z_1D(1:Nz_segm(1))  = linspace(z_segm(1),z_segm(2)-dz_segm(1)/2,Nz_segm(1));
% for si = 2:N_segm-1
%     z_1D(Nz_cumsum(si-1)+1:Nz_cumsum(si-1)+Nz_segm(si))  = linspace(z_segm(si)+dz_segm(si)/2,z_segm(si+1)-dz_segm(si)/2,Nz_segm(si));
% end
% si = N_segm;
% z_1D(Nz_cumsum(si-1)+1:Nz_cumsum(si-1)+Nz_segm(si))  = linspace(z_segm(si)+dz_segm(si)/2,z_segm(si+1),Nz_segm(si));
%
% % smooth grid in z (only in the interior domain)
% for i=1:Nz_smooth, z_1D(Nz_cumsum(1)+1:Nz_cumsum(end-1)) = mvg(z_1D(Nz_cumsum(1)+1:Nz_cumsum(end-1))); end
z_1D            = segmspace(z_segm,Nz_segm,Nz_smooth);

z_sponge1       = z_segm(2);
z_sponge2       = z_segm(end-1);
z1              = z_segm(2);
z2              = z_segm(end-1);

% apply exponential stretching in damping region
Nsponge_top      = sum(r_1D>=r_sponge);
Nsponge_left     = Nz_segm(1);
Nsponge_right    = Nz_segm(end);
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
        z_1D(i)   = z_1D(i-1)+stretch_fact(count)*dz_segm(end);
        count     = count+1;
    end
end
% inlet
if aSponge_left>0
    stretch_fact  = exp(aSponge_left*linspace(0,1,Nsponge_left+1)'); stretch_fact  = stretch_fact(2:end);
    count         = 1;
    for i=Nsponge_left:-1:1
        z_1D(i)   = z_1D(i+1)-stretch_fact(count)*dz_segm(1);
        count     = count+1;
    end
end

% radial derivative
[Dr_1D,D2r_1D,~,~,~,~]      = Dmats_SBP(Nr,dr,FDorder);
[Dr_1D,D2r_1D]              = gridtrans1D(r_1D_equidist,Dr_1D,D2r_1D,r_1D);
% axial derivative
[Dz_1D,D2z_1D,~,~,~,~]      = Dmats_SBP(Nz,dz,FDorder); % fixed to 4th order for now
[Dz_1D,D2z_1D]              = gridtrans1D(z_1D_equidist,Dz_1D,D2z_1D,z_1D);

Dz          = sparse(kron(Dz_1D,speye((Nr))));
Dr          = sparse(kron(speye((Nz)),Dr_1D));
D2z         = sparse(kron(D2z_1D,speye(Nr)));
D2r         = sparse(kron(speye(Nz),D2r_1D));

% build mesh
[z, r]  = meshgrid(z_1D, r_1D);
DR      = kron(speye(5,5),Dr);
DZ      = kron(speye(5,5),Dz);
D2R     = kron(speye(5,5),D2r);
D2Z     = kron(speye(5,5),D2z);

time = toc;
disp(['    elapsed time - Differentiation matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);

tic
% flate plate mean flow
file = [root_data '/' eas_file_mean '.eas'];
% [dim1, dim2, dim3, ndim1, ndim2, ndim3, data, Ma, Re, Pr, T_inf, kappa, ~, ~] = eas2mat( file, ATTRLEN, UDEFLEN);
[ks, data] = eas2mat( file, ATTRLEN, UDEFLEN);
dim1        = [ks.gf(3).dat(1)-100 ks.gf(3).dat ks.gf(3).dat(end)+100];
dim2        = [ks.gf(4).dat ks.gf(4).dat(end)+100];
[dim1,dim2] = meshgrid(dim1,dim2);
for fieldi=1:length(ks.cf)
    switch  strtrim(ks.cf(fieldi,:))
        case('Ma')
            Ma = ks.rf(fieldi);
        case('Re')
            if exist('Re','var')
                disp(['--> USING EFFECTIVE REYNOLDS NUMBER OF Re=' num2str(Re) ' INSTEAD OF KENNSATZ VALUE OF ' num2str(ks.rf(fieldi))])
            else
                Re = ks.rf(fieldi);
            end
        case('kappa')
            kappa = ks.rf(fieldi);
        case('Pr')
            Pr = ks.rf(fieldi);
        case('T_inf')
            T_0 = ks.rf(fieldi);
    end
end
nVars       = 6; % rho,u,v,w,T,p
Q_mean      = zeros(nVars,ks.ndim1+2,ks.ndim2+1);
Q_mean(:,2:end-1,1:end-1)   = squeeze(mean(data(1,1:nVars,:,:,:),5));
Q_mean(:,1,:)       = Q_mean(:,2,:);
Q_mean(:,end,:)     = Q_mean(:,end-1,:);
Q_mean(:,:,end)     = Q_mean(:,:,end-1);


RHO = zeros(Nr,Nz); U = zeros(Nr,Nz); V = zeros(Nr,Nz); W = zeros(Nr,Nz); T = zeros(Nr,Nz); P = zeros(Nr,Nz);

if parallelFlow
    disp('--> PARALLEL FLOW');
    xi    = findnearest(dim1(1,:),x_parallel);
    RHO   = repmat(interp1(dim2(:,1), squeeze(Q_mean(1,xi,:)).',r_1D,'linear'),1,Nz);
    U     = repmat(interp1(dim2(:,1), squeeze(Q_mean(4,xi,:)).',r_1D,'linear'),1,Nz);
    V     = repmat(interp1(dim2(:,1), squeeze(Q_mean(3,xi,:)).',r_1D,'linear'),1,Nz);
    W     = repmat(interp1(dim2(:,1), squeeze(Q_mean(2,xi,:)).',r_1D,'linear'),1,Nz);
    T     = repmat(interp1(dim2(:,1), squeeze(Q_mean(5,xi,:)).',r_1D,'linear'),1,Nz);
else
    RHO   = interp2(dim1, dim2, squeeze(Q_mean(1,:,:)).',z,r,'linear');
    U     = interp2(dim1, dim2, squeeze(Q_mean(4,:,:)).',z,r,'linear');
    V     = interp2(dim1, dim2, squeeze(Q_mean(3,:,:)).',z,r,'linear');
    W     = interp2(dim1, dim2, squeeze(Q_mean(2,:,:)).',z,r,'linear');
    T     = interp2(dim1, dim2, squeeze(Q_mean(5,:,:)).',z,r,'linear');
end
P     = RHO.*T/kappa/Ma^2;%interp2(dim1, dim2, squeeze(Q_mean(6,:,:)).',z,r,'linear');

if strcmp(forcingType,'tke')&&strcmp(caseName,'250ReTheta400')
    file = [root_data '/' eas_file_tke '.eas'];
    [~,~,~,~,~,~,tke,~,~,~,~,~,~,~] = eas2mat( file, ATTRLEN, UDEFLEN);
    TKE                     = zeros(ndim1+2,ndim2+1);
    TKE(2:end-1,1:end-1)    = squeeze(tke);
    TKE                     = interp2(dim1, dim2, TKE',z,r,'linear');
elseif strcmp(forcingType,'tke')&&strcmp(caseName,'ReTheta2200')
    file = [root_data '/' eas_file_tke '.eas'];
    [dim1, dim2,~,~,~,~,data_rms,~,~,~,~,~,~,~] = eas2mat( file, ATTRLEN, UDEFLEN);
    data_rms    = squeeze(data_rms);
    u_rms     = interp2(dim1, dim2, squeeze(data_rms(1,:,:)).',z,r,'linear',0);
    v_rms     = interp2(dim1, dim2, squeeze(data_rms(2,:,:)).',z,r,'linear',0);
    w_rms     = interp2(dim1, dim2, squeeze(data_rms(3,:,:)).',z,r,'linear',0);
    TKE       = 0.5*(u_rms.^2+v_rms.^2+w_rms.^2);
else
    TKE     = ones(Nr,Nz);
end
if parallelFlow
    xi      = findnearest(z_1D,z_segm(end-1));
    TKE     = repmat(TKE(:,xi),1,Nz);
end
TKE     = sparse(TKE);

%
% if strcmp(forcingType,'tke')
%     file = [root_data '/' eas_file_tke '.eas'];
%     [~,~,~,~,~,~,data_rms,~,~,~,~,~,~,~] = eas2mat( file, ATTRLEN, UDEFLEN);
%     data_rms    = squeeze(data_rms);
%     if length(size(data_rms))==2
%         TKE = interp1(dim2, 0.5*squeeze(data_rms(x0_idx,:)).^2,y,'linear',0);
%     else
%         u_rms       = interp1(dim2, squeeze(data_rms(1,x0_idx,:)),y,'linear',0);
%         v_rms       = interp1(dim2, squeeze(data_rms(2,x0_idx,:)),y,'linear',0);
%         w_rms       = interp1(dim2, squeeze(data_rms(3,x0_idx,:)),y,'linear',0);
%         TKE         = 0.5*(u_rms.^2+v_rms.^2+w_rms.^2);
%     end
% else
%     TKE         = zeros(Ny,1);
% end
% TKE     = sparse(TKE);

clear dim1 dim2 dim3 ndim1 ndim2 ndim3 data Q_mean tke

time = toc;
disp(['    elapsed time - Base-flow interpolation: ' datestr(time/24/3600, 'HH:MM:SS')]);

if (vizBaseAndRes)&&(calcCount==1)
    fig_grid = figure('name','Grid');
    pcolor_h=pcolor(z, r, U); set(pcolor_h,'EdgeColor','none'), hold on
    mesh(z,r,zeros(size(z)),'FaceColor','none','edgecolor','k'), text(z(1), r(end), ' $$u$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
    axis equal tight
    
    fig_base = figure('name','Base-state');
    subplot(3,3,1)
    pcolor_h=pcolor(z, r, RHO); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\rho$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,2)
    pcolor_h=pcolor(z, r, W); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$u$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,3)
    pcolor_h=pcolor(z, r, V); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$v$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,4)
    pcolor_h=pcolor(z, r, T); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$T$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
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
U   = 0*reshape(U,     NrNz, 1);
V   = reshape(V,     NrNz, 1);
W   = reshape(W,     NrNz, 1);
T   = reshape(T,     NrNz, 1);
P   = reshape(P,     NrNz, 1);
MU  = reshape(MU,    NrNz, 1);
dMUdT       = reshape(dMUdT, NrNz, 1);
d2MUdT2     = reshape(d2MUdT2, NrNz, 1);
TKE = reshape(TKE,    NrNz, 1);

tic
% Base flow derivatives
dUdr    = 0*Dr*U;     d2Udr2    = 0*D2r*U;
dVdr    = Dr*V;     d2Vdr2    = D2r*V;
dWdr    = Dr*W;     d2Wdr2    = D2r*W;
dRHOdr  = Dr*RHO;   d2RHOdr2  = D2r*RHO;
dTdr    = Dr*T;     d2Tdr2    = D2r*T;
dUdz    = 0*Dz*U;     d2Udz2    = 0*D2z*U;
dVdz    = Dz*V;     d2Vdz2    = D2z*V;
dWdz    = Dz*W;     d2Wdz2    = D2z*W;
dRHOdz  = Dz*RHO;   d2RHOdz2  = D2z*RHO;
dTdz    = Dz*T;     d2Tdz2    = D2z*T;
d2Udrz  = Dz*dUdr;
d2Vdrz  = Dz*dVdr;
d2Wdrz  = Dz*dWdr;
dPdr    = Dr*P;
dPdz    = Dz*P;

% set up sponge region
spongeFun       = zeros(Nr,Nz);
for i = 1:Nr
    for j = 1:Nz
        % left
        if j<=Nsponge_left
            loc     = 1-(z_1D(j)-z_1D(1))/(z_segm(2)-z_1D(1));% 1-(j-1)/(Nsponge_left);
            if strcmp(spongeType,'poly5')
                poly    = -(-6*loc.^5+15*loc.^4-10*loc^3);
            elseif strcmp(spongeType,'pow')
                poly    = loc.^sponge_pow;
            end
            if spongeFun(i,j)<spongeEps_left*poly; spongeFun(i,j) = spongeEps_left*poly; end
        end
        % right
        if j>Nz-Nsponge_right
            %(z_1D(j)-z_segm(end-1))/(z_1D(end)-z_segm(end-1))
            loc     = (z_1D(j)-z_segm(end-1))/(z_1D(end)-z_segm(end-1));% 1+(j-Nz)/Nsponge_right;
            if strcmp(spongeType,'poly5')
                poly    = -(-6*loc.^5+15*loc.^4-10*loc^3);
            elseif strcmp(spongeType,'pow')
                poly    = loc.^sponge_pow;
            end
            if spongeFun(i,j)<spongeEps_right*poly; spongeFun(i,j) = spongeEps_right*poly; end
        end
        % top
        if i>Nr-Nsponge_top
            loc     = (r_1D(i)-r_sponge)/(r_1D(end)-r_sponge);%1+(i-Nr)/Nsponge_top;
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

time    = toc;
disp(['    elapsed time - Base flow derivatives: ' datestr(time/24/3600, 'HH:MM:SS')]);


%%
if (vizBaseAndRes)&&(calcCount==1)
    figure(fig_base)
    subplot(3,3,5)
    pcolor_h=pcolor(z, r, reshape(spongeFun,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' sponge', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,6)
    pcolor_h=pcolor(z, r, reshape(W,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' $$\bar{u}$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,7)
    pcolor_h=pcolor(z, r, reshape(dWdr,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$\bar{u}$$/d$$r$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,8)
    pcolor_h=pcolor(z, r, reshape(dWdz,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$\bar{u}$$/d$$z$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    subplot(3,3,9)
    pcolor_h=pcolor(z, r, reshape(d2Wdrz,Nr,Nz)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ' d$$^2\bar{u}$$/d$$r$$d$$z$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
    set(gcf,'Position',[171         103        1324         666])
    
    drawnow
end

%%
tic

% Coefficient matrices
getCoeffsLaminar_Cartesian

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

% farfield
rho_farf_i  = ti_rho;
u_farf_i  = ti_u;
v_farf_i  = ti_v;
w_farf_i  = ti_w;
T_farf_i  = ti_T;
% bottom
rho_bi_i  = bi_rho;
u_bi_i  = bi_u;
v_bi_i  = bi_v;
w_bi_i  = bi_w;
T_bi_i  = bi_T;
% inlet
rho_in_i    = li_rho;
u_in_i    = li_u;
v_in_i    = li_v;
w_in_i    = li_w;
T_in_i    = li_T;
% outlet
rho_out_i = ri_rho;
u_out_i   = ri_u;
v_out_i   = ri_v;
w_out_i   = ri_w;
T_out_i   = ri_T;

size_LHS = [5*NrNz 5*NrNz];

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
    DR_i    = zeros(FDorder*length(index_set),1); % 5-point stencil with 4 entries
    DR_j    = DR_i;
    DR_val  = DR_i;
    for row_i=1:length(index_set)
        DR_i((row_i-1)*FDorder+1:row_i*FDorder)    = ones(FDorder,1)*index_set(row_i);%[index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
        [~,jj,val] = find(DR(index_set(row_i),:));
        DR_j((row_i-1)*FDorder+1:row_i*FDorder)    = jj;
        DR_val((row_i-1)*FDorder+1:row_i*FDorder)  = val;
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
    Diag(diag_index_set)= 1;
    LHS = LHS + Diag;
    RHS = RHS + Diag;
end

% inlet
if strcmpi(inletBC,'zero_gradient')
    index_set = [u_in_i v_in_i w_in_i rho_in_i T_in_i];
    LHS(index_set, :)   = 0;
    RHS(index_set, :)   = 0;
    diag_index_set      = sub2ind(size_LHS, index_set, index_set);
    Diag                = speye(5*NrNz)*0;
    Diag(diag_index_set)= 1;
    LHS = LHS + Diag;
    
    [RHS_i,RHS_j,RHS_val] = find(RHS);
    DZ_i    = zeros(FDorder*length(index_set),1); % 5-point stencil with 4 entries
    DZ_j    = DZ_i;
    DZ_val  = DZ_i;
    for row_i=1:length(index_set)
        DZ_i((row_i-1)*FDorder+1:row_i*FDorder)    = ones(FDorder,1)*index_set(row_i);%[index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
        [~,jj,val] = find(DZ(index_set(row_i),:));
        DZ_j((row_i-1)*FDorder+1:row_i*FDorder)    = jj;
        DZ_val((row_i-1)*FDorder+1:row_i*FDorder)  = val;
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
    Diag(diag_index_set)= 1;
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
    DZ_i    = zeros(FDorder*length(index_set),1); % 5-point stencil with 4 entries
    DZ_j    = DZ_i;
    DZ_val  = DZ_i;
    for row_i=1:length(index_set)
        DZ_i((row_i-1)*FDorder+1:row_i*FDorder)    = ones(FDorder,1)*index_set(row_i);%[index_set(row_i); index_set(row_i); index_set(row_i); index_set(row_i)];
        [~,jj,val] = find(DZ(index_set(row_i),:));
        DZ_j((row_i-1)*FDorder+1:row_i*FDorder)    = jj;
        DZ_val((row_i-1)*FDorder+1:row_i*FDorder)  = val;
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
    Diag(diag_index_set)= 1;
    LHS = LHS + Diag;
    RHS = RHS + Diag;
end

% wall
if strcmpi(bottomBC,'wall')
    pseudo_EV = -1000-1000i;
    index_set = [bi_u bi_v bi_w bi_T];
    % Dirichlet BCs but for rho @ wall (Bjoern)
    LHS(index_set,:)    = 0;
    diag_index_set      = sub2ind(size_LHS, index_set, index_set);
    Diag                = speye(5*NrNz)*0;
    Diag(diag_index_set)= pseudo_EV;
    LHS                 = LHS + Diag;
else
    error('bottomBC other than ''wall'' not implemented')
end

time    = toc;
disp(['    elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

clear D2* Dr* Dz* DR DZ Diag DB* RHS_* DZ_* DR_*

%%%%%%%%%%%%%
% Solve EVP %
%%%%%%%%%%%%%

% undo Fourier ansatz in time if dq/dt = L0*q is needed, else keep (L0-omega*I)*q = 0 form as EVP
if strcmpi(calcMethod,'timestepper')||strcmpi(calcMethod,'iterate_optimal')||strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'les_forcingFFT_response')||strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'pseudospectrum')||strcmpi(calcMethod,'calculate_forcings')
    RHS = RHS/1i;
end

% reduce to generalized to standard EVP
L0               = RHS\LHS;
clear LHS RHS
%%
for SIGMA =  reshape(SIGMA_vec,1,size(SIGMA_vec,1)*size(SIGMA_vec,2)) % makes sure that SIGMA_vec can be a matrix (pseudospectrum)
    % set up save file name
    [saveFile] = getSaveFileName_turbBL(saveFolder,custom_comment,calcMethod,Re,kz,Nr,Nz,r_farf,z1,z2,SIGMA);
    % check if case is already calculating or already exists, and skip if so
    if (~exist(saveFile,'file')&&~exist([saveFile '.lock'],'file'))||(~checkLock)||strcmpi(calcMethod,'pseudospectrum')
        if (checkLock)&&~strcmpi(calcMethod,'pseudospectrum'), save([saveFile '.lock'],'SIGMA'); end
        
  %%      
        if strcmpi(calcMethod,'les_forcingPOD_response')
            fi = findnearest(freq_SPOD,SIGMA/2/pi);
            file_load  = [tfft_root '/SPOD_freq' num2str(fi,'%.4i')];
            if ~isempty(comment_SPOD), file_load = [file_load '_' comment_SPOD]; end
            file_load  = [file_load '.mat'];
            load(file_load,'q_SPOD');
            noEigs      = size(q_SPOD,1);
            kzi         = findnearest(kz_SPOD,kz);
            q_SPOD      = squeeze(q_SPOD(:,:,:,:,kzi));
            %u       = squeeze(real(q_SPOD(modei,2,:,:,kzi)));            
            U_out       = zeros(5*NrNz,noEigs);
            %[r_SPOD,x_SPOD] = meshgrid(r_SPOD,x_SPOD);
            for eig_i = 1:noEigs
                U_out(0*NrNz+1:1*NrNz,eig_i) = reshape(interp2(x_SPOD,r_SPOD,squeeze(q_SPOD(eig_i,1,:,:)).',z,r,'spline',0),NrNz,1);
                U_out(1*NrNz+1:2*NrNz,eig_i) = reshape(interp2(x_SPOD,r_SPOD,squeeze(q_SPOD(eig_i,4,:,:)).',z,r,'spline',0),NrNz,1);
                U_out(2*NrNz+1:3*NrNz,eig_i) = reshape(interp2(x_SPOD,r_SPOD,squeeze(q_SPOD(eig_i,3,:,:)).',z,r,'spline',0),NrNz,1);
                U_out(3*NrNz+1:4*NrNz,eig_i) = reshape(interp2(x_SPOD,r_SPOD,squeeze(q_SPOD(eig_i,2,:,:)).',z,r,'spline',0),NrNz,1);
                U_out(4*NrNz+1:5*NrNz,eig_i) = reshape(interp2(x_SPOD,r_SPOD,squeeze(q_SPOD(eig_i,5,:,:)).',z,r,'spline',0),NrNz,1);
            end
            
            if (vizBaseAndRes)&&(calcCount==1)
                figure('name','response SPOD modes')
                if noEigs<5, n_row = noEigs; else n_row = 5; end
                n_col = ceil(noEigs/5);
                w_q_all     = zeros(Nr, Nz);
                for i = 1:noEigs
                    subplot(n_row,n_col,i)                
                    pcolor_h = pcolor(z, r, real(reshape(U_out(3*NrNz+1:4*NrNz,i),Nr,Nz)));
                    set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1]); shading interp; colormap(jetblack)
                    xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge])             
                end
%                 subplot(2,3,1)
%                 pcolor_h=pcolor(z, r, real(reshape(U_out(0*NrNz+1:1*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['${\rho}''}^{(' num2str(1) ')} (kz='  num2str(kz) ', f=' num2str(SIGMA/2/pi,'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge]), colormap(jetblack)
%                 subplot(2,3,2)
%                 pcolor_h=pcolor(z, r, real(reshape(U_out(1*NrNz+1:2*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['${\hat{v}_{z}''}^{(' num2str(1) ')} (kz='  num2str(kz) ', f=' num2str(SIGMA/2/pi,'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge]), colormap(jetblack)
%                 subplot(2,3,3)
%                 pcolor_h=pcolor(z, r, real(reshape(U_out(2*NrNz+1:3*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['${\hat{v}_{y}''}^{(' num2str(1) ')} (kz='  num2str(kz) ', f=' num2str(SIGMA/2/pi,'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge]), colormap(jetblack)
%                 subplot(2,3,4)
%                 pcolor_h=pcolor(z, r, real(reshape(U_out(3*NrNz+1:4*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['${\hat{v}_{x}''}^{(' num2str(1) ')} (kz='  num2str(kz) ', f=' num2str(SIGMA/2/pi,'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge]), colormap(jetblack)
%                 subplot(2,3,5)
%                 pcolor_h=pcolor(z, r, real(reshape(U_out(4*NrNz+1:5*NrNz,1),Nr,Nz))); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(z(1), r(end), ['${\hat{T}''}^{(' num2str(1) ')} (kz='  num2str(kz) ', f=' num2str(SIGMA/2/pi,'%.3f') ')$'], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge]), colormap(jetblack)
                drawnow
            end
        end
%%        
        if (strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'les_forcingPOD_response'))
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Input-output analysis %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % see Nichols, CTR report 2014 "Input-output analysis of high-speed jet noise"
            disp(['    calling intput-output analysis with SIGMA=' num2str(SIGMA)]);
            tic
            
            % weights induced by compressible energy norm
            intWeight   = reshape(trapzWeights2(r_1D,z_1D),NrNz,1);
            vol         = 1;
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
            
            %%%%%%%%%%%
            % FORCING %
            %%%%%%%%%%%
            switch forcingType
                case {'global','tke'}
                    B_win       = ones(Nr,Nz);
                    B_win(z<INPUT_z_min|z>INPUT_z_max|r<INPUT_r_min|r>INPUT_r_max) = 0;
                case {'strip'}
                    %%
                    strip_fun   = zeros(1,Nz);
                    L_strip     = INPUT_z_max-INPUT_z_min;
                    for j=1:Nz
                        if z_1D(j)>=INPUT_z_min&&z_1D(j)<=INPUT_z_min+L_strip/2
                            % left
                            loc             = (z_1D(j)-INPUT_z_min)/(L_strip/2);
                            poly            = -(-6*loc.^5+15*loc.^4-10*loc^3);
                            strip_fun(j)    = poly;
                        elseif z_1D(j)>=INPUT_z_min+L_strip/2&&z_1D(j)<=INPUT_z_max
                            % right
                            loc             = (INPUT_z_max-z_1D(j))/(L_strip/2);
                            poly            = -(-6*loc.^5+15*loc.^4-10*loc^3);
                            strip_fun(j)    = poly;
                        else
                            strip_fun(j)    = 0;
                        end
                    end
                    B_win                               = repmat(strip_fun,Nr,1);
                    B_win(r<INPUT_r_min|r>INPUT_r_max)  = 0;                    
                otherwise
                    error('Unexpected forcingType!')
            end
            % smoothing
            B_win       = B_win+1;
            for si = 1:nSmooth_B_x
                for ri = 1:Nr
                    B_win(ri,:)     = mvg(B_win(ri,:));
                end
            end
            for si = 1:nSmooth_B_y
                for zi = 1:Nz
                    B_win(:,zi)     = mvg(B_win(:,zi));
                end
            end
            B_win       = B_win-1;
            
            if vizBaseAndRes, figure('name','forcing mask')
                pcolor(z,r,B_win), shading interp, axis equal tight
            end
            
            B_win       = reshape(sparse(B_win),NrNz,1);
            Z = 0*speye(NrNz);
            switch forcingType
                case {'global'}
                    B = [ ...
                        Z                Z                Z                Z                            Z
                        Z                diag(B_win).*speye(NrNz)          Z                Z           Z
                        Z                Z                diag(B_win).*speye(NrNz)          Z           Z
                        Z                Z                Z                diag(B_win).*speye(NrNz)     Z
                        Z                Z                Z                Z                            Z
                        ];
                case {'tke'}
                    B = [ ...
                        Z                Z                Z                Z                            Z
                        Z                diag(B_win).*diag(TKE)            Z                Z           Z
                        Z                Z                diag(B_win).*diag(TKE)            Z           Z
                        Z                Z                Z                diag(B_win).*diag(TKE)       Z
                        Z                Z                Z                Z                            Z
                        ];
                case   {'strip'}
                     B = [ ...
                        Z                Z                Z                Z                            Z
                        Z                diag(B_win).*speye(NrNz)          Z                Z           Z
                        Z                Z                diag(B_win).*speye(NrNz)          Z           Z
                        Z                Z                Z                diag(B_win).*speye(NrNz)     Z
                        Z                Z                Z                Z                            Z
                        ];                           
                otherwise
                    error('Unexpected forcingType!')
            end
            
            %%%%%%%%%%%%
            % RESPONSE %
            %%%%%%%%%%%%           
            C_win       = ones(Nr,Nz);
            C_win(z<OUTPUT_z_min|z>OUTPUT_z_max|r<OUTPUT_r_min|r>OUTPUT_r_max) = 0;
            C_win       = C_win+1;
            for si = 1:nSmooth_C_x
                for ri = 1:Nr
                    C_win(ri,:)     = mvg(C_win(ri,:));
                end
            end
            for si = 1:nSmooth_C_y
                for zi = 1:Nz
                    C_win(:,zi)     = mvg(C_win(:,zi));
                end
            end
            C_win       = C_win-1;
            C_win       = reshape(sparse(C_win),NrNz,1);
            C = [ ...
                0*diag(C_win).*speye(NrNz)                Z                Z                Z                            Z
                Z                diag(C_win).*speye(NrNz)          Z                Z           Z
                Z                Z                diag(C_win).*speye(NrNz)          Z           Z
                Z                Z                Z                diag(C_win).*speye(NrNz)     Z
                Z                Z                Z                Z                            0*diag(C_win).*speye(NrNz)
                ];
            
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
                FZ  = sparse([ ...parallelparallel
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
                %Htr     = @(v) invF*(   FZ_ct*(MZ_ct\(FR_ct*(MR_ct\(   B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(v))))))))))));% CM^-1BS (f only)
                %H       = @(v) F*(C*(qq*(UU\(LL\(pp*(rr\(B*    (MR\(FR*(MZ\(FZ*    v)))))))))));% CM^-1BS (f only)
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
            if (strcmpi(calcMethod,'les_forcingPOD_response'))
                V_in = zeros(5*Nz*Nr,noEigs);
                for eig_i = 1:noEigs
                    % normalize given output to get gain right
                    U_out(:,eig_i)  = U_out(:,eig_i)/sqrt(U_out(:,eig_i)'*F*U_out(:,eig_i));
                    % get optimal forcing for given output
                    V_in(:,eig_i)   = LsI*U_out(:,eig_i);
                end
            else
                [V_in, OMEGA]       = eigs(HtrH, 5*Nz*Nr, noEigs, 'lm', opts);
                gain                = sqrt(real(diag(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))
                U_out               = zeros(size(V_in));
            end
            
            time = toc;
            disp(['    elapsed time - SVD via ''eigs'' using matrix-free Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);            %%
            % output
            for eig_i = 1:noEigs
                % normalization in inner product norm
                V_in(:,eig_i)  = V_in(:,eig_i)/sqrt(V_in(:,eig_i)'*F*V_in(:,eig_i));
                if strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response')
                    gain(eig_i)     = real(sqrt(V_in(:,eig_i)'*(B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(F*(C*(qq*(UU\(LL\(pp*(rr\(B*V_in(:,eig_i)))))))))))))))))); % real() because of O(eps) imaginary error
                end
                % % without filter
                if filter_q
                    U_out(:,eig_i) = MR\(FR*(MZ\(FZ*(C*(qq*(UU\(LL\(pp*(rr\(B*(V_in(:,eig_i)/gain(eig_i))))))))))));
                else
                    U_out(:,eig_i) = C*(qq*(UU\(LL\(pp*(rr\(B*(V_in(:,eig_i)/gain(eig_i))))))));
                end
            end
            
            u_max_x     = zeros(Nz,noEigs);
            u_max_y     = zeros(Nr,noEigs);
            u_int_y     = zeros(Nr,noEigs);
            x_maxu      = zeros(1,noEigs);
            for iii=1:noEigs
                [~,idx]         = max(abs(U_out(w_j,iii)));
                x_maxu(iii)     = z(idx);
                tmp             = reshape(U_out(w_j,iii),Nr,Nz);
                u_max_y(:,iii)  = tmp(:,z_1D==x_maxu(iii));
                u_int_y(:,iii)  = sum(tmp,2);
                u_max_x(:,iii)  = max(tmp,[],1);
            end
            clear tmp yWeight
            %figure, subplot(1,2,1), semilogx(r_1D,abs(u_max_y).*gain'), subplot(1,2,2), plot(z_1D,abs(u_max_x).*gain')
            
            % % check orthogonality
            % abs(V_in'*V_in)         % I/O should be orthonormal in 2-norm
            % abs(V_in'*F*V_in)       % frequency response should be orthonormal in (induced) energy norm
            % abs(U_out'*U_out)
            % abs(U_out'*F*U_out)
            
            clear HtrH LsI LL  UU pp qq rr C B LL_ct UU_ct pp_ct qq_ct rr_ct C_ct B_ct opts B_win C_win B C M F invF MR_ct MR FR_ct FR MZ_ct MZ FZ_ct FZ
            clear diag_index_set dMUdT dPdr dPdz index_set intWeight mu MU P p_f p_q spongeFun T_bi_i t_f T_farf_i  ...
                T_in_i T_inf T_out_i t_q T_ref ti ti_rho ti_T ti_u ti_v ti_w time Tstar u_bi_i u_farf_i u_in_i ...
                u_out_i v_bi_i v_farf_i v_in_i v_out_i vol w_bi_i w_f w_farf_i w_in_i w_out_i w_q w_q_all ...
                bi bi_rho bi_T bi_u bi_v bi_w li li_rho li_T li_u li_v li_w ri_w ri_v ri_u ri_T ri_rho ri rho_out_i rho_in_i rho_farf_i rho_f rho_bi_i
            
            %%
            if (vizBaseAndRes)
                %%
                figure('name','gain factors');
                stem(gain);
                xlabel('mode number N'); ylabel('singular value $$\sigma$$');
                
                %                 fig_f_p   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - forcing pressure'],'fileName',['resolvent_f_pressure_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                fig_f_w   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - forcing u_x'],'fileName',['resolvent_f_u_x_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                %                 fig_q_p   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - response pressure'],'fileName',['resolvent_q_pressure_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                fig_q_w   = figure('name',['St=' num2str(SIGMA/2/pi,'%.2g') ' - response u_x'],'fileName',['resolvent_q_u_x_St' num2str(SIGMA/2/pi,'%.2g') '.fig']);
                
                if noEigs<5, n_row = noEigs; else n_row = 5; end
                n_col = ceil(noEigs/5);
                w_q_all     = zeros(Nr, Nz);
                for i = 1:noEigs
                    rho_f   = reshape(V_in(rho_j,i), Nr, Nz);
                    w_f     = reshape(V_in(w_j,i), Nr, Nz);
                    t_f     = reshape(V_in(T_j,i), Nr, Nz);
                    p_f     = (reshape(RHO, Nr, Nz).*t_f + rho_f.*reshape(T, Nr, Nz))/kappa/Ma^2;
                    rho_q   = reshape(U_out(rho_j,i), Nr, Nz);
                    w_q     = reshape(U_out(w_j,i), Nr, Nz);
                    t_q     = reshape(U_out(T_j,i), Nr, Nz);
                    p_q     = (reshape(RHO, Nr, Nz).*t_q + rho_q.*reshape(T, Nr, Nz))/kappa/Ma^2;
                    w_q_all = w_q_all + gain(i)*abs(w_q);
                    
                    %                     figure(fig_f_p)
                    %                     subplot(n_row,n_col,i)
                    %                     pcolor_h=pcolor(z,r,real(p_f)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    %                     caxis(0.5*max(abs(caxis))*[-1 1]); shading interp; colormap(jetblack)
                    %                     hold on
                    %                     ylabel('$r$');
                    
                    figure(fig_f_w)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(w_f)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1]); shading interp;
                    xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge])
                    hold on
                    ylabel('$r$');
                    
                    %                     figure(fig_q_p)
                    %                     subplot(n_row,n_col,i)
                    %                     pcolor_h=pcolor(z,r,real(p_q)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    %                     caxis(0.5*max(abs(caxis))*[-1 1]); shading interp; colormap(jetblack)
                    %                     hold on
                    %                     ylabel('$r$');
                    
                    figure(fig_q_w)
                    subplot(n_row,n_col,i)
                    pcolor_h=pcolor(z,r,real(w_q)); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
                    caxis(0.5*max(abs(caxis))*[-1 1]); shading interp;
                    xlim([z_segm(2) z_segm(end-1)]); ylim([0 r_sponge])
                    hold on
                    ylabel('$r$');
                end
%                 figure(fig_q_w)
%                 subplot(n_row+1,n_col,1)
%                 pcolor_h=pcolor(z,r,w_q_all.^2); set(pcolor_h,'EdgeColor','none'); axis equal, axis tight
%                 caxis(0.5*max(abs(caxis))*[-1 1]); shading interp; colormap(jetblack)
                
                this_fig = {fig_f_w,fig_q_w};%{fig_f_p,fig_f_w,fig_q_p,fig_q_w};
                for figi = 1:2%4
                    figure(this_fig{figi})
                    xlabel('$x$');
                    fig_pos = get(gcf,'position'); set(gcf,'Position',[fig_pos(1) fig_pos(2) n_col*fig_pos(3) fig_pos(3)*n_row*((1.3*r_sponge)/z_sponge2)]);
                    fig_w   = fig_pos(3); fig_h = fig_pos(4);
                end
                %%clear fig_f_p fig_f_w fig_q_p fig_q_w rho_f w_f t_f p_f rho_q w_q t_q p_q
            end
            time = toc;
            disp(['    elapsed time - intput-output analysis: ' datestr(time/24/3600, 'HH:MM:SS')]);
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
                %save(saveFile, '-regexp', '^(?!(L0)$).', '-v7.3');
                V_in    = single(V_in);
                U_out   = single(U_out);
                save(saveFile,'gain','V_in','U_out', '-v7.3');
            else
                save(saveFile, '-regexp', '^(?!(L0|eigVecs|V_in|U_out)$).', '-v7.3');
            end
            saveFile_info   = getSaveFileName_turbBL(saveFolder,[custom_comment '.info'],calcMethod,Re,0,Nr,Nz,r_farf,z1,z2,0);
            vari_count      = 1;            
            if ~exist(saveFile_info,'file')             
                for vari = 1:length(varList)
                    var     = char(varList(vari));
                    if exist(var,'var')
                        if vari_count>1
                            save(saveFile_info,var,'-append');
                        else
                            save(saveFile_info,var,'-v7.3');
                        end
                        vari_count = vari_count+1;
                    end
                end
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
% %%
% if (vizBaseAndRes)
%     figure('name','spectrum');
%     plot(real(OMEGA),imag(OMEGA),'go')
%     hold on
%     for i=1:noEigs
%         text(real(OMEGA(i)),imag(OMEGA(i)),[' ' num2str(i)])
%     end
%     plot(real(SIGMA),imag(SIGMA),'rx');
%     
%     nrows = floor(sqrt(noEigs));
%     ncols = ceil(noEigs/nrows);
%     fig_overview_w = figure('name', 'w''');
%     fig_overview_p = figure('name', 'p''');
%     for i=1:noEigs
%         figure(fig_overview_w)
%         subplot(nrows,ncols,i)
%         w = real(reshape(eigVecs(w_j,i), Nr, Nz));
%         pcolor_h=pcolor(z(1,:), r(:,1), w); set(pcolor_h,'EdgeColor','none'), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[]), axis equal, axis tight
%         colormap(jetblack)
%         text(z(1), r(end), ['mode ' num2str(i)], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
%         sub_pos = get(gca,'position'); % get subplot axis position
%         set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
%         [val idx] = max(abs(w(:)));
%         caxis([-max(abs(w(idx))) max(abs(w(idx)))]);
%         hold on
%         figure(fig_overview_p)
%         subplot(nrows,ncols,i)
%         rho = reshape(eigVecs(rho_j,i), Nr, Nz);
%         t   = reshape(eigVecs(T_j,i), Nr, Nz);
%         p = (reshape(RHO, Nr, Nz).*t + rho.*reshape(T, Nr, Nz))/kappa/Ma^2;
%         pcolor_h=pcolor(z(1,:), r(:,1), real(p)); set(pcolor_h,'EdgeColor','none'), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[]), axis equal, axis tight
%         colormap(jetblack)
%         text(z(1), r(end), ['mode ' num2str(i)], 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left')
%         sub_pos = get(gca,'position'); % get subplot axis position
%         set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
%         [val idx] = max(abs(p(:)));
%         caxis([-max(abs(p(idx))) max(abs(p(idx)))]);
%         hold on
%     end
%     
%     while 1
%         i = input('Eigenvector # (or 0 to quit)?');
%         if i==0, break, end
%         
%         % detailed individual mode plot
%         figure('name',['mode#' num2str(i)],'position',[2         452        1346         508])
%         rho = reshape(eigVecs(rho_j,i), Nr, Nz);
%         w = reshape(eigVecs(w_j,i), Nr, Nz);
%         u = reshape(eigVecs(u_j,i), Nr, Nz);
%         v = reshape(eigVecs(v_j,i), Nr, Nz);
%         t = reshape(eigVecs(T_j,i), Nr, Nz);
%         p = (reshape(RHO, Nr, Nz).*t + rho.*reshape(T, Nr, Nz))/kappa/Ma^2;
%         sponge = reshape(spongeFun, Nr, Nz);
%         
%         [~,idx]       = max(abs(w(:)));
%         [z_idx,r_idx]   = ind2sub(size(w),idx);
%         % contourplot streamwise plane
%         subplot(2,4,1:3)
%         plotvar = real(p);
%         pcolor_h=pcolor(z, r, plotvar); set(pcolor_h,'EdgeColor','none'), axis tight, hold on
%         plot([z(1,r_idx) z(1,r_idx)],[r(1) r(end)],'r:')
%         plot([z(1) z(end)],[r(z_idx,1) r(z_idx,1)],'r:')
%         contour(z, r, sponge, [10*eps 10*eps], 'k:');
%         caxis([-0.1*(norm(plotvar(:),Inf)) 0.1*(norm(plotvar(:),Inf))]);
%         colormap(jetblack)
%         c = colorbar('east');
%         c.Label.String = 'real(p'')';
%         xlabel('x'); ylabel('r');
%         hold on
%         [c1,h1]=contour(z, r, real(w), linspace(-norm(w(:),Inf),norm(w(:),Inf),10), 'edgecolor','r');
%         contour_negDashedposSolid(c1,h1)
%         % radial profiles
%         subplot(2,4,[4 8])
%         plot(abs(  w(:,r_idx)), r_1D, 'r-'), hold on
%         plot(abs(  u(:,r_idx)), r_1D, 'k:')
%         plot(abs(rho(:,r_idx)), r_1D, 'k--')
%         plot(abs(  t(:,r_idx)), r_1D, 'k-.')
%         plot(abs(  p(:,r_idx)), r_1D, 'b-')
%         % plot(abs(v(:,r_idx)), r_1D, 'k.')
%         legend('|u''|','|v''|','|\rho''|','|T''|','|p''|')
%         xlabel('r'), ylabel('amplitude'), axis tight
%         % streamwise profiles
%         subplot(2,4,5:7)
%         plot(z_1D, abs(  w(z_idx,:)), 'r-'), hold on
%         plot(z_1D, abs(  real(w(z_idx,:))), 'k-', 'color', [.5 .5 .5])
%         plot(z_1D, abs(  u(z_idx,:)), 'k:')
%         plot(z_1D, abs(rho(z_idx,:)), 'k--')
%         plot(z_1D, abs(  t(z_idx,:)), 'k-.')
%         plot(z_1D, abs(  p(z_idx,:)), 'b-')
%         legend('|u''|','real(u'')','|v''|','|\rho''|','|T''|','|p''|','location','northeast')
%         xlabel('x'), ylabel('amplitude'), axis tight
%         
%         % pressure and velocity contour plots
%         figure('name',['mode#' num2str(i)],'position',[575           2        1345         453])
%         subplot(2,1,1)
%         plotvar = real(p(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right)); varname = 'p''';
%         plotvar = plotvar./norm(plotvar(:),Inf);
%         pcolor_h=pcolor(z(1,Nsponge_left+1:end-Nsponge_right), r(1:end-Nsponge_top,1), plotvar); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on
%         caxis([-0.5 0.5]);
%         colormap(jetblack)
%         c = colorbar('eastoutside');
%         c.Label.String = varname;
%         xlabel('x'); ylabel('r'); title('pressure');
%         subplot(2,1,2)
%         plotvar = real(w(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right)); varname = 'u_x''';
%         plotvar = plotvar./norm(plotvar(:),Inf);
%         pcolor_h=pcolor(z(1,Nsponge_left+1:end-Nsponge_right), r(1:end-Nsponge_top,1), plotvar); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on
%         caxis([-0.5 0.5]);
%         colormap(jetblack)
%         c = colorbar('eastoutside');
%         c.Label.String = varname;
%         xlabel('x'); ylabel('r'); title('streamwise velocity');
%         drawnow
%     end
% end

% profile viewer
% p = profile('info');














