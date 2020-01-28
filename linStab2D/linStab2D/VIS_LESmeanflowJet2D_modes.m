clear all
% close all
clc
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
addpath('../aux');

useTanhTrafo = 0;

% % M=0.9 & M=0.4
% r_farf          = 7.0923;
% z_end           = 31.0087;
% custom_comment  = 'M04_global'
% % custom_comment  = 'M04_global_ReStudyConstSt'
% % custom_comment  = 'M04_global_leftForcing10'
% % custom_comment  = 'M04_global_rightForcing10'

% LESmeanflowJet2D_frequency_response_Re=30000_m=1_Nr=195_r=7.1804_Nz=950_-1>z>31.0677_ReSIGMA=2.2907_ImSIGMA=0.mat
% B118 & B122
r_farf          = 7.1804;
z_end           = 31.0677;
custom_comment  = ''%'B118_global'

calcMethod      = 'frequency_response';
Re              = 30000%3727.5937%30000;
Nr              = 195;
Nz              = 950;
z_start         = -1;
m               = 0;
St_cmp          = 0.6;

% % M=0.4 LES grid calcs for SPOD paper
% r_farf          = 7.0739;
% z_end           = 31.0087;
% custom_comment  = 'M04_global_LESgrid'
% 
% calcMethod      = 'frequency_response';
% Re              = 10000;
% Nr              = 143;
% Nz              = 869;
% z_start         = -1;
% m               = 0;
% St_cmp          = 1;

% get file
files   = getSaveFileName('results',custom_comment,calcMethod,Re,m,Nr,Nz,r_farf,z_start,z_end,'*');
files   = dir(files);
St    = zeros(length(files),1);
for ii=1:length(files)
    load(['results/' files(ii).name],'SIGMA');
    St(ii)      = SIGMA/2/pi;
end
St_i    = findnearest(St,St_cmp);
load(['results/' files(St_i).name]);

% set some addition visualization parameters
r_line_pos      = 0.5;  % z position of line along r or 'auto' for max. streamwise velocity position    
z_line_pos      = 0.0;  % r position of line along z or 'auto' for max. streamwise velocity position    
filter_vis      = true;
interp_vis      = true;
z_vis_0         = 0;
z_vis_1         = z_sponge2;
n_vis_z         = 1000;
r_vis_0         = 0;
r_vis_1         = r_sponge;
n_vis_r         = ceil((r_vis_1-r_vis_0)/((z_sponge2-z_sponge1)/n_vis_z))
render3D        = false;
Nr_interp       = 200; % even to remove singularity
Nphi_interp     = 180;

%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%

if strcmpi(calcMethod,'input_output')||strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response')
    figure('name','gain factors');
    stem(gain);
    xlabel('mode number N'); ylabel('singular value $$\sigma$$');
else
    figure('name','spectrum');
    plot(real(OMEGA)/2/pi,imag(OMEGA),'go')
    hold on
    for i=1:noEigs
        text(real(OMEGA(i))/2/pi,imag(OMEGA(i)),[' ' num2str(i)])
    end
    plot(real(SIGMA)/2/pi,imag(SIGMA),'rx');
end

i = 1;
while 1    
    if noEigs>1, i = input('mode # (ENTER to quit)? '); else, i = 1; end
    if isempty(i), break, end
    
    if strcmpi(calcMethod,'input_output')||strcmpi(calcMethod,'pseudospectrum')||strcmpi(calcMethod,'frequency_response')||strcmpi(calcMethod,'les_forcingPOD_response')||strcmpi(calcMethod,'les_forcingFFT_response')
        io = input('input (1) or output (2)? ');
        if io==1
            eigVecs = V_in;
        else
            eigVecs = U_out;
        end
    end
    
    w               = reshape(eigVecs(w_j,i), Nr, Nz);
    [val,idx]       = max(abs(w(:)));
    if strncmp(r_line_pos,'auto',4)
        [z_idx,~]   = ind2sub(size(w),idx);
    else
        z_idx = findnearest(r_line_pos,r(:,1));
    end
    
    if strncmp(z_line_pos,'auto',4)
        [~,r_idx]   = ind2sub(size(w),idx);
    else
        r_idx = findnearest(z_line_pos,z(1,:));
    end
    
    % nice plot
    figure('name',['Isocontours of streamwise velocity for mode ' num2str(i) ', St=' num2str(real(SIGMA)/2/pi)],'position',[2         452        1346         508])
    rho = reshape(eigVecs(rho_j,i), Nr, Nz);
    w = reshape(eigVecs(w_j,i), Nr, Nz);
    u = reshape(eigVecs(u_j,i), Nr, Nz);
    v = reshape(eigVecs(v_j,i), Nr, Nz);
    t = reshape(eigVecs(T_j,i), Nr, Nz);
    p = (reshape(RHO, Nr, Nz).*t + rho.*reshape(T, Nr, Nz))/kappa/Ma^2;
    %dilatation  = reshape(Dr*(eigVecs(rho_j,i).*eigVecs(w_j,i))+Dz*(eigVecs(rho_j,i).*eigVecs(u_j,i)), Nr, Nz);
    %vort_theta  = reshape(Dz*eigVecs(u_j,i)-Dr*eigVecs(w_j,i), Nr, Nz);
    sponge = reshape(spongeFun, Nr, Nz);
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
    c.Label.String = '$$Re(p'')$$'; c.TickLabelInterpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    xlabel('$x$'); ylabel('$$r$$');
    hold on
    isoVals = linspace(0,norm(w(:),Inf),5); isoVals = isoVals(2:end-1);
    contour(z, r, real(w), isoVals,  'edgecolor','r');
    contour(z, r, real(w), -isoVals, 'edgecolor','b');
    
    % find wavepackets 
    f = mvg(mvg(mvg(filter_O10(abs(p(z_idx,:)),0.4)))); % use smoothed pressure -> find valleys in abs 
    [pks,locs] = findpeaks(1./f,z_1D);
    for pi = 1:length(locs)
       plot([locs(pi) locs(pi)],[r(1) r(end)],'g--') 
    end
%     if pipeBC, rectangle('position',[z(pipe_idx_bottom(1)),r(pipe_idx_bottom(1)),z(pipe_idx_bottom(end))-z(pipe_idx_bottom(1)),r(pipe_idx_top(1))-r(pipe_idx_bottom(1))],'edgecolor','b','facecolor','b'), end
    % radial profiles
    subplot(2,4,[4 8])
    plot(abs(  w(:,r_idx)), r_1D, 'r-'), hold on
    plot(abs(  u(:,r_idx)), r_1D, 'k-')
    plot(abs(rho(:,r_idx)), r_1D, 'k--')
    plot(abs(  t(:,r_idx)), r_1D, 'k-.')
    plot(abs(  p(:,r_idx)), r_1D, 'b-')
    % plot(abs(v(:,r_idx)), r_1D, 'k:')
    legend('$$|u''_{x}|$$','$$|u''_{r}|$$','$$|\rho''|$$','$$|T''|$$','$$|p''|$$')
    ylabel('$$r$$'), xlabel('$$|\cdot|$$'), axis tight
    % streamwise profiles
    subplot(2,4,5:7)
    plot(z_1D, abs(  w(z_idx,:))/max(abs(w(z_idx,:))), 'r-'), hold on
    plot(z_1D, abs(  real(w(z_idx,:)))/max(abs(w(z_idx,:))), 'k-', 'color', [.5 .5 .5])
    plot(z_1D, abs(  u(z_idx,:))/max(abs(w(z_idx,:))), 'k-')
    plot(z_1D, abs(rho(z_idx,:))/max(abs(rho(z_idx,:))), 'k--')
    plot(z_1D, abs(  t(z_idx,:))/max(abs(t(z_idx,:))), 'k-.')
    plot(z_1D, abs(  p(z_idx,:))/max(abs(p(z_idx,:))), 'b-')
    plot(z_1D, max(abs(p(:,:)),[],1)/max(abs(p(:))), 'g-')
    plot(z_1D, max(abs(w(:,:)),[],1)/max(abs(w(:))), 'g--')    
    legend('$$|u''_{x}|$$','$$Re(u'')$$','$$|u''_{r}|$$','$$|\rho''|$$','$$|T''|$$','$$|p''|$$','$$\max_r|p''|$$','$$\max_r|u_x''|$$','location','northeast')
    xlabel('$$x$$'), ylabel('$$|\cdot|, Re(\cdot)$$'), axis tight

    max_p   = norm(p(:),Inf);
    max_w   = norm(w(:),Inf);
    max_val = max([max_p max_w]);
    % make nice plots and normalize bei max. p or w for comparability
    figure('name',['Polished plot for mode ' num2str(i)],'position',[575           2        1345         453])
%     subplot(2,1,2)
    plotvar = w;%(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
    varname = '$$\hat{u}$$'; 
    r_new   = r;%(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
    z_new   = z;%(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);    
    plotvar = plotvar./max_val;
    if filter_vis, plotvar = filter_x(real(plotvar),0.4); plotvar = filter_y(real(plotvar),0.4); end
    if interp_vis
        z_new           = linspace(z_vis_0,z_vis_1,n_vis_z);
        r_new           = linspace(r_vis_0,r_vis_1,n_vis_r);
        [z_new,r_new]   = meshgrid(z_new,r_new);
        plotvar         = interp2(z, r, plotvar, z_new, r_new);
        z_plot = z_new(1,:);
        r_plot = r_new(:,1);
    else
        z_plot  = z_1D(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
        r_plot  = r_1D(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right); 
        plotvar = plotvar(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
    end   
    pcolor_h=pcolor(z_plot, r_plot, real(plotvar)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on        
%     contour(z, r, real(w), isoVals,  'edgecolor','r','linewidth',0.5);
%     contour(z, r, real(w), -isoVals, 'edgecolor','b','linewidth',0.5);
    % % for acoustics
    %farfield_max = max(abs(plotvar(r_new==r_vis_1)));
    %caxis(1.5*[-farfield_max farfield_max]);
    farfield_max = max(abs(plotvar(:)));
    caxis(0.1*[-farfield_max farfield_max]);
    colormap(gray)
%     c = colorbar('eastoutside');
    c.Label.String = varname; c.TickLabelInterpreter = 'latex';
    %caxis([-.5 .5])
    xlabel('$$x$$'); ylabel('$$r$$');
    hold on
    if exist('pipeBC','var'), if pipeBC, rectangle('position',[z(pipe_idx_bottom(1)),r(pipe_idx_bottom(1)),z(pipe_idx_bottom(end))-z(pipe_idx_bottom(1)),r(pipe_idx_top(1))-r(pipe_idx_bottom(1))],'edgecolor','b','facecolor','b'), end, end
    
    xlim([z_vis_0 z_vis_1]); ylim([r_vis_0 r_vis_1]); makeFigurePublishable; makeFigurePublishable; set(gca, 'Position', get(gca, 'OuterPosition') -get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    %isoVals = linspace(0,norm(plotvar(:),Inf),8); isoVals = isoVals(2:end-1);
    %contour(z_plot, r_plot, real(plotvar), 'edgecolor','r');     
%     [p_max,p_max_i] = max(abs(p(:)));
%     plot(z(p_max_i), r(p_max_i),'ro','MarkerFaceColor','r');
%     ff_line_i       = findnearest(r(:,1), 6);
%     %plot(z(ff_line_i,:), r(ff_line_i,:),'r:');
%     [p_max_ff,p_max_ff_i] = max(abs(p(ff_line_i,:)));
%     plot(z(ff_line_i,p_max_ff_i), r(ff_line_i,p_max_ff_i),'ro');
%     p_max_ff/p_max
    
%     % figure('name',['Polished streamwise velocity plot for mode ' num2str(i)],'position',[575           2        1345         453])
%     subplot(2,1,1)
%     plotvar = real(w(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right)); varname = '$$\hat{u}_x$$';
%     r_new   = r(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
%     z_new   = z(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);    
%     plotvar = plotvar./max_val;
%     if filter_vis, plotvar = filter_x(real(plotvar),0.4); plotvar = filter_y(real(plotvar),0.4); end
%     if interp_vis
%         z_new           = linspace(z_new(1),z_new(end),5*length(z_1D));
%         r_new           = linspace(r_new(1),r_new(end),5*length(r_1D));
%         [z_new,r_new]   = meshgrid(z_new,r_new);
%         plotvar         = interp2(z(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right), r(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right), plotvar, z_new, r_new);
%         z_plot = z_new(1,:);
%         r_plot = r_new(:,1);
%     else
%         z_plot  = z_1D(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
%         r_plot  = r_1D(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right); 
%         plotvar = plotvar(1:end-Nsponge_top,Nsponge_left+1:end-Nsponge_right);
%     end   
%     pcolor_h=pcolor(z_plot, r_plot, plotvar); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on
%     colormap(gray)
%     c = colorbar('eastoutside');
%     c.Label.String = varname; c.TickLabelInterpreter = 'latex';
%     caxis([-.5 .5])
%     xlabel('$$x$$'); ylabel('$$r$$');    

    
    if render3D
        phi_1D                  = linspace(0,2*pi,Nphi_interp+1);
        phi_1D                  = phi_1D(1:end-1);
        
        [rmesh,phimesh,zmesh]   = meshgrid(r_1D,phi_1D,z_1D);
        [xmesh,ymesh,zmesh]     = pol2cart(phimesh,rmesh,zmesh);
        
        w = permute(repmat(w,[1 1 Nphi_interp]),[3 1 2]);
        w = w.*exp(1).^(1i*m.*phimesh);
        w = real(w);
        
        [xg, yg, zg]    = meshgrid(linspace(-r_1D(end),r_1D(end),Nr_interp),linspace(-r_1D(end),r_1D(end),Nr_interp),z_1D);
        gdata           = zeros(Nr_interp,Nr_interp,Nz);
        for j=1:Nz
            gdata(:,:,j)     = griddata(squeeze(xmesh(:,:,j)),squeeze(ymesh(:,:,j)),squeeze(w(:,:,j)),squeeze(xg(:,:,j)),squeeze(yg(:,:,j)));
        end
        gdata           = gdata/max(abs(gdata(:)));
        
        % plot transverse plane
        subplot(1,3,3)
        contourf(squeeze(xg(:,:,z_idx)),squeeze(yg(:,:,z_idx)),squeeze(gdata(:,:,z_idx))), axis equal, axis tight
        caxis([-1 1]);
        colormap(gray)
        
        % 3-D plot
        figure('name',['+/- 0.2*max(real(w'')) for mode ' num2str(i)])
        p = patch(isosurface(xg,yg,zg,gdata,0.05));
        isonormals(xg,yg,zg,gdata,p)
        set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5,'FaceLighting','phong');
        hold on
        % negative contour
        p = patch(isosurface(xg,yg,zg,gdata,-0.05));
        isonormals(xg,yg,zg,gdata,p)
        set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.5,'FaceLighting','phong');
        
        W0 = reshape(W, Nr, Nz);%real(reshape(eigVecs(w_j,i), Nr, Nz));
        surface(zeros(size(r)) , r , z, W0,'FaceAlpha',0.5,'LineStyle','none','FaceColor','interp')
        %     dilatation = real(reshape(Dz*eigVecs(w_j,i) + Dr*eigVecs(u_j,i) + m*eigVecs(v_j,i), Nr, Nz));
        %     p = real(reshape(eigVecs(T_j,i).*eigVecs(T_j,i)/kappa/Ma^2, Nr, Nz));
        imagw = imag(reshape(eigVecs(w_j,i), Nr, Nz));
        surface(zeros(size(r)) ,-r , z, W0,'FaceAlpha',0.5,'LineStyle','none','FaceColor','interp')
        colormap(hsv)
        
        daspect([1,1,1])
        view(3); %axis tight
        camlight
        lighting gouraud
        box on
    end
    
    drawnow
    if noEigs==1, break, end
end

