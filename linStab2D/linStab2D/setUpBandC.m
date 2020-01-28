
function [B,C] = setUpBandC(r_farf,Nr,z1,z2,Nz,r,z,Nsponge_left,Nsponge_right,Nsponge_top,spongeFun,RHO,U,V,W,T,MU,T_0,rho_0,mu_0,p_0,nu_0,R,kappa,Pr,c,Re,Ma,u_0,OUTPUT_z_min,OUTPUT_z_max,OUTPUT_r_min,OUTPUT_r_max,INPUT_z_min,INPUT_z_max,INPUT_r_min,INPUT_r_max,forcingType,responseType,baseState,vizBaseAndRes)

Nzr         = Nr*Nz;

%%%%%%%%%%%%%%%%%%
% FORCING MATRIX %
%%%%%%%%%%%%%%%%%%

switch forcingType
    case {'core','jet','shearlayer','farfield','global'}
        u_x_mean    = reshape(W,Nr,Nz);
        delta       = 0.05;
        switch forcingType
            case 'global'
                B_win       = ones(Nr,Nz);
            case 'core'
                B_win       = u_x_mean>(1-delta);
            case 'jet'
                B_win       = u_x_mean>delta;
            case 'shearlayer'
                B_win       = u_x_mean<(1-delta)&u_x_mean>delta;
            case 'farfield'
                B_win       = u_x_mean<delta;
        end
        % smoothing
        nSmooth     = 3;
        B_win       = B_win+1;
        for si = 1:nSmooth
            for ri = 1:Nr
                B_win(ri,:)     = mvg(B_win(ri,:));
            end
            for zi = 1:Nz
                B_win(:,zi)     = mvg(B_win(:,zi));
            end
        end
        B_win       = B_win-1;
        B_win(z<INPUT_z_min|z>INPUT_z_max|r<INPUT_r_min|r>INPUT_r_max) = 0;
        
        if vizBaseAndRes, figure('name','forcing mask')
            pcolor(z,r,B_win), shading interp, colormap(flipud(gray)), axis equal tight
        end
        
        B_win       = reshape(sparse(B_win),Nzr,1);
        B = [ ...
            diag(B_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(B_win).*speye(Nzr)
            ]; 
    case 'global_momentum'
        B_win       = z>INPUT_z_min&z<INPUT_z_max&r>INPUT_r_min&r<INPUT_r_max;
        
        % smoothing
        nSmooth     = 0;
        B_win       = B_win+1;
        for si = 1:nSmooth
            for ri = 1:Nr
                B_win(ri,:)     = mvg(B_win(ri,:));
            end
            for zi = 1:Nz
                B_win(:,zi)     = mvg(B_win(:,zi));
            end
        end
        B_win       = B_win-1;
        
        B_win       = reshape(sparse(B_win),Nzr,1);
        % momentum only input
        B = [ ...
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(B_win).*speye(Nzr)     0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            ];
    case 'diag_TKE'
        % TKE sqared (so that dimensions match)
        if strcmp(baseState,'118B')
            % new B118 case with TKE as mask
            root    = '/mnt/SSD_RAID_0/B118/CylGrid_B118';
            [x_mean,r_mean,~,~,~,nx,nr,nth]    = read_CylGrid_698_136_128_xyz(root);
            [~,~,~,~,~,~,R_rr,R_rt,R_rz,R_tt,R_tz,R_zz] = read_CylGrid_698_136_128_stats(root,nx,nr,nth);
            k = 1/2*(R_rr+R_tt+R_zz);
            k = k./max(abs(k(:)));
            k = interp2(x_mean,r_mean,k,z,r,'spline',0);
        elseif strcmp(baseState,'M09')
            load('baseFlows/CylGrid_656_138_longTimeFavreMean.mat');
            reynolds_xx = interp2(x_mean,r_mean,reynolds_xx,z,r,'spline',0);
            reynolds_tt = interp2(x_mean,r_mean,reynolds_tt,z,r,'spline',0);
            reynolds_rr = interp2(x_mean,r_mean,reynolds_rr,z,r,'spline',0);
            k = (reynolds_rr+reynolds_tt+reynolds_xx);
        elseif strcmp(baseState,'M04')  
            [~,~,~,~,~,~,TKE] = get_LES_meanFlow_M04(z,r,5,5,0,[],[]);
            k   = TKE;
        elseif strcmp(baseState,'M09laminar')
            [~,~,~,~,~,~,~,TKE] = get_LES_meanFlow_M09laminar(z,r,5,5,0,[],[]);
            k   = TKE;
        else
            error('TKE forcing on diagonal not available for this base-state!');
        end
        
        % smoothing
        nSmooth     = 3;
        k           = k+1;
        for si = 1:nSmooth
            for ri = 1:Nr
                k(ri,:)     = mvg(k(ri,:));
            end
            for zi = 1:Nz
                k(:,zi)     = mvg(k(:,zi));
            end
        end
        k   = k-1;
        k   = k/max(k(:));
        
        k(z<INPUT_z_min|z>INPUT_z_max|r<INPUT_r_min|r>INPUT_r_max) = 0;
        
        % restrict to >1% to completely avoid far-field forcing 
        k(k<0.01)   = 0;
        
        if vizBaseAndRes, figure('name','forcing mask')
            pcolor(z,r,k), shading interp, colormap(flipud(gray)), axis equal tight
        end
        
        % f proportional to k
        B = [ ...
            0*diag(sparse(k(:)))        0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)
            0*speye(Nzr)                diag(sparse(k(:)))              0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                    diag(sparse(k(:)))              0*speye(Nzr)                    0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                    0*speye(Nzr)                    diag(sparse(k(:)))              0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)                    0*diag(sparse(k(:)))
            ];
    case 'RST'  % Reynolds stress tensor
        if strcmp(baseState,'M09')
            load('baseFlows/CylGrid_656_138_longTimeFavreMean.mat');
            reynolds_xx = interp2(x_mean,r_mean,reynolds_xx,z,r,'spline',0);
            % reynolds_xt = interp2(x_mean,r_mean,reynolds_xt,z,r,'spline',0);
            reynolds_xr = interp2(x_mean,r_mean,reynolds_xr,z,r,'spline',0);
            reynolds_tt = interp2(x_mean,r_mean,reynolds_tt,z,r,'spline',0);
            % reynolds_rt = interp2(x_mean,r_mean,reynolds_rt,z,r,'spline',0);
            reynolds_rr = interp2(x_mean,r_mean,reynolds_rr,z,r,'spline',0);
        else
            error('TKE forcing on diagonal not available for this base-state!');
        end
        % f proportional to RST
        B = [ ...
            0*speye(Nzr)                0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)
            0*speye(Nzr)                diag(sparse(reynolds_rr(:)))    0*speye(Nzr)                    diag(sparse(reynolds_xr(:)))    0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                    diag(sparse(reynolds_tt(:)))    0*speye(Nzr)                    0*speye(Nzr)
            0*speye(Nzr)                diag(sparse(reynolds_xr(:)))    0*speye(Nzr)                    diag(sparse(reynolds_xx(:)))    0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)                    0*speye(Nzr)
            ];
    otherwise
        error('Unexpected forcingType!')
end

%%%%%%%%%%%%%%%%%%%%%
% OBSERVABLE MATRIX %
%%%%%%%%%%%%%%%%%%%%%

nSmooth     = 1;
C_win       = z>OUTPUT_z_min&z<OUTPUT_z_max&r>OUTPUT_r_min&r<OUTPUT_r_max;
C_win       = C_win+1;
for si = 1:nSmooth
    for ri = 1:Nr
        C_win(ri,:)     = mvg(C_win(ri,:));
    end
    for zi = 1:Nz
        C_win(:,zi)     = mvg(C_win(:,zi));
    end
end
C_win       = C_win-1;
C_win       = reshape(sparse(C_win),Nzr,1);

switch responseType
    case 'global'
        C = [ ...
            diag(C_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                diag(C_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                diag(C_win).*speye(Nzr)     0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(C_win).*speye(Nzr)     0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(C_win).*speye(Nzr)
            ];
    case 'pressure'
        % % pressure output
        T   = sparse(T);
        RHO = sparse(RHO);
        C = [ ...
            diag(C_win).*diag(T)*R      0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                diag(C_win).*diag(RHO)*R
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)                0*speye(Nzr)
            ];
    otherwise
        error('Unexpected responseType!')
end