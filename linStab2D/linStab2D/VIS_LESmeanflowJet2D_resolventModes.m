% OTS, 2018
clear variables
close all
clc

%%%%%%%%%%%%%%
% user input %
%%%%%%%%%%%%%%
%file        = 'results/LESmeanflowJet2D_frequency_response_Re=10000_m=14_Nr=215_r=7.0923_Nz=950_-1>z>41.0087_ReSIGMA=1.2566_ImSIGMA=0_M04_global_streak.mat';
file        = 'results/LESmeanflowJet2D_frequency_response_Re=500_m=0_Nr=200_r=28_Nz=600_-1>z>53_ReSIGMA=6.2832_ImSIGMA=0_M04_diag_TKE_pressureFF.mat';
var_forcing = 'u_x'; % valid choices are 'p', 'u_x', 'u_r', 'u_t', 'rho', 'T'
var_response= 'p';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(file)

%%%%%%%%%%%%%%%%%%%%%%
% plot gain spectrum %
%%%%%%%%%%%%%%%%%%%%%%
figure
stem(gain)
xlabel('mode')
ylabel('gain')
set(gca,'XTick',1:noEigs), xlim([0.5 noEigs+0.5])

%%%%%%%%%%%%%%
% plot modes %
%%%%%%%%%%%%%%
figure
noEigs  = size(V_in,2);
count   = 1;
for i = 1:noEigs
    switch var_forcing
        case 'p'
            rho_f   = reshape(V_in(rho_j,i), Nr, Nz);
            t_f     = reshape(V_in(T_j,i), Nr, Nz);
            forcing = (reshape(RHO, Nr, Nz).*t_f + rho_f.*reshape(T, Nr, Nz))/kappa/Ma^2;
        case 'u_x'
            forcing = reshape(V_in(w_j,i), Nr, Nz);
        case 'u_r'
            forcing = reshape(V_in(u_j,i), Nr, Nz);           
        case 'u_t'
            forcing = reshape(V_in(v_j,i), Nr, Nz);        
        case 'T'
            forcing = reshape(V_in(T_j,i), Nr, Nz);           
        case 'rho'
            forcing = reshape(V_in(rho_j,i), Nr, Nz);         
        otherwise
            disp('please choose valid variable'); break
    end
    switch var_response
        case 'p'
            rho_q   = reshape(U_out(rho_j,i), Nr, Nz);
            t_q     = reshape(U_out(T_j,i), Nr, Nz);
            response= (reshape(RHO, Nr, Nz).*t_q + rho_q.*reshape(T, Nr, Nz))/kappa/Ma^2;
        case 'u_x'
            response= reshape(U_out(w_j,i), Nr, Nz);
        case 'u_r'
            response= reshape(U_out(u_j,i), Nr, Nz);            
        case 'u_t'
            response= reshape(U_out(v_j,i), Nr, Nz);            
        case 'T'
            response= reshape(U_out(T_j,i), Nr, Nz);            
        case 'rho'
            response= reshape(U_out(rho_j,i), Nr, Nz);            
        otherwise
            disp('please choose valid variable'); break
    end
    
    subplot(noEigs,2,count)
    pcolor(z,r,real(forcing)); axis equal, axis tight
    caxis(0.5*max(abs(caxis))*[-1 1])
    hold on
    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
    xlim([z_sponge1 z_sponge2]);
    ylim([r(1) r_sponge]);
    text(z_sponge1,r_sponge,['forcing mode ' num2str(i)],'VerticalAlignment','top')
    colormap(gray), shading interp
    ylabel('r')
    if i<noEigs, set(gca,'xticklabel',[]); else, xlabel('x'), end
    count   = count + 1;
    
    subplot(noEigs,2,count)
    pcolor(z,r,real(response)); axis equal, axis tight
    caxis(0.5*max(abs(caxis))*[-1 1])
    hold on
    contour(z,r,reshape(W,Nr,Nz),[0.99 0.99],'edgecolor','r','LineStyle','-');
    contour(z,r,reshape(W,Nr,Nz),[0.1  0.1], 'edgecolor','r','LineStyle','--');
    xlim([z_sponge1 z_sponge2]);
    ylim([r(1) r_sponge]);
    text(z_sponge1,r_sponge,['response mode ' num2str(i)],'VerticalAlignment','top')
    colormap(gray), shading interp
    ylabel('r')
    if i<noEigs, set(gca,'xticklabel',[]); else, xlabel('x'), end
    count   = count + 1;
end