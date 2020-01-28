clear all
% close all
clc
set (0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot, 'defaulttextInterpreter','latex')

% % % M=0.9 & M=0.4
% r_farf          = 7.0923;
% z_end           = 31.0087;
% % custom_comment  = 'M09_diag_TKE'
% % custom_comment  = 'M09_shearlayer'
% % custom_comment  = 'M09_global'
% % custom_comment  = 'M09_global_noFilter'
% custom_comment  = 'M04_global'
% % custom_comment  = 'M04_global_leftForcing10'
% % custom_comment  = 'M04_global_rightForcing10'
% % custom_comment  = 'M04_diag_TKE'%_ReStudy'
% 
% % % B118 & B122
% % r_farf          = 7.1804;
% % z_end           = 31.0677;
% % custom_comment  = '118B_global'
% % 
% calcMethod      = 'frequency_response';
% Re              = 30000;
% Nr              = 195;
% Nz              = 950;
% z_start         = -1;
% m               = 0;

% % domain truncation case: LESmeanflowJet2D_frequency_response_Re=30000_m=0_Nr=195_r=7.0923_Nz=620_-1>z>16_ReSIGMA=10.1243_ImSIGMA=0_M04_global_domainTruncation.mat
% % M=0.4
% r_farf          = 7.0923;
% z_end           = 16;
% custom_comment  = 'M04_global_domainTruncation';
% calcMethod      = 'frequency_response';
% Re              = 30000;
% Nr              = 195;
% Nz              = 620;
% z_start         = -1;
% m               = 0;

% [saveFile] = getSaveFileName('results',custom_comment,calcMethod,Re,m,Nr,Nz,r_farf,z_start,z_end,'*');

% saveFile    = './results/LESmeanflowJet2D_frequency_response_Re=10000_m=0_Nr=150_r=11_Nz=1000_-1>z>101_ReSIGMA=*_ImSIGMA=0_Cambridge_global.mat';
% Re          = 10000;

saveFile    = './results/LESmeanflowJet2D_frequency_response_Re=500_m=0_Nr=290_r=26_Nz=1100_-1>z>51_ReSIGMA=*_ImSIGMA=0_M04_diag_TKE_pressureFF.mat';
Re          = 500;

files = dir(saveFile);

load(['results/' files(1).name],'noEigs');
St    = zeros(length(files),1);
gains = zeros(length(files),noEigs);
for ii=1:length(files)
    files(ii).name;
    load(['results/' files(ii).name],'SIGMA','m','gain');
    if length(gain)>noEigs; gain=gain(1:noEigs); end
    St(ii)      = SIGMA/2/pi;
    gains(ii,:) = sort(gain,'descend');   
end
[St,sort_i]     = sort(St,'ascend');
gains           = gains(sort_i,:);

%%
figure('name',['m=' num2str(m)])
for i=noEigs:-1:1
    % line_color  = [noEigs+2 i i]/(noEigs+2); % shades of red
    % line_color  = [i i noEigs+2]/(noEigs+2); % shades of blue
    line_color  = [i i i]/(noEigs+2); % shades of gray
    subplot(1,3,1)
    loglog(St, gains(:,i).^2,'color',line_color,'linewidth',1); hold on  % /sqrt(Re) to normalize
    ylabel('$$\sigma_1^2$$'); xlabel('\emph{St}');
    subplot(1,3,2)
    loglog(St, gains(:,i).^2/Re,'color',line_color,'linewidth',1); hold on  % /sqrt(Re) to normalize 
    ylabel('$$\sigma_1^2\mbox{\emph{Re}}^{-1}$$'); xlabel('\emph{St}');
    subplot(1,3,3)
    loglog(St, gains(:,i).^2/(Re^(2/3)),'color',line_color,'linewidth',1); hold on  % /sqrt(Re) to normalize 
    ylabel('$$\sigma_1^2\mbox{\emph{Re}}^{-\frac{1}{2}}$$'); xlabel('\emph{St}');
end
%%
subplot(1,3,1)
set(gca,'xtick',[0.1 0.2 0.5 1 2 3])
area(St, gains(:,2).^2,'facecolor','w','linestyle','none'); hold on
area(St, gains(:,1).^2,'facecolor','r','linestyle','none'); hold on
chH = get(gca,'Children');
set(gca,'Children',flipud(chH))
xlim([St(1) St(end)])
set(findall(gcf,'type','Axes'),'layer','top');
grid on
grid minor

% subplot(1,2,1), ylabel('$$\lambda$$'); xlabel('\emph{St}');
% subplot(1,2,2), ylabel('$$\lambda$$'); xlabel('\emph{St}');