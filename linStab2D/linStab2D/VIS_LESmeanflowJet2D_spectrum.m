clear all
close all
clc
set (0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
addpath('../aux');

line_specs      = {'-'; '--'; ':'; '-.'};
symbol_specs    = {'^'; 'o'; '+'; '*';  'v'; 'x'; 's'; 'd'; '>'; 'p'; 'h'; '<'};
color_specs     = {'b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'};

fig1 = figure;
fig2 = figure;

OMEGA_all   = [];
% files       = dir(['./results/*_domainStudyR_constStep.mat']);
% files       = dir(['./results/*higherAzInstab.mat']);
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=100000_m=1_Nr=250_r=2.5_Nz=1150_-0.5>z>25.5_ReSIGMA=*_ImSIGMA=0.mat']);
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=100000_m=0_Nr=250_r=*_Nz=1150_-0.5>z>*.5_ReSIGMA=6*_ImSIGMA=0*.mat']);
% files       = dir(['./results/LESmeanflowJet2D_adjoint_Re=100000_m=2_Nr=250_r=2.5_Nz=1150_-0.5>z>25.5_ReSIGMA=*_ImSIGMA=0.mat']);
% files       = dir(['./results/*direct*1200*_acoustic_branch00.mat']); files = files(1:5);
% files       = dir(['./results/*direct*_wavemaker_branch10.mat']);
% files       = dir(['./results/*sensitivity*.mat']);

% % Ma=1.5 (B118) @ Re=10,000
% m=0
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=10000_m=0_Nr=150_r=4_Nz=850_-1>z>26_ReSIGMA=*_ImSIGMA=0_B118.mat']);
% m=1
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=10000_m=1_Nr=150_r=4_Nz=850_-1>z>26_ReSIGMA=*_ImSIGMA=0_B118.mat']);
% Ma=0.9 @ Re=10,000
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=10000_m=1_Nr=150_r=4_Nz=850_-1>z>31_ReSIGMA=*_ImSIGMA=0.mat']);

% files       = dir(['./results/LESmeanflowJet2D_direct_Re=150000_m=0_Nr=250_r=2.5_Nz=1150_-0.5>z>25.5_ReSIGMA=*_ImSIGMA=0.mat']);

% files       = dir(['./results/LESmeanflowJet2D_direct_Re=100000_m=0_Nr=250_r=2.5_Nz=1150_-0.5>z>30.5_ReSIGMA=*_ImSIGMA=0.mat']);
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=100000_m=1_Nr=250_r=2.5_Nz=1150_-0.5>z>25.5_ReSIGMA=*_ImSIGMA=0.mat']);
% files       = dir(['./results/LESmeanflowJet2D_direct_Re=100000_m=0_Nr=350_r=3.5_Nz=1150_-0.5>z>25.5_ReSIGMA=*_ImSIGMA=0.mat']);

%files       = dir(['./results/LESmeanflowJet2D_direct_Re=30000_m=0_Nr=175_r=7.1804_Nz=850_-1>z>31.0677_ReSIGMA=*_ImSIGMA=0_B118.mat']);
files       = dir(['./results/LESmeanflowJet2D_direct_Re=30000_m=0_Nr=195_r=7.0923_Nz=950_-1>z>31.0087_ReSIGMA=*_ImSIGMA=0_M09.mat']);

xa  = 0;
ya  = 0;
for ii=1:length(files)
    disp([' --> processing LST file ' num2str(ii) '/' num2str(length(files))]);
    load(['results/' files(ii).name],'OMEGA','SIGMA','m');
    rand_color(ii)  = ceil(rand(1)*length(color_specs));
    rand_symbol(ii) = ceil(rand(1)*length(symbol_specs));
    if ii==1
        symbol      = {[char(strcat(color_specs(rand_color(ii)), symbol_specs(rand_symbol(ii))))]};
        lgnd        = {['m = ' num2str(m)]};
    else
        symbol(ii)  = {[char(strcat(color_specs(rand_color(ii)), symbol_specs(rand_symbol(ii))))]};
        lgnd(ii)    = {['m = ' num2str(m)]};
    end
    
    figure(fig1) % colorful plot for interpretation
    plot(real(SIGMA),imag(SIGMA),char(symbol(ii)),'LineWidth',2,'MarkerSize',10); hold on
    legend(lgnd)
    % plot spectrum
    plot(real(OMEGA),imag(OMEGA),char(symbol(ii)))
    % draw Arnoldi convergence area
    dist_SIGMA_max = max(abs(OMEGA-SIGMA));
    pos = [(real(SIGMA)-dist_SIGMA_max) imag(SIGMA)-dist_SIGMA_max 2*dist_SIGMA_max 2*dist_SIGMA_max];
    rectangle('Position',pos,'Curvature',[1 1],'LineWidth',0.5, 'LineStyle','--','EdgeColor',char(color_specs(rand_color(ii))))
    set(gca,'Layer','top')
    text(real(SIGMA),imag(SIGMA),[num2str(real(SIGMA)) '+' num2str(imag(SIGMA)) 'i']);
    xlabel('$$\omega_r$$'), ylabel('$$\omega_i$$')
    axis tight
    this_xlims = xlim;
    
    figure(fig2) % plain style plot
    OMEGA_all = [OMEGA_all; OMEGA];
    pos = [(real(SIGMA)-dist_SIGMA_max)/2/pi imag(SIGMA)-dist_SIGMA_max 2*dist_SIGMA_max/2/pi 2*dist_SIGMA_max];
    rectangle('Position',pos,'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold on
    
    theta   = linspace(0, 2*pi, 100);
    x1      = dist_SIGMA_max/2/pi*cos(theta)  + real(SIGMA)/2/pi;
    y1      = -dist_SIGMA_max*sin(theta)      + imag(SIGMA);
       
    [xa, ya] = polybool('union', x1, y1, xa, ya);  
end
OMEGA_all = unique(round(OMEGA_all,5));
%%
figure(fig2) % plain style plot
% h = rectangle('Position',[0 0 10 10],'LineStyle','none','FaceColor',[.8 .8 .8]); %uistack(h,'bottom')%
plot(real(OMEGA_all)/2/pi,imag(OMEGA_all),'.','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
set(gca,'Layer','top')
xlabel('St'), ylabel('$$\omega_i$$')
axis tight
ylim([min(imag(OMEGA_all)) max(imag(OMEGA_all))+0.05])
xlim([0 max(real(OMEGA_all)/2/pi)])
%set(gca,'Color',[0.7 0.7 0.7]);
hp = patch([0;10;10;0],[-10;-10;10;10],1); uistack(hp,'bottom')
hp = hatchfill(hp);uistack(hp,'bottom')
patch(xa, ya, 1, 'FaceColor', 'none', 'EdgeColor', 'k')
plot(xlim,[0 0],'k--')
box on



















