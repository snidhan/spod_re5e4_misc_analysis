%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting the bar plots

clear;
addpath('./aux_plots/');
dirout = './';

%% Reading the filename

filename = './aux_plots/files/eigvalues_similarity_diff_loc.mat';
load(filename);

%% Setting up the figure environment

close all;
figure;
x0=0;
y0=0;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.15 .02]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% Similarity of m = 2, St = 0

axes(ha(1));
C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';


count = 1;
Nf_sampled = 100;
hold on;

for i = 4:2:20
   disp(i);
%    plot(f(1:Nf_sampled)*(LK_TKE_loc_planes(i,2)/TKE_centerline_loc_planes(i,2)^0.5), eigenspectra_allm(:,1,3,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count});
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,3,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count}, 'markersize', 10);

   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0  0.5])
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'0','0.1', '0.2', '0.3', '0.4', '0.5'});
ylim([0 0.2]);
yticks([0 0.05 0.1 0.15 0.20]);
yticklabels({'0','0.05', '0.10', '0.15', '0.20'});

box on;

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\lambda^{(1)}(m=2, \mbox{\textit{St}})/(K_{o}^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.16, 0.98, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_m2_st0_eigvalue_x_D_klk', '.eps'),'-depsc2','-r600');

%% Similarity analysis of the m=1, St=0.136 leading order SPOD mode

axes(ha(2));

C = {'bo','r*','kd','c>','g<','b+','r^','kv','gp','ch'}; % Cell array of markerstyle

count = 1;

Legend = cell(9,1);
Legend{1} = 'x/D = 20';
Legend{2} = 'x/D = 30';
Legend{3} = 'x/D = 40';
Legend{4} = 'x/D = 50';
Legend{5} = 'x/D = 60';
Legend{6} = 'x/D = 70';
Legend{7} = 'x/D = 80';
Legend{8} = 'x/D = 90';
Legend{9} = 'x/D = 100';


count = 1;
Nf_sampled = 100;
hold on;
for i = 4:2:20
   disp(i);
   plot(f(1:Nf_sampled), eigenspectra_allm(:,1,2,i)/((TKE_centerline_loc_planes(i,2)^0.5)*LK_TKE_loc_planes(i,2))^2, C{count}, 'markersize', 10);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0  0.5])
xticks([0 0.1 0.2 0.3 0.4 0.5]);
xticklabels({'0','0.1', '0.2', '0.3', '0.4', '0.5'});
ylim([0 0.2]);
yticks([0 0.05 0.1 0.15 0.20, 0.25]);
yticklabels({'0','0.05', '0.10', '0.15', '0.20', '0.25'});


box on;

hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$\lambda^{(1)}(m=1,\mbox{\textit{St}})/(K_{o}^{1/2}L_{k})^{2}$','interpreter','latex','fontsize',20,'Units', 'normalized', 'Position', [-0.15, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.16, 0.98, 0]);

% hXLabel = xlabel('$\mbox{\textit{St}}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$\lambda^{(1)}(m=1, \mbox{\textit{St}}=0.135)/(U_{d}L_{d})^{2}$','interpreter','latex','fontsize',15);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%% Saving the plots
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'similarity_m1_st0136_m2_st0_eigvalue_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_m1_st0136_m2_st0_eigvalue_x_D', '.eps'),'-depsc2','-r600');
