%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting TKE profiles

clear;
addpath('./aux_plots/');
dirout = './';

%% File of TKE profiles

close all;
load('./aux_plots/files/similarity_tke.mat');

%% Generating figure environment

close all;
figure;
x0=0;
y0=0;
width=15;
height=5;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height])

[ha, pos] = tight_subplot(1,2,.1,[.2,.05],[.08 .02]);

%% Plotting unnormalized TKE profiles

axes(ha(2));
hold on;

C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(11,@srgb_to_Jab);
count = 1;


count = 1;
for i = 4:2:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), smoothdata(tke_1d(1:nr,i)/max(abs(tke_1d(:,i))), 'loess',2), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), smoothdata(tke_1d(1:nr,i)/max(abs(tke_1d(:,i))),'loess',2), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 16; 

ylim([0 1.2]);
xlim([0.05 5]);
xticks([0 1 2 3 4 5]);
xticklabels({'0','1','2', '3', '4', '5'});


box on;

hXLabel = xlabel('$r/L_{k}$','interpreter','latex','fontsize',20);
hTitle  = title('(b)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);

%% Plotting normalized TKE profiles

axes(ha(1)); 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(11,@srgb_to_Jab);
count = 1;

Legend = cell(11,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 30$';
Legend{3} = '$x/D = 40$';
Legend{4} = '$x/D = 50$';
Legend{5} = '$x/D = 60$';
Legend{6} = '$x/D = 70$';
Legend{7} = '$x/D = 80$';
Legend{8} = '$x/D = 90$';
Legend{9} = '$x/D = 100$';
Legend{10} = '$x/D = 110$';
Legend{11} = '$x/D = 120$';


for i = 4:2:20
   disp(i);
   plot(rc, smoothdata(tke_1d(1:nr,i)/max(abs(tke_1d(:,i))),'loess',2), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc, smoothdata(tke_1d(1:nr,i)/max(abs(tke_1d(:,i))),'loess',2), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 20; 

xlim([0.05 10]);
xticks([0 2 4 6 8 10]);
xticklabels({'0','2','4', '6', '8', '10'});

ylim([0 1.2]);
yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',20);
hYLabel = ylabel('$K$/$K_{o}$','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize',20, 'Units', 'normalized', 'Position', [-0.1, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Units = 'normalized';
hLegend.Position = [-0.1 0.1 1 1];

%% Saving the plot
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,strcat(dirout, 'k_func_r_D_r_Lk_x_D', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'k_func_r_D_r_Lk_x_D', '.eps'),'-depsc2','-r600');