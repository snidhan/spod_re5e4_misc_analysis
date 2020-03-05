%% Name - Sheel Nidhan
%  Date - 12 Dec 2019

% Code for plotting TKE profiles

clear;
addpath('./aux_plots/');
dirout = './';

%% File of TKE profiles

close all;
load('./aux_plots/files/similarity_uxur.mat');
reystress_uw_1d_smooth = smoothdata(reystress_uw_1d,'loess',2);
reystress_uw_1d = reystress_uw_1d_smooth;

%% Generating the figure object for the plot

x0=0;
y0=0;
width=10;
height=10;
set(gcf, 'units', 'inches', 'position',[x0,y0,width,height]);

[ha, pos] = tight_subplot(3,2,[.05, 0.05],[.1,.05],[.1 .02]);
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(ha(1),'Position',[new,pos{1}(2:end)])
delete(ha(2));

%% Plotting <u_x u_r> as a function of r/D from x/D = 20 to 120

axes(ha(1));

%C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
%color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

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

hold on;
for i =  4:2:20
   disp(i);
   plot(rc, -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc, -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 15;

xlim([0 10]);
xticks([0 2 4 6 8 10]);
xticklabels({'0','2','4', '6', '8', '10'});

ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(a)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%% Plotting <u_xu_r> as a function of r/Ld from x/D = 20 to 60

axes(ha(3));

hold on;
C = {'-','-','-','-','-','-','-','-','-','-'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

%lineStyles = linspecer(10);
lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 25$';
Legend{3} = '$x/D = 30$';
Legend{4} = '$x/D = 35$';
Legend{5} = '$x/D = 40$';
Legend{6} = '$x/D = 45$';
Legend{7} = '$x/D = 50$';
Legend{8} = '$x/D = 55$';
Legend{9} = '$x/D = 60$';

count = 1;

for i = 4:12
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

% for i =  21:22
%    disp(i);
%    plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end

ax = gca;
ax.FontSize = 15;

xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca, 'Xticklabel', []);
% xticks([0 1 2 3 4]);
% xticklabels({'0','1','2', '3', '4'});
ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

%hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$-\langle u_{x}u_{r}\rangle$/max$(-\langle u_{x}u_{r}\rangle)_{r}$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(b)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_20_60', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_20_60', '.eps'),'-depsc2','-r600');


%% Plotting <u_xu_r> as a function of r/Lk from x/D = 20 to 60

axes(ha(4));

hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 20$';
Legend{2} = '$x/D = 25$';
Legend{3} = '$x/D = 30$';
Legend{4} = '$x/D = 35$';
Legend{5} = '$x/D = 40$';
Legend{6} = '$x/D = 45$';
Legend{7} = '$x/D = 50$';
Legend{8} = '$x/D = 55$';
Legend{9} = '$x/D = 60$';

count = 1;
for i = 4:12
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

% for i =  21:22
%    disp(i);
%    plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
%    count = count + 1;
% end

xlim([0 4]);
xlim([0 4]);
xticks([0 1 2 3 4]);
set(gca, 'Xticklabel', []);
ylim([0 1.2]);

% ax = gca;
% ax.FontSize = 16; 

box on;

% hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',15);
% hYLabel = ylabel('$-\langle u_{x}u_{r}\rangle$/max$(-\langle u_{x}u_{r}\rangle)_{r}$','interpreter','latex','fontsize',15);
hTitle  = title('(c)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);


hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk_20_60', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_lk_20_60', '.eps'),'-depsc2','-r600');

%% Plot of <u_{x}u_{r}> using r/Ld from x/D = 70 to 120

axes(ha(5));

hold on;
C = {'-','-','-','-','-','-','-','-','-','-'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

%lineStyles = linspecer(10);
lineStyles = maxdistcolor(9,@srgb_to_Jab);
count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 70$';
Legend{2} = '$x/D = 75$';
Legend{3} = '$x/D = 80$';
Legend{4} = '$x/D = 85$';
Legend{5} = '$x/D = 90$';
Legend{6} = '$x/D = 95$';
Legend{7} = '$x/D = 100$';
Legend{8} = '$x/D = 110$';
Legend{9} = '$x/D = 120$';

count = 1;

for i = 14:20
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_mean_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end


ax = gca;
ax.FontSize = 15;

xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4'});

ylim([0 1.2]);
yticks([0 .2 .4 .6 .8 1]);
yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$\eta_{d} = r/L_{d}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
hYLabel = ylabel('$-\langle u_{x}''u_{r}''\rangle$/max$(-\langle u_{x}''u_{r}''\rangle)_{r}$','interpreter','latex','fontsize',15,...
    'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(d)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_70_120', '.png'),'-dpng2','-r600');  
% print(gcf,strcat(dirout, 'similarity_uxur_x_D_ld_70_120', '.eps'),'-depsc2','-r600');

%% Plot of <u_{x}u_{r}> using r/Lk from x/D = 70 to 120

axes(ha(6));
 
hold on;
C = {'-','-','-','-','-','--','--','--','--','--'}; % Cell array of linestyle
color = {'k','b','r','g','m','m','g','r','b','k'}; % Cell array of linestyle

lineStyles = maxdistcolor(9,@srgb_to_Jab);

count = 1;

Legend = cell(9,1);
Legend{1} = '$x/D = 70$';
Legend{2} = '$x/D = 75$';
Legend{3} = '$x/D = 80$';
Legend{4} = '$x/D = 85$';
Legend{5} = '$x/D = 90$';
Legend{6} = '$x/D = 95$';
Legend{7} = '$x/D = 100$';
Legend{8} = '$x/D = 110$';
Legend{9} = '$x/D = 120$';

count = 1;
for i = 14:20
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

for i =  21:22
   disp(i);
   plot(rc/LK_TKE_loc_planes(i,2), -reystress_uw_1d(1:nr,i)/max(abs(reystress_uw_1d(:,i))), '-', 'Color', lineStyles(count,:), 'Linewidth',2);
   count = count + 1;
end

ax = gca;
ax.FontSize = 15;
xlim([0 4]);
xticks([0 1 2 3 4]);
xticklabels({'0','1','2', '3', '4'});
ylim([0 1.2]);
% yticks([0 .2 .4 .6 .8 1]);
% yticklabels({'0','0.2','0.4', '0.6', '0.8', '1'});

box on;

hXLabel = xlabel('$\eta_{k} = r/L_{k}$','interpreter','latex','fontsize',15,'Units', 'normalized', 'Position', [0.5, -0.10, 0]);
% hYLabel = ylabel('$-\langle u_{x}^{''}u_{r}^{''}\rangle$/max$(-\langle u_{x}^{''}u_{r}^{''}\rangle)_{r}$','interpreter','latex','fontsize',10,...
%     'Units', 'normalized', 'Position', [-0.08, 0.45, 0]);
hTitle  = title('(e)','interpreter','latex','fontsize', 15, 'Units', 'normalized', 'Position', [-0.08, 0.95, 0]);

hLegend = legend(Legend);
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 15;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

%% Saving the figure
set(gcf, 'PaperPositionMode', 'auto');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_complete', '.png'),'-dpng2','-r600');  
print(gcf,strcat(dirout, 'similarity_uxur_x_D_complete', '.eps'),'-depsc2','-r600');
