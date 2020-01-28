function vis_baseFlow(grid,sponge,vars)

figure('name','Grid')
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.W,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), hold on
mesh(grid.x,grid.r,zeros(grid.Nr,grid.Nx),'FaceColor','none','edgecolor','k'); axis equal tight
plot([sponge.x_inlet sponge.x_inlet],[0 sponge.r_farfield],'r--');
plot([sponge.x_inlet sponge.x_outlet],[sponge.r_farfield sponge.r_farfield],'r--');
plot([sponge.x_outlet sponge.x_outlet],[0 sponge.r_farfield],'r--');

figure('name','Base flow')
subplot(3,2,1)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.NU,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\nu$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','cap','HorizontalAlignment','left'), colorbar
subplot(3,2,2)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.W,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$u_x$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,3)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.U,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$u_r$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,4)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.P,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$p$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,5)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.T,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$T$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,6)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.spongeFun,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), 'sponge', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar

figure('name','Base flow derivatives')
subplot(3,2,1)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.DIL,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' dilatation', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','cap','HorizontalAlignment','left'), colorbar
subplot(3,2,2)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.LAP_NU,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\nabla^2\nu$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,3)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.dWdr,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\partial u_x/\partial r$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,4)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.dWdx,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\partial u_x/\partial x$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,5)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.d2Wdr2,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\partial^2 u_x/\partial r^2$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar
subplot(3,2,6)
pcolor_h=pcolor(grid.x_1D, grid.r_1D, reshape(vars.d2Wdx2,grid.Nr,grid.Nx)); set(pcolor_h,'EdgeColor','none'), axis equal, axis tight, hold on, text(grid.x_1D(1), grid.r_1D(end), ' $$\partial^2 u_x/\partial x^2$$', 'BackgroundColor',[.7 .9 .7],'VerticalAlignment','Top','HorizontalAlignment','left'), colorbar

drawnow

