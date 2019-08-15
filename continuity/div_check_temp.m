Uc     = 1.0;
ur     = u_eigenmode(:,26,6);
utheta = v_eigenmode(:,26,6);
ux     = w_eigenmode(:,26,6);
for i = 2:nr
    mult1 = 2/(rc(i,1) + rc(i-1,1));
    mult2 = 1/(rc(i,1) - rc(i-1,1));
    mult3 = rc(i,1)*ur(i,1) - rc(i-1,1)*ur(i-1,1);
    divur(i-1,1) = mult1*mult2*mult3; 
end

for i = 2:nr
    mult1 = sqrt(-1)*mode;
    mult2 = 0.5*(utheta(i,1) + utheta(i-1,1));
    mult3 = 2/(rc(i,1) + rc(i-1,1));
    divutheta(i-1,1) = mult1*mult2*mult3;
end

for i = 2:nr
    mult1 = sqrt(-1)*2*pi*(-f(5))/(Uc);
    mult2 = 0.5*(ux(i,1) + ux(i-1,1));
    divux(i-1,1) = mult1*mult2;
end

total_div = abs(divur + divutheta + divux);
tota_div_mean = mean(total_div(2:end,1))
absr = abs(divur);
abst = abs(divutheta);
absx = abs(divux);

hold on
h1 = plot(rc(3:end,1), total_div(2:end,1), 'k*');
h2 = plot(rc(3:end,1), absr(2:end,1), 'rs');
h3 = plot(rc(3:end,1), abst(2:end,1), 'bo');
h4 = plot(rc(3:end,1), absx(2:end,1), 'yd');
xlim([0 5]);
ylim([0 10]);

hXLabel = xlabel('$r/D$','interpreter','latex','fontsize',15);
hYLabel = ylabel('$|D^{n=1, St = 0.135, m=1}|$','interpreter','latex','fontsize',15);
hTitle = title('Divergence as a function of r at $x/D = 80$','interpreter','latex','fontsize',15);

hLegend = legend([h1,h2,h3,h4],'$D_{total}$', '$D_{r}$', '$D_{\theta}$', '$D_{x}$');
hLegend.Interpreter = 'Latex';
hLegend.FontSize = 10;
hLegend.FontWeight = 'bold';
hLegend.Position = [0 0 1 1];

set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'divergence_pod_1_st_0.135_m1_x_D_80.png','-dpng','-r300');  
print(gcf,'divergence_pod_1_st_0.135_m1_x_D_80.eps','-depsc','-r600');
