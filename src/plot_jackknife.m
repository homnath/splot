function plot_jackknife

%fin=fopen('py_3d_lookup_strike30_test_more_monte_carlo.txt','r');
%[dumc]=textread('./output/py_3d_lookup_strike30_test_more_monte_carlo.txt','%s\n',1);
x0=260.5;
y0=330.0;
z0=380.1;
t0=0.0;

[xs,ys,zs,to,mean_xs,mean_ys,mean_zs,mean_to,sd_space,sd_time,ngen,nfun]= ...
    textread('../output/2006_041_14_23_51_3_lookup_new_depar_new_optimization_result_jack_knife.txt','%f %f %f %f %f %f %f %f %f %f %d %d\n',18); 

fprintf(1,'Mean location: %.1f %.1f %.1f %.5f\n',mean(xs),mean(ys),mean(zs),mean(to));
fprintf(1,'Standard deviation: %.1f %.1f %.1f %.5f\n',std(xs),std(ys),std(zs),std(to));
errs=sqrt((xs-x0).^2+(ys-y0).^2+(zs-z0).^2);
errt=sqrt((to-t0).^2);

%figure
%plot(errs,'ok','MarkerFaceColor','k');
%hold on
%plot(sd_space,'-');

figure
set(gca,'Fontsize',14)
set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])

errorbar(errs,sd_space,'ok','MarkerSize',3,'MarkerFaceColor','k')
xlabel('Iteration number');
ylabel('Error in space location (m)');
%ylim([-200 600])
savefigure('../output/jackknife');

figure
set(gca,'Fontsize',14)
set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])
plot3(x0,y0,-z0,'o','markersize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);    
plot3(xs,ys,-zs,'o','markersize',8,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);    
xlim=([0 622]);
ylim=([0 622]);
zlim=([-622 0]);
set(gca,'box','on','xlim',xlim,'ylim',ylim,'zlim',zlim);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
%xlabel('Iteration number');
%ylabel('Error in space location (m)');
%ylim([-200 600])
savefigure('./output/jackknife3D');
return