function plot_monte

%fin=fopen('py_3d_lookup_strike30_test_more_monte_carlo.txt','r');
%[dumc]=textread('./output/py_3d_lookup_strike30_test_more_monte_carlo.txt','%s\n',1);
x0=260.5;
y0=330.0;
z0=380.1;
t0=0.0;

[xs,ys,zs,to,mean_xs,mean_ys,mean_zs,mean_to,sd_space,sd_time,ngen,nfun]= ...
    textread('./output/py_3d_lookup_strike30_new_depar_monte_monte_carlo.txt','%f %f %f %f %f %f %f %f %f %f %d %d\n',100); 

errs=sqrt((xs-x0).^2+(ys-y0).^2+(zs-z0).^2);
errt=sqrt((to-t0).^2);

allmem=zeros(100,4);
p=0.95; %0 to 1.0
outf=fopen('./output/2006_041_14_23_51_3_lookup_new_depar1_optimization_result_jack_knife.txt','w');
for i=1:18
    optf=strcat('./output/2006_041_14_23_51_3_lookup_new_depar1_optimization_result_jack_knife',num2str(i),'.mat');
    load(optf);
    
    ind=allval<=p*min(allval);
    sd=std(allmem(ind,:));
    r=range(allmem(ind,:)).*0.5;
    
    sd_p(i)=sqrt(sd(1)^2+sd(2)^2+sd(3)^2);
    r_p(i)=sqrt(r(1)^2+r(2)^2+r(3)^2);
    
    boundx=max(allmem(ind,1))-min(allmem(ind,1));
    boundy=max(allmem(ind,2))-min(allmem(ind,2));
    boundz=max(allmem(ind,3))-min(allmem(ind,3));
    bound(i)=sqrt(boundx^2+boundy^2+boundz^2);
end

figure
plot(errs,'-ok');
%hold on
%plot(sd_p,'--g');
%hold on
%plot(2*sd_p,'-g');
hold on
plot(r_p,'-ob');
legend('Actual error','scaled range');
error completed

figure
plot(errs,'-ok','MarkerFaceColor','k');
hold on
plot(sd_space,'-');


figure
set(gca,'Fontsize',14)
set(gcf,'Paperunits','centimeters','Paperposition',[1 1 14 10])

errorbar(errs,sd_space,'ok','MarkerSize',3,'MarkerFaceColor','k')
xlabel('Iteration number');
ylabel('Error in space location (m)');
%ylim([-200 600])
savefigure('./output/monte_carlo');

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
savefigure('./output/monte_carlo_3D');
return