function compute_error(infname,p)

% infname: input file name where the optimization result is stored
% p: percentage of minimum objective function

x0=260.5;
y0=330.0;
z0=380.1;
t0=0.0;


allmem=zeros(100,4);

load(infname);

ind=allval<=0.01*p*min(allval);
sd=std(allmem(ind,:));
r=range(allmem(ind,:)).*0.5;

sd_p=sqrt(sd(1)^2+sd(2)^2+sd(3)^2);
r_p=sqrt(r(1)^2+r(2)^2+r(3)^2);
errs=sqrt((x0-xs(1))^2+(y0-xs(2))^2+(z0-xs(3))^2);
errt=abs(t0-xs(4));

fprintf(1,'Actual error: %f %f\n',errs,errt);
fprintf(1,'%.0f%% of objective function locaton uncertainty: %f %f %f %f\n',p,sd);
fprintf(1,'%.0f%% of objective function space uncertainty: %f\n',p,sd_p);
fprintf(1,'%.0f%% of objective function range: %f %f %f %f\n',p,r);
fprintf(1,'%.0f%% of objective function space range: %f\n',p,r_p);
return