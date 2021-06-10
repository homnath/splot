% Open file
inpf=fopen('../data/DB.X10.FXP.semp','r');
data_pair=fscanf(inpf,'%f %f',[2 inf])';
t=data_pair(:,1);
data=data_pair(:,2);
fclose(inpf);

dt=t(2)-t(1);

% Frequency filter parameters
%   forder  : order of butterworth  
%   ffreq   : cut off frequency(ies)
%   dt      : time increment
%   ftype   : filter type (1: low, 2: band, 3: high) 
%   fcausal : causal (1),  acausal (2)
forder=3;
ffreq=[0 0.5];
ftype=2;
fcasual=1;

data_filt=butterworth(forder,ffreq,dt,ftype,fcasual,data);

figure
hold on
plot(t,data,'k');
plot(t,data_filt,'b');