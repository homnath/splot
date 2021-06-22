

% Open file
inpf=fopen('/Users/madisonreed/Work/Projects/Examples/36mesh_25000steps_f0=0.01/DB.X10.FXP.semp','r');
data_pair=fscanf(inpf,'%f %f',[2 inf])';
t=data_pair(:,1);
data=data_pair(:,2);
fclose(inpf); 

dt=t(2)-t(1);


% Frequency filter parameters
%   forder  : order of butterworth  
%   ffreq   : cut off frequency(ies)
%             Low pass  : ffreq(2)
%             Band pass : ffreq
%             High Pass : ffreq(1)
%   dt      : time increment
%   ftype   : filter type (1: low, 2: band, 3: high) 
%   fcausal : causal (1),  acausal (2)
<<<<<<< HEAD

forder=3;
=======
forder=-3;
>>>>>>> 5b367941d6443f60a26babdcfd7fd452d65d9ee5
fcasual=1;


% Low pass
ftype=1;
ffreq=[0 500];
data_low=butterworth(forder,ffreq,dt,ftype,fcasual,data);

% Band pass
ftype=2;
ffreq=[10 500];
data_band=butterworth(forder,ffreq,dt,ftype,fcasual,data);

% High pass
ftype=3;
ffreq=[10 0];
data_high=butterworth(forder,ffreq,dt,ftype,fcasual,data);

figure
hold on
plot(t,data,'k');
%plot(t,data_low,'b');
%plot(t,Pressure1)
plot(t,data_band,'g');
<<<<<<< HEAD
%plot(t,data_high,'r');
=======
plot(t,data_high,'r');

>>>>>>> 5b367941d6443f60a26babdcfd7fd452d65d9ee5
xlabel('Time (s)')
ylabel('Pressure (Pa)')
legend('Original Data','Low Pass','Band Pass','High Pass')