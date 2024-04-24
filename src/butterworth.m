%==========================================================================
% function filtered_data = butterworth(order,freq,dt,filtype,causal_flag,data)
% lowpass butterworth filter
% INPUT : 
%         order         : order of butterworth  
%         freq          : cut off frequency(ies)
%         dt            : time increment
%         filtype       : filter type (1: low, 2: band, 3: high) 
%         causal_flag   : causal (1),  acausal (2)
%         data                  : matrix with seismograms as columns
% OUTPUT: 
%         filtered_data         : matrix with filtered seismograms 
%
% NOTE!  acausal filtering increases order by factor 2 
%
% 27.01.04 MRO
% maximum upper cutoff frequency is set to 0.9 * f_nyquist 
%
% Michael Roth 25.08.00
%==========================================================================  

function filtered_data = butterworth(order,freq,dt,filtype,causal_flag,data)

% This will ensures the normalized frequency bounds to [0 1]
% we take the values closer to 0 and 1 to avoid numerical instability
minf=0.0000001;
maxf=0.999999;

freq=freq*dt*2.;
freq(1)=max(freq(1),minf);
freq(2)=min(freq(2),maxf);
if filtype==1
  tfreq=freq(2); 
  [B,A]=butter(order,tfreq);
elseif filtype==2
  tfreq=freq;
  [B,A]=butter(order,tfreq);
elseif filtype==3
  tfreq=freq(1);
  [B,A]=butter(order,tfreq,'high');
end

if (causal_flag==1)
 filtered_data=filter(B,A,data);
elseif (causal_flag==2)
 filtered_data=filtfilt(B,A,data);  
end


