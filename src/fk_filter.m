function [newdata,mask] = fk_filter(segy_file,rec,trange,freq,norm)
%INPUT
%   segy_file: segy file name
%   rec: range of receiver, e.g., [136 161] will process receivers from 136 to
%       161 only
%   trange: time range, e.g., [0 3] will process only the time segment of 0 to
%       3s
%   freq: bandpass filter frequency, e.g., [5 20] will filter the data
%       in the band 5-20Hz
%   norm: normalization scheme
%       0: no normalization, 1: trace normalization, 2: absolute normalization, 3:
%       normalization to trace mean and absolute normalization
% function [newdata,mask] = fkdesign(data,t,offset)
% Program to plot a set of data in the f-k domain then design a filter
% which is manually picked using the mouse.
% data are the data in column oriented format
% t is the timebase in seconds, this is a vector of the same length as each trace.
% offset is the offset in metres and is a vector equal to the number of data
% 
% You will probably need to run this a number of times iteratively updating
% what you have surgically removed.
%
% The f-k image shown may first be zoomed if required, then picking begins
% when any key is hit.
%
% User may begin selection of area immediately then finish polgonal region
% picking by a shift-click,  a right-click, or a double-click adds
% final vertex to the selection and then starts making the mask which
% is then applied to the data.  The program then outputs the filtered
% data plus the mask used.
%
% A hamming taper is used to mask the data.  The data are further padded with 
% zeros to the next highest power of 2 - this implementation runs more 
% quickly than earlier ones.
%
% Students should feel free to modify the program to their own desires.
% One improvement centres on the tapering of the filter mask.

% Written originally by D. Schmitt, last modified January 27, 2000.
% Modified by Hom N Gharti while at NORSAR

[Data,SegyTraceHeaders,SegyHeader]=ReadSegy(segy_file);

dt=SegyHeader.time(2)-SegyHeader.time(1);
ntot=SegyHeader.ns;
tmin=SegyHeader.time(1);

nrec=rec(2)-rec(1)+1;
if nrec<1
    error('there is no receiver!\n');
end
v=3000;

% receiver
dx=50;
offset=(0:nrec-1)'*dx;

n1=round((trange(1)-tmin)/dt)+1;
if(n1<1)
    n1=1;
end
n2=round((trange(2)-tmin)/dt)+1;
if n2>ntot
    n2=ntot;
end
ntot=n2-n1+1;
t=(0:ntot-1)'*dt+trange(1);

data=Data(n1:n2,rec(1):rec(2));
clear Data

% Filtering
for i_rec = 1:nrec
    data(:,i_rec)=buttwo(3,[freq(1) freq(2)],dt,2,1,data(:,i_rec));
end


% normalization
if norm==1
    % trace normalization
    for i_rec=1:nrec
        data(:,i_rec)=data(:,i_rec)/max(abs(data(:,i_rec)));
    end
end
if norm==3
    % normalization to trace mean
    for i_rec=1:nrec
        data(:,i_rec)=data(:,i_rec)/mean(abs(data(:,i_rec)));
    end
end
   
if norm==2 || norm==3
    % absolute normalization
    data=data/max(max(abs(data)));
end

%t=0:dt:(ntot-1)*dt;t=t'+min(min(tr0)); %,min(tr1));
plot_off=20;
roff=(0:nrec-1)*plot_off+10;
plot_amp=15;
geolab=rec(1):rec(2);

figure
hold on
for i_rec=1:nrec
    plot(t,data(:,i_rec)*plot_amp+roff(i_rec),'-');
end
hold off
xlabel('Time (s)');
ylabel('Receivers');
set(gca,'ytick',roff,'yticklabel',geolab);
set(gca,'xlim',[trange(1) trange(2)],'ylim',[roff(1)-min(data(:,1))*plot_amp roff(nrec)+max(data(:,nrec))*plot_amp]);

% fk analysis starts here
%--------------------------

[m,n] = size(data); 


[x,y] = meshgrid(gausswin(n,1),gausswin(m,1));  
tracemask = x.*y; clear x; clear y;
data = tracemask.*data;  % smooth out some rough edges for ffting later, applied twice.

M = pow2(nextpow2(m)); 
N = pow2(nextpow2(n));

delt = t(2)-t(1);
fnyq = 1/(2*delt); 
delf = 2*fnyq/M;
freqs = -fnyq:delf:fnyq - delf;
delx = abs(offset(2)-offset(1)); 
%knyq = 2*pi/(2*delx); 
knyq = 1/(2*delx); 
delk = knyq/N;
ks = -knyq:delk:knyq-delk;

newdata = zeros(M,N);  

% moveout = 1000000000 %1850;  % Moveout velocity in m/s if want to flatten
% intercept =  0 % 0.008 % Time intercept for linear moveout% ind = round((intercept + offset/moveout)/delt)+1;
% for i = 1:n
% newdata(m/2+(1:m-ind(i)+1),n/2+i) = data(ind(i):m,i);
% end
% newdata(m/2+1:3*m/2,n/2+1:3*n/2) = tracemask.*newdata(m/2+1:3*m/2,n/2+1:3*n/2);
% imagesc(newdata); This section may be reincluded if one wants to preferentially 
% flatten on a given events as may be done in the vsp data.
% clear data ind
indxy = [(M-m)/2+1  m+(M-m)/2  (N-n)/2+1  n+(N-n)/2 ]; % Index corners of data within newdata
%newdata((M-m)/2+1:m+(M-m)/2,(N-n)/2+1:n+(N-n)/2) = data;
newdata(indxy(1):indxy(2),indxy(3):indxy(4)) = data;

clear data

z = fft2(newdata); z = fftshift(z); 
clear newdata
figure
%imagesc(ks,freqs(1:M/2+1),abs(z(1:M/2+1,:)));
imagesc(ks,fnyq+freqs(1:M/2+1),abs(flipud(z(1:M/2+1,:))));
max_fk = max(max(20*log10(abs(z))));
min_fk = min(min(20*log10(abs(z))));
colorbar
caxis ([min_fk max_fk]);
ylabel('frequency (Hz)'); xlabel('Wavenumber (radians/m)')
title('fk plot')

%select first region to filter out
[mask0,xi,yi] = roipoly;
%select second region to filter out
[mask1,xi,yi] = roipoly;

mask=mask0+mask1;
mask = flipud(mask);
mask = not(mask);  % change area with 1's to 0's and vice-versa
%mask(1:m+1,n+1:2*n) = 0*mask(1:m+1,n+1:2*n);  % This includes zero frequency
%mask(m+2:2*m,1:n+1) = 0*mask(m+2:2*m,1:n+1); % Also includes zero here.
mask(M/2+1:M,1:N) = fliplr(flipud(mask(1:M/2,1:N)));

kernel = fspecial('average',[30 3]);  %Smooth the edges of the filter
mask = conv2(single(mask),single(kernel),'same');

z = mask.*z;  
z = ifftshift(z); newdata = real(ifft2(z)); 
clear z
newdata = newdata(indxy(1):indxy(2),indxy(3):indxy(4))./tracemask;

figure
hold on
for i_rec=1:nrec
    plot(t,newdata(:,i_rec)*plot_amp+roff(i_rec),'-');
end
hold off
xlabel('Time (s)');
ylabel('Receivers');
set(gca,'ytick',roff,'yticklabel',geolab);
set(gca,'xlim',[trange(1) trange(2)],'ylim',[roff(1)-min(newdata(:,1))*plot_amp roff(nrec)+max(newdata(:,nrec))*plot_amp]);


[x,y] = meshgrid(gausswin(n,1),gausswin(m,1));  
tracemask = x.*y; clear x; clear y;
data = tracemask.*newdata;  % smooth out some rough edges for ffting later, applied twice.

M = pow2(nextpow2(m)); 
N = pow2(nextpow2(n));

delt = t(2)-t(1);
fnyq = 1/(2*delt); 
delf = 2*fnyq/M;
freqs = -fnyq:delf:fnyq - delf;
delx = abs(offset(2)-offset(1)); 
%knyq = 2*pi/(2*delx); 
knyq = 1/(2*delx); 
delk = knyq/N;
ks = -knyq:delk:knyq-delk;

newdata = zeros(M,N);  

% moveout = 1000000000 %1850;  % Moveout velocity in m/s if want to flatten
% intercept =  0 % 0.008 % Time intercept for linear moveout% ind = round((intercept + offset/moveout)/delt)+1;
% for i = 1:n
% newdata(m/2+(1:m-ind(i)+1),n/2+i) = data(ind(i):m,i);
% end
% newdata(m/2+1:3*m/2,n/2+1:3*n/2) = tracemask.*newdata(m/2+1:3*m/2,n/2+1:3*n/2);
% imagesc(newdata); This section may be reincluded if one wants to preferentially 
% flatten on a given events as may be done in the vsp data.
% clear data ind
indxy = [(M-m)/2+1  m+(M-m)/2  (N-n)/2+1  n+(N-n)/2 ]; % Index corners of data within newdata
%newdata((M-m)/2+1:m+(M-m)/2,(N-n)/2+1:n+(N-n)/2) = data;
newdata(indxy(1):indxy(2),indxy(3):indxy(4)) = data;

clear data

z = fft2(newdata); z = fftshift(z); 
clear newdata
figure
%imagesc(ks,freqs(1:M/2+1),abs(z(1:M/2+1,:)));
imagesc(ks,fnyq+freqs(1:M/2+1),abs(flipud(z(1:M/2+1,:))));
max_fk = max(max(20*log10(abs(z))));
min_fk = min(min(20*log10(abs(z))));
colorbar
caxis ([min_fk max_fk]);
ylabel('frequency (Hz)'); xlabel('Wavenumber (radians/m)')
title('New fk plot')




