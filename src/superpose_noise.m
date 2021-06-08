function superpose_noise(dheader,nheader,nrec,mg)

% sensitivity of the noise sensors
nsens=22; %sensitivity of the noise sensors [V/m/s]
nchan=nrec*1; % we will use only vertical component
    
% initialize DATA structure    
DATA=repmat(struct('chnum',0,'comp',' ','data',0.0,'dt',0.0, ...
'geonum',0,'nsamp',0,'t0',0.0,'xyz',[NaN NaN NaN],'filename',[],'subset',[],'flip',[],'statname',[],'plotmask',[],'procmask',[],'datestring',[]),1,nchan);

% Number of digits and extension format
ndig=ceil(log10(nrec+1));%count_digit(nrec);    

form=sprintf('%%0%dd',ndig);

for i_rec=1:nrec
    count=sprintf(form,i_rec);
    %fnamex=strcat(rdata.fheader,'.',count,'.x');
    %fnamey=strcat(rdata.fheader,'.',count,'.y');
    fnamez=strcat(dheader,'.',count,'.z');
    if i_rec==1
        [dataz,header] = readsac(fnamez);
        dt=header(1);
        t0=header(8);
        % Initilize arrays
        nsamp=size(dataz,2);
        %sx=zeros(nsamp,rdata.nrec);
        %sy=zeros(nsamp,rdata.nrec);
        sz=zeros(nsamp,nrec);
        %[datay] = readsac(fnamey);
        %[dataz] = readsac(fnamez);
    else
        %[datax] = readsac(fnamex);
        %[datay] = readsac(fnamey);
        [dataz] = readsac(fnamez);
    end

    % Coversion to nm/s
    %sx(:,i_rec)=datax(1,:);
    %sy(:,i_rec)=datay(1,:);
    sz(:,i_rec)=-dataz(1,:); % Convert to Z upward from Z downward in e3d computation    
end
dt=4e-4; %only for this case 
%cross=zeros(nrec,nrec);
% test cross correlation
%for i_rec=1:nrec
%    for j_rec=1:nrec
%        cross(i_rec,j_rec)=sum(xcorr(sz(:,i_rec),sz(:,j_rec),'coeff'));
%    end
%end
%figure
%imagesc(cross);

% Open noise
inhead=fopen(nheader,'r');
fgetl(inhead);
fgetl(inhead);
tag=textscan(inhead,'%*s %*f %*f %*f %*f %f %d %s\n',1);
dt_noise=tag{1};
nsamp_noise=tag{2};

dpath=fileparts(nheader);    

file_noise=strcat(dpath,'/',tag{3});

nsamp_noise=1+ceil((nsamp-1)*dt/dt_noise);
sz_noise=zeros(nsamp_noise,nrec);

indata=fopen(file_noise{1},'r');
sz_noise(1:nsamp_noise,1)=fscanf(indata,'%f',nsamp_noise);
fclose(indata);


for i_rec=2:nrec
  tag=textscan(inhead,'%*s %*f %*f %*f %*f %*f %*d %s\n',1);
  file_noise=strcat(dpath,'/',tag{1});
  indata=fopen(file_noise{1},'r');
  sz_noise(1:nsamp_noise,i_rec)=fscanf(indata,'%f',nsamp_noise);
  fclose(indata);
end
% correct the noise data with sensitivity
sz_noise=sz_noise/nsens;
%max_data=max(abs(sz));
%max_noise=max(abs(sz_noise));

%mfac=max_noise./(nsr*max_data);

%fprintf(1,'Expected NSR: %f, Multiplication factor: %f\n',nsr,mfac);
% original magnitude
m0=1e15;
mg0=2*log10(m0)/3-10.7;
m1=10^(1.5*(mg+10.7)); % Seismic moment for the required magnitude
mfac=m1/m0;

fprintf(1,'Original magnitude: %f\n',mg0);
fprintf(1,'Desied magnitude: %f\n',mg);
fprintf(1,'Multiplication factor: %f\n',mfac);


% resample synthetic data
[p,q]=rat(dt/dt_noise,0.00001);
sz=mfac*resample(sz,p,q);

nsamp=size(sz,1);

maxsz=max(abs(sz));
maxsz_noise=max(abs(sz_noise));
figure
plot(1:nrec,maxsz,'-og');
hold on
plot(1:nrec,maxsz_noise,'-or');
legend('Max signal', 'Max noise');

figure
plot(1:nrec,maxsz_noise./maxsz,'-og')

mean_signal=mean(abs(sz));
mean_noise=mean(abs(sz_noise));
figure
plot(1:nrec,mean_signal,'-og');
hold on
plot(1:nrec,mean_noise,'-or');
legend('Mean signal', 'Mean noise');

figure
plot(1:nrec,mean_noise./mean_signal,'-og')


% superpose
for i_rec=1:nrec
    sz_noise(1:nsamp,i_rec)=sz_noise(1:nsamp,i_rec)+sz(1:nsamp,i_rec);
end

% Convert to structure    
for i_rec=1:nrec
    for i_chan=1:1
        ichan=(i_rec-1)*1+i_chan;
        DATA(ichan).chnum=ichan;
        DATA(ichan).geonum=i_rec;
        DATA(ichan).t0=t0;
        DATA(ichan).dt=dt_noise;
        DATA(ichan).nsamp=nsamp_noise;
        %if i_chan==1
        %    DATA(ichan).comp='E';
        %    DATA(ichan).data=sx(:,i_rec);                
        %elseif i_chan==2
        %    DATA(ichan).comp='N';
        %    DATA(ichan).data=sy(:,i_rec);
        %elseif i_chan==3
            DATA(ichan).comp='Z';
            DATA(ichan).data=sz_noise(:,i_rec);
        %else
        %    error('illegal i_chan: %d',i_chan);
        %end
        DATA(ichan).xyz=[NaN NaN NaN];
    end
end
mfile=strcat(dpath,'/','superpose_data.mat');
save(mfile,'DATA');