function [rec,hdata] = process_mimo(fheader,DATA,process,compute,iloop)

if compute.njack > 0    
    % Throw other geophones' data
    rec_throw=iloop;
    for i_throw=1:inf
        [tf,ichan_throw]=ismember(rec_throw,[DATA.geonum]); % This gives only the maximum index therefore I need to do recursively

        if isempty(ichan_throw) || max(ichan_throw)==0
            break;
        end
        chan_throw=nonzeros(ichan_throw)';    

        DATA(chan_throw)=[];
    end
end

% Add noise
if process.noise>0
    % Find absolute maximum
    nchan=length(DATA);
    maxval=zeros(nchan,1);
    for i_chan=1:length(DATA)
        maxval(i_chan)=max(abs(DATA(i_chan).data));
    end
    absmax=max(maxval);
    noise_max = absmax*0.01*process.noise;       
            
    rec_id=unique([DATA.geonum]); % ID of live receiver
    for i_rec=1:length(rec_id)
        ichan=find([DATA.geonum]==rec_id(i_rec));
        nsamp=max([DATA(ichan).nsamp]);
        noise=(rand(1,nsamp)*2-1)'*noise_max;
        if compute.nmonte>0
            outnoise=strcat(fheader,'_noise_monte',num2str(iloop+1),'_geophone',num2str(rec_id(i_rec)));
        else
            outnoise=strcat(fheader,'_noise_geophone',num2str(rec_id(i_rec)));
        end
        save(outnoise,'noise'); 
        
        for i_chan=1:length(ichan)           
            DATA(ichan(i_chan)).data=DATA(ichan(i_chan)).data+noise;
        end
    end   
end
     
% Retain only the necessary component data
if isempty(strfind(process.comp,'E'))
    ichan_off=[DATA.comp]=='E';    
    DATA(ichan_off)=[];
end

if isempty(strfind(process.comp,'N'))
    ichan_off=[DATA.comp]=='N';
    DATA(ichan_off)=[];
end

if isempty(strfind(process.comp,'Z'))
    ichan_off=[DATA.comp]=='Z';
    DATA(ichan_off)=[];
end

if length(DATA)<1
    error('everything is thrown! check your process options.');
end

clear ichan_off
tdatum=min([DATA.t0]); % Always take minimum of origin time
nchan=length(DATA);
for i_chan=1:nchan
    DATA(i_chan).chnum=i_chan;
    DATA(i_chan).t0=DATA(i_chan).t0-tdatum;
end
      
dt_chan=[DATA.dt];
% Clip time if necessary
% Start time
if process.tclip(1)>0
    for i_chan=1:nchan
        if process.tclip(1)>DATA(i_chan).t0
            ind0=ceil((process.tclip(1)-DATA(i_chan).t0)/dt_chan(i_chan)); % This gives an index just before process.tclip(1)
            DATA(i_chan).data(1:ind0)=[];
            DATA(i_chan).t0=DATA(i_chan).t0+ind0*dt_chan(i_chan); % New origin shifts to the index ind0+1
            DATA(i_chan).nsamp=DATA(i_chan).nsamp-ind0;
        end
    end
end
% End time
tend_chan=([DATA.nsamp]-1).*[DATA.dt]+[DATA.t0];
if process.tclip(2)>0
    for i_chan=1:nchan
        if process.tclip(2)<tend_chan(i_chan)
            ind0=floor((tend_chan(i_chan)-process.tclip(2))/dt_chan(i_chan)); % This gives the (number of samples-1) after process.tclip(2)
            DATA(i_chan).data(end-ind0+1:end)=[];
            DATA(i_chan).nsamp=DATA(i_chan).nsamp-ind0;
        end
    end
end

% Filtering
if max(process.ffreq) > 0
    for i_chan = 1:nchan
        DATA(i_chan).data=buttwo(3,process.ffreq,[DATA(i_chan).dt],2,1,DATA(i_chan).data);
    end
end

%chan_comp=[DATA.comp];
tend_chan=([DATA.nsamp]-1).*[DATA.dt]+[DATA.t0];
t0_chan=[DATA.t0];

% Extract information only for active receivers
rec.id=unique([DATA.geonum]); % ID of live receiver
rec.nrec=length(rec.id);
rec.t0=zeros(rec.nrec,1);
rec.tend=zeros(rec.nrec,1);
%rec.chan=cell(rec.nrec,1);
%rec.comp=cell(rec.nrec,1);
for i_rec=1:rec.nrec
    irec=rec.id(i_rec);
    ichan=find([DATA.geonum]==irec); % Channel numbers corresponding to each live receiver
    rec.chan{i_rec}=ichan;
    rec.comp{i_rec}=[DATA(ichan).comp];
    rec.t0(i_rec)=max(t0_chan(ichan));
    rec.tend(i_rec)=max(tend_chan(ichan));
    rec.dt(i_rec)=max(dt_chan(ichan));
end
[rec_id,rec_chan]=unique([DATA.geonum]); % rec_chan gives the maximum channel index of each receiver
rec.x=reshape([DATA(rec_chan).xyz],3,rec.nrec)'; % reshape fill the elements column-wise
rec.tend=rec.tend-rec.t0; % Scale relative to each origin time

% Compute Hilbert transform
hdata=cell(1,nchan);
if compute.proj % Compute only the transform
    for i_chan=1:nchan    
        hdata{i_chan} = hilbert(DATA(i_chan).data);
    end
else % Compute also the envelop
    for i_chan=1:nchan    
        hdata{i_chan} = abs(hilbert(DATA(i_chan).data));
    end
end

clear DATA