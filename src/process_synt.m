function [rec,hdata]=process_synt(rdata,process,proj)

t_start=clock; % Record starting clock

kmtom=1000.0;

% Number of digits and extension format
ndig=ceil(log10(rdata.nrec+1));%count_digit(nrec);
form=sprintf('%%0%dd',ndig);

for i_rec=1:nrec
    count=sprintf(form,i_rec+offset);
    syntfnamex=strcat(syntfname,'.',count,'.x');
    syntfnamey=strcat(syntfname,'.',count,'.y');
    syntfnamez=strcat(syntfname,'.',count,'.z');
    if i_rec==1
        [datax,header] = readsac(syntfnamex);
        dt=header(1);
        % Initilize arrays
        nsamp=size(datax,2);
        sx=zeros(nsamp,nrec);
        sy=zeros(nsamp,nrec);
        sz=zeros(nsamp,nrec);
        [datay] = readsac(syntfnamey);
        [dataz] = readsac(syntfnamez);
    else
        [datax] = readsac(syntfnamex);
        [datay] = readsac(syntfnamey);
        [dataz] = readsac(syntfnamez);
    end

    % Coversion to nm/s
    sx(:,i_rec)=datax(1,:)*mtonano;
    sy(:,i_rec)=datay(1,:)*mtonano;
    sz(:,i_rec)=-dataz(1,:)*mtonano; % Convert to Z upward from Z downward in e3d computation    
end

[nrow,ncol]=size(sx);

if nrec ~= ncol
    error('Number of receivers mismatched!');
end

% Add synthetic (random) noise on data
% Change here to bandpass filter the data 
% freq(1)=100; freq(2)=400;

% Noise is same for all geophones
if noise_percent>0
    noise_max = max(max(max(max(abs(sx), abs(sy)), abs(sz))))*0.01*noise_percent; 
    
    
    for i_rec=1:nrec
        noise=(rand(1,nsamp)*2-1)*noise_max;
        
        outnoise=strcat(fheader,'_noise','_monte',num2str(nmonte),'_geophone',num2str(i_rec));            
        save(outnoise,'noise');            
                
%         if njack==1
%             disp('Wrong');
%             outnoise=strcat(fheader,'_noise_geophone',num2str(i_rec));            
%             save(outnoise,'noise');            
%         elseif njack>1
%             % Read existing noise written for njack=1
%             disp('Wrong');
%             inpnoise=strcat(fheader,'_noise_geophone',num2str(i_rec),'.mat');
%             load(inpnoise);            
%         end
        %load(outnoise);
            
        sx(:,i_rec)=sx(:,i_rec)+noise';
        sy(:,i_rec)=sy(:,i_rec)+noise';
        sz(:,i_rec)=sz(:,i_rec)+noise';        
    end
end
clear noise;

% Filter
if max(freq)>0
    sx=buttwo(3,freq,dt,2,1,sx);
    sy=buttwo(3,freq,dt,2,1,sy);
    sz=buttwo(3,freq,dt,2,1,sz);
end

% Extract information only for active receivers
rec.id=1:nrec;
rec.t0=0;
rec.dt=dt;
rec.tend=(nsamp-1)*dt;

% Compute Hilbert transform
if proj % Compute only the transform
    hdatax = hilbert(sx);
    hdatay = hilbert(sy);
    hdataz = hilbert(sz);
else % Compute also the envelop
    hdatax = abs(hilbert(sx));
    hdatay = abs(hilbert(sy));
    hdataz = abs(hilbert(sz));
end 

% Receivers information
% Receivers information
infr=fopen(rdata.rfile);
fgetl(infr);
[recx]=fscanf(infr,'%f %f %f',[3 inf]);
recx=recx';
rec.x=kmtom*recx(rec_id,:); % in (m)    
fclose(infr);