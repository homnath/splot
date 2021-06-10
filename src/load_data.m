function [DATA]=load_data(rdata,varargin)
% Revision
%   HNG, Jan 25,2010: added seg2 file loading
%   HNG, Nov 11,2009

if nargin>1 && varargin{1}
    load('cache_file','DATA');
    return;
end

% Receivers information
if rdata.type==0 % Synthetic data
    
    nchan=rdata.nrec*3; % we will generate always 3C data in synthetic
    
    % initialize DATA structure    
    DATA=repmat(struct('chnum',0,'comp',' ','data',0.0,'dt',0.0, ...
    'geonum',0,'nsamp',0,'t0',0.0,'xyz',[NaN NaN NaN]),1,nchan);

    % Number of digits and extension format
    ndig=ceil(log10(rdata.nrec+1));%count_digit(nrec);    

    form=sprintf('%%0%dd',ndig);

    for i_rec=1:rdata.nrec
        count=sprintf(form,i_rec);
        fnamex=strcat(rdata.fheader,'.',count,'.x');
        fnamey=strcat(rdata.fheader,'.',count,'.y');
        fnamez=strcat(rdata.fheader,'.',count,'.z');
        if i_rec==1
            [datax,header] = readsac(fnamex);
            dt=header(1);
            t0=header(8);
            % Initilize arrays
            nsamp=size(datax,2);
            sx=zeros(nsamp,rdata.nrec);
            sy=zeros(nsamp,rdata.nrec);
            sz=zeros(nsamp,rdata.nrec);
            [datay] = readsac(fnamey);
            [dataz] = readsac(fnamez);
        else
            [datax] = readsac(fnamex);
            [datay] = readsac(fnamey);
            [dataz] = readsac(fnamez);
        end

        % Coversion to nm/s
        sx(:,i_rec)=datax(1,:);
        sy(:,i_rec)=datay(1,:);
        sz(:,i_rec)=-dataz(1,:); % Convert to Z upward from Z downward in e3d computation    
    end
    
    %[nrow,ncol]=size(sx);
    ncol=size(sx,2);

    if rdata.nrec ~= ncol
        error('number of receivers mismatched!');
    end 
    
    % Convert to structure    
    for i_rec=1:rdata.nrec
        for i_chan=1:3
            ichan=(i_rec-1)*3+i_chan;
            DATA(ichan).chnum=ichan;
            DATA(ichan).geonum=i_rec;
            DATA(ichan).t0=t0;
            DATA(ichan).dt=dt;
            DATA(ichan).nsamp=nsamp;
            if i_chan==1
                DATA(ichan).comp='E';
                DATA(ichan).data=sx(:,i_rec);                
            elseif i_chan==2
                DATA(ichan).comp='N';
                DATA(ichan).data=sy(:,i_rec);
            elseif i_chan==3
                DATA(ichan).comp='Z';
                DATA(ichan).data=sz(:,i_rec);
            else
                error('illegal i_chan: %d',i_chan);
            end
            DATA(ichan).xyz=[NaN NaN NaN];
        end
    end
elseif rdata.type==1 % MIMO data
    %DATA=[];
    load(rdata.mfile,'DATA');
    % remove unnecessary fields of DATA
    DATA=rmfield(DATA,{'filename','subset','flip','statname','plotmask','procmask','datestring'});
    % change Z up coordinates to down
    nchan=length(DATA);
    for i_chan=1:nchan
        DATA(i_chan).xyz(3)=-DATA(i_chan).xyz(3);
    end
    
elseif rdata.type==2 % SEG2 data
    % Michael Roth 09.05.03
    % HNG, 25.01.2010
    pivotyear=2000; % Reference year

    inpf_rset=fopen(rdata.rsetfile,'r');
    if inpf_rset<0
        error('file %s cannot be opened!',rdata.rsetfile);
    end
    fgetl(inpf_rset);
    bulk=textscan(inpf_rset,'%*d %d %s');
    geonum=bulk{1};
    comp=char(bulk{2});
    fclose(inpf_rset);
    
    %nchan=length(chnum);    

    %check for accessibility of datafile
    inpf_seg2=fopen(rdata.seg2file,'r');
    if (inpf_seg2==-1)
        error('file %s cannot be opened!',rdata.seg2file);
    end
    fclose(inpf_seg2);

    [rstat,fdb,tdb,data] = read_seg2(rdata.seg2file);
    % if there was a problem loading the data
    if rstat==0
        error('reading seg2 file %s failed!',rdata.seg2file);
    end
    try
      nchan=size(data,2);
      
      %get the date from the file descriptor block
      [~,rightstr]=strtok(fdb(strmatch('ACQUISITION_DATE',fdb),:));


      %%%%changes MRO 5.12.05
      ac_date=strtrim(rightstr); %fliplr(deblank(fliplr(deblank(rightstr))));
      i_test=isstrprop(ac_date,'alpha');
      if sum(i_test)==0
        % for format dd/mm/yy and dd/mm/yyyy
        dummy=ac_date(1:2);
        ac_date(1:2)=ac_date(4:5);
        ac_date(4:5)=dummy;
      %elseif sum(i_test)>0
        % for format dd/mmm/yy, dd/mmm/yyyy use the follwing line
      %  ac_date=ac_date ;
      end

      %get the time from the file descriptor block
      %and update the date vector dv (y,m,d,h,m,s)
      [~,rightstr]=strtok(fdb(strmatch('ACQUISITION_TIME',fdb),:));
      ac_time = rightstr;
      dv = datevec([ac_date ac_time]);

      %get the second fraction from the file descriptor block
      %and update the date vector dv (y,m,d,h,m,s.frac)
      [~,rightstr]=strtok(fdb(strmatch('ACQUISITION_SECOND_FRACTION',fdb),:));
      if ~isempty(rightstr)
        dv(6)=dv(6)+str2double(rightstr);
        sec_from_fdb=dv(6);
        %%%%changes MRO 5.12.05
      else
        sec_from_fdb=0;
      end
        
      % initialize DATA structure    
      DATA=repmat(struct('chnum',0,'comp',' ','data',0.0,'dt',0.0, ...
      'geonum',0,'nsamp',0,'t0',0.0,'xyz',[NaN NaN NaN]),1,nchan);

      for i=1:nchan
        [~,rightstr]=strtok(tdb{i}(strmatch('DELAY' ,tdb{i}),:));
        if ~isempty(rightstr)
          delay=str2double(rightstr);
          dv(6)=sec_from_fdb+delay;
        end
        [~,rightstr]=strtok(tdb{i}(strmatch('CHANNEL_NUMBER',tdb{i}),:));
        chnum=str2double(rightstr);
        DATA(i).chnum=chnum;
        DATA(i).geonum=geonum(chnum);
        DATA(i).t0=(datenum(dv)-datenum(pivotyear,1,1,0,0,0))*86400;
        [~,rightstr]=strtok(tdb{i}(strmatch('SAMPLE_INTERVAL' ,tdb{i}),:));
        DATA(i).dt=str2double(rightstr);
        % [leftstr,rightstr]=strtok(tdb{i}(strmatch('CHANNEL_NUMBER',tdb{i}),:));
        % chnum(i)=str2double(rightstr);
        % [leftstr,rightstr]=strtok(tdb{i}(strmatch('RECEIVER_LOCATION',tdb{i}),:));
        % rloc(i,:)=str2double(rightstr);
        % [leftstr,rightstr]=strtok(tdb{i}(strmatch('SAMPLE_INTERVAL' ,tdb{i}),:));
        % dtt(i)=str2double(rightstr);
        %DATA(i).filename=filename;
        %DATA(i).subset=1;
        %DATA(i).chnum=chnum(i);        
        %DATA(i).dt=dt;
        
        DATA(i).nsamp=length(data{i});
        %DATA(i).xyz=xyz{i};
        %DATA(i).flip=flip(i);
        DATA(i).comp=comp(DATA(i).chnum);
        %DATA(i).statname  = statname{i};
        %DATA(i).plotmask=plotmask(i);
        %DATA(i).procmask=procmask(i);
        DATA(i).data=data{i};
        DATA(i).xyz=[NaN NaN NaN];
        %tt=datenum(0,0,0,0,0,floor(DATA(i).t0))+datenum(pivotyear,1,1,0,0,0);
        %ttt= num2str(DATA(i).t0-floor(DATA(i).t0),'%0.4f');
        %DATA(i).datestring= strcat(datestr(tt),ttt(2:end));       
      end      
    catch exception
      str=strcat('error loading SEG2 data:',exception.message);
      disp(str);
      return
    end
elseif rdata.type==3 % Synthetic data produced by SPECFEM3D
    % Read station file
    inpf=fopen(rdata.sfile,'r');
    bulk=textscan(inpf,'%s %s %f %f %*f %f');
    fclose(inpf);
    rdata.nrec=length(bulk{1});
    if rdata.fheader(end) ~= '/' && rdata.fheader(end) ~= '\'
        rdata.fheader(end+1)='/';
    end    
    
    comp_read=['X','Y','Z','P']; % only for SPECFEM3D CARTESIAN
    comp=['E','N','Z','P']; % but we always process as ENZ components
    
    % Read only the two time values of first file
%     fnamex=strcat(rdata.fheader,char(bulk{1}(1)),'.',char(bulk{2}(1)),'.','FX',comp_read(1),'.',rdata.ext);
%     if ~exist(fnamex,'file')
%         error('file ''%s'' not found!',fnamex);
%     end
%     inpf=fopen(fnamex,'r');
%     t_info=fscanf(inpf,'%f %*f',[1 inf]);
%     fclose(inpf);
%     t0=t_info(1);
%     dt=t_info(2)-t_info(1);
%     nsamp=length(t_info);
%     clear t_info;
    initialized=false;
    ichan=0;
    for i_rec=1:rdata.nrec
        for i_chan=1:4
            ichan=ichan+1;
            fnamex=strcat(rdata.fheader,char(bulk{2}(i_rec)),'.',char(bulk{1}(i_rec)),'.','FX',comp_read(i_chan),'.',rdata.ext);
            if ~exist(fnamex,'file')
                warning('file ''%s'' not found!',fnamex);
                DATA (ichan).chnum=ichan;
                DATA(ichan).geonum=i_rec;            
                DATA(ichan).t0=[];
                DATA(ichan).dt=[];
                DATA(ichan).nsamp=[];
                DATA(ichan).comp=comp(i_chan);
                DATA(ichan).xyz=[];
                DATA(ichan).data=[];
                continue;
            end
            
            inpf=fopen(fnamex,'r');            
            %DATA(ichan).data=fscanf(inpf,'%*f %f',[1 inf])';
            % Extract t0, dt, and nsamp for the first time.
            if ~initialized      
                data_tmp=fscanf(inpf,'%f %f',[2 inf])';
                initialized = true;
                t0=data_tmp(1,1);
                dt=data_tmp(2,1)-data_tmp(1,1);
                nsamp=length(data_tmp);
                DATA(ichan).data=data_tmp(:,2);
            else
                DATA(ichan).data=fscanf(inpf,'%*f %f',[1 inf])';
            end
            DATA (ichan).chnum=ichan;
            DATA(ichan).geonum=i_rec;            
            DATA(ichan).t0=t0;
            DATA(ichan).dt=dt;
            DATA(ichan).nsamp=nsamp;
            DATA(ichan).comp=comp(i_chan);
            DATA(ichan).xyz=[bulk{3}(i_rec) bulk{4}(i_rec) bulk{5}(i_rec)];
            fclose(inpf);  
        end 
    end
    clear bulk;
elseif rdata.type==32 % Synthetic data produced by SPECFEM3D and take the difference
    % Read station file
    inpf=fopen(rdata.sfile,'r');
    bulk=textscan(inpf,'%s %s %f %f %*f %f');
    fclose(inpf);
    rdata.nrec=length(bulk{1});
    if rdata.fheader(end) ~= '/' && rdata.fheader(end) ~= '\'
        rdata.fheader(end+1)='/';
    end
    if rdata.fheader1(end) ~= '/' && rdata.fheader1(end) ~= '\'
        rdata.fheader1(end+1)='/';
    end
    
    comp_read=['X','Y','Z','P']; % only for SPECFEM3D CARTESIAN
    comp=['E','N','Z','P']; % but we always process as ENZ components
    
%     % Read only the two time values of first file
%     fnamex=strcat(rdata.fheader,char(bulk{1}(1)),'.',char(bulk{2}(1)),'.','FX',comp_read(1),'.',rdata.ext);
%     inpf=fopen(fnamex,'r');
%     t_info=fscanf(inpf,'%f %*f',[1 inf]);
%     fclose(inpf);
%     nsamp1=length(t_info);    
%     
%     fnamex1=strcat(rdata.fheader1,char(bulk{1}(1)),'.',char(bulk{2}(1)),'.','FX',comp_read(1),'.',rdata.ext);
%     inpf1=fopen(fnamex1,'r');
%     t_info1=fscanf(inpf1,'%f %*f',[1 inf]);
%     fclose(inpf1);
%     nsamp2=length(t_info1);
%     
%     nsamp=nsamp1;
%     if nsamp1 ~= nsamp2
%        fprintf(1,'WARNING: number of sampling points mismatched: %d, %d!\n',nsamp1,nsamp2);
%        fprintf(1,'some of the recordings will be truncated!\n');
%        nsamp=min([nsamp1, nsamp2]);
%     end
%     if t_info(1:nsamp) ~= t_info1(1:nsamp)
%        error('time vectors mismatched!'); 
%     end
%     t0=t_info(1);
%     dt=t_info(2)-t_info(1);  
%     
%     clear t_info t_info1;
    
    initialized=false;
    ichan=0;
    for i_rec=1:rdata.nrec
        for i_chan=1:4
            ichan=ichan+1;
            fnamex=strcat(rdata.fheader,char(bulk{2}(i_rec)),'.',char(bulk{1}(i_rec)),'.','FX',comp_read(i_chan),'.',rdata.ext);
            fnamex1=strcat(rdata.fheader1,char(bulk{2}(i_rec)),'.',char(bulk{1}(i_rec)),'.','FX',comp_read(i_chan),'.',rdata.ext);
            if ~exist(fnamex,'file')
                warning('file ''%s'' not found!',fnamex);
                DATA (ichan).chnum=[];
                DATA(ichan).geonum=[];            
                DATA(ichan).t0=[];
                DATA(ichan).dt=[];
                DATA(ichan).nsamp=[];
                DATA(ichan).comp=comp(i_chan);
                DATA(ichan).xyz=[];
                DATA(ichan).data=[];
                continue;
            end
            inpf=fopen(fnamex,'r'); 
            inpf1=fopen(fnamex1,'r');
            %DATA(ichan).data=fscanf(inpf,'%*f %f',[1 inf])';
            % Extract t0, dt, and nsamp for the first time.
            if ~initialized      
                data_tmp=fscanf(inpf,'%f %f',[2 inf])';
                data_tmp1=fscanf(inpf1,'%f %f',[2 inf])';
                initialized = true;
                t0=data_tmp(1,1);
                dt=data_tmp(2,1)-data_tmp(1,1);
                nsamp=length(data_tmp);
                DATA(ichan).data=data_tmp(:,2)-data_tmp1(:,2);
            else
                data_tmp=fscanf(inpf,'%*f %f',[1 inf])';
                data_tmp1=fscanf(inpf1,'%*f %f',[1 inf])';
                DATA(ichan).data=data_tmp-data_tmp1;
            end
            DATA (ichan).chnum=ichan;
            DATA(ichan).geonum=i_rec;            
            DATA(ichan).t0=t0;
            DATA(ichan).dt=dt;
            DATA(ichan).nsamp=nsamp;
            DATA(ichan).comp=comp(i_chan);
            DATA(ichan).xyz=[bulk{3}(i_rec) bulk{4}(i_rec) bulk{5}(i_rec)];
            fclose(inpf);
            fclose(inpf1);            
                 
%             fnamex1=strcat(rdata.fheader1,char(bulk{1}(i_rec)),'.',char(bulk{2}(i_rec)),'.','FX',comp_read(i_chan),'.',rdata.ext);
%             inpf=fopen(fnamex,'r');
%             inpf1=fopen(fnamex1,'r');
%             DATA(ichan).chnum=ichan;
%             DATA(ichan).geonum=i_rec;            
%             DATA(ichan).t0=t0;
%             DATA(ichan).dt=dt;
%             DATA(ichan).nsamp=nsamp;
%             DATA(ichan).comp=comp(i_chan);
%             % data recordings
%             data1=fscanf(inpf,'%*f %f',[1 inf])';
%             data2=fscanf(inpf1,'%*f %f',[1 inf])';
%             DATA(ichan).data=data1(1:nsamp)-data2(1:nsamp);            
%             DATA(ichan).xyz=[bulk{3}(i_rec) bulk{4}(i_rec) bulk{5}(i_rec)];
%             fclose(inpf);  
%             fclose(inpf1);
        end 
    end
    clear bulk;
elseif rdata.type==4 % SEGY data
     % Scan the header
    %[asc_h,bin_h,n_traces,nsamp,dt,dat_form,form_code,trace_length,extra_bytes] = segy_header(rdata.segyfile);
    [~,~,n_traces,nsamp,dt,dat_form,form_code,trace_length,extra_bytes] = segy_header(rdata.segyfile);
  %% Quick and dirty bug-fix for NGI-lab data where some header
  %% values are corrupted!!! VO: 31.7.2008
        if dt==0
          disp('error in file header, set sampling interval to 1e-7 sec! = 10MHz')
          dt=1e-7;
        end
  %% end of bug-fix
    nchan=n_traces;
    num_subsets=1; %n_traces/nchan;
    subset=1;
    if num_subsets~=fix(num_subsets)
      disp('routine load_data: wrong number of channels per subset or wrong number of traces in the SEGY file')
      return
    end
    
    if subset<1 || subset> num_subsets
      str=strcat('error invalid subset number:',num2str(subset),'max. number of sets:',num2str(num_subsets));
      disp(str)
      return;
    end

    try

      %get the tracenumbers to be extracted from the SEGY file
      trace_num = (1:nchan)+(subset-1)*nchan;

      %load the trace headers and the traces
      [trace_h,data] = read_segy(rdata.segyfile,trace_num,trace_length,extra_bytes,nsamp,form_code,dat_form);


      %rloc  = [trace_h(75,:)' trace_h(74,:)' trace_h(73,:)'];

      % initialize DATA structure    
      DATA=repmat(struct('chnum',0,'comp',' ','data',0.0,'dt',0.0, ...
      'geonum',0,'nsamp',0,'t0',0.0,'xyz',[NaN NaN NaN]),1,nchan);

      for i=1:nchan
        %wrap-around for absolute time
        %number of subset-1 times duration of traces
        %tabs=(subset-1)*dt*nsamp;
        
        year = trace_h(60,1);
        doy = trace_h(61,1);
        hour = trace_h(62,1);
        minute = trace_h(63,1);
        second = trace_h(64,1);
% time_mode_code = trace_h(65,1);
% 1 : 'Local';
% 2 : 'GMT';
% 3 : 'Other';
% 4 : 'UTC';
        [day, month] = doy2date(doy,year);
        t_absolute = datenum(year, month, day, hour, minute, second);
        %[f1 f2]=fileparts(ffilename);
        %DATA(i).filename=f2;
        %DATA(i).subset=subset;
        DATA(i).chnum=i;
        DATA(i).geonum=i;
  %% Quick and dirty bug-fix for NGI_lab data where some header
  %% values are corrupted!!! VO: 31.7.2009
        if dt==0
          disp('error in file header, set sampling interval to 0.1 microsec!')
          dt=0.0000001;
        end
  %% end of bug-fix
        DATA(i).dt=dt;
        %DATA(i).t0=tabs;
        DATA(i).nsamp=nsamp;
        DATA(i).xyz=[trace_h(22,i) trace_h(23,i) 0];
        %DATA(i).flip=1;
        DATA(i).comp='Z';
        %DATA(i).statname  = [];
        %DATA(i).plotmask=1;
        %DATA(i).procmask=1;
        DATA(i).data=data(:,i);
%        tt=datenum(0,0,0,0,0,floor(DATA(i).t0))+datenum(pivotyear,1,1,0,0,0);
%        ttt= num2str(DATA(i).t0-floor(DATA(i).t0),'%0.4f');
%        DATA(i).datestring= strcat(datestr(tt),ttt(2:end));
        DATA(i).t0 = (t_absolute - datenum(1980,1,1,0,0,0))*86400;
        %DATA(i).datestring = datestr(t_absolute);
      end
      %stat=1;
    catch exception
      str=strcat('error loading SEGY data:',exception.message);
      disp(str);
      return
    end
else
    error('wrong data type: %d!',rdata.type);
end

DATA=orderfields(DATA);

if ~isempty(rdata.rfile)
    % Receivers information
    infr=fopen(rdata.rfile);
    fgetl(infr);
    [recx]=fscanf(infr,'%f %f %f',[3 inf]);
    recx=recx'; % in (m)    
    fclose(infr);
    
    for i_rec=1:size(recx,1)
        ichan=find(i_rec==[DATA.geonum]);
        for i_chan=1:length(ichan)
            DATA(ichan(i_chan)).xyz=recx(i_rec,:);
        end        
    end
    fprintf('\nWARNING: existing receiver coordinates (if any) replaced with ''rfile'' coordinates!\n');
end
 
% [pathstr, name, ext] = fileparts(rdata.segyfile);
% ext=ext(2:length(ext));
% header_file=sprintf('%s_header',ext);
% out_header=fopen(header_file,'w');
% fprintf(out_header,'Total number of channels: %d\n',nchan-4);
% fprintf(out_header,'Channel name, x, y, z, start time, dt, nsamp, data file (32-bit floating point)\n');
% % Number of digits and extension format
% ndig=ceil(log10(nchan-4+1));%count_digit(nrec);
% form=sprintf('CHAN%%0%dd %%f %%f %%f %%f %%f %%d %%s\\n',ndig);
% form_data=sprintf('%s_CHAN%%0%dd_data',ext,ndig);
% ichan=0;
% for i_chan=5:nchan
%     ichan=ichan+1;
%     out_file=sprintf(form_data,ichan);
%     % write header file
%     %fprintf(out_header,'CHAN%
%     fprintf(out_header,form,ichan,DATA(i_chan).xyz,DATA(i_chan).t0,DATA(i_chan).dt,DATA(i_chan).nsamp,out_file);       
%         
%     % write data file    
%     out_data=fopen(out_file,'wb');
%     count=fwrite(out_data,DATA(i_chan).data,'float32');
%     if count ~= DATA(i_chan).nsamp
%         error('nsamp mismatched!');
%     end
%     fclose(out_data);
% end
% fclose(out_header);
save cache_file DATA


    
