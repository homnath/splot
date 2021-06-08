function segy2for(segyfile)

     % Scan the header
    [asc_h,bin_h,n_traces,nsamp,dt,dat_form,form_code,trace_length,extra_bytes] = segy_header(segyfile);
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
    if num_subsets~=fix(num_subsets);
      disp('routine load_data: wrong number of channels per subset or wrong number of traces in the SEGY file')
      return
    end
    
    if (subset<1 || subset> num_subsets);
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
    out_header=fopen('out_header','w');

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
        % write header file
        
        
        % write data file
        out_file=sprintf('data_chan%d',i);
        out_data=fopen(out_file,'w');
        count=fwrite(out_data,data(:,i),'float32');
      end
      %stat=1;
    catch
      str=strcat('error loading SEGY data');
      disp(str);
      return
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
save cache_file DATA


    
