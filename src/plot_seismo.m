function plot_seismo(preinfo,model,process,seisplot,DATA)
% This program plots seismograms, envelopes, and SNR profiles with and without
% onsets, and polarization analysis. See fo example, self-explanatory input 

%last revision: HNG Apr 29, 2010

% Code below is only for plotting the figure(fig_prop1,fig_val1)
fig_prop1='visible';
fig_val1='on';
if preinfo.plot==0
    fig_val1='off';
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
    %nsamp=max([DATA(1).nsamp]);
    %noise=(rand(1,nsamp)*2-1)'*noise_max;        
    rec_id=unique([DATA.geonum]); % ID of live receiver
    for i_rec=1:length(rec_id)
        ichan=find([DATA.geonum]==rec_id(i_rec));
        nsamp=max([DATA(ichan).nsamp]);
        noise=(rand(1,nsamp)*2-1)'*noise_max;        
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

if isempty(strfind(process.comp,'P'))
    ichan_off=[DATA.comp]=='P';
    DATA(ichan_off)=[];
end

if length(DATA)<1
    error('everything is thrown! check your process options.');
end

clear ichan_off
tdatum=min([DATA.t0]); % Always take minimum of origin time
nchan=length(DATA);
geolab=cell(nchan,1);
for i_chan=1:nchan
    DATA(i_chan).chnum=i_chan;
    DATA(i_chan).t0=DATA(i_chan).t0-tdatum;
    geolab{i_chan}=strcat(num2str(DATA(i_chan).geonum),DATA(i_chan).comp);
end

n0=10;
dn=20;
rr_comb=(0:nchan-1)'*dn+n0; % Plot bottom to top
     
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

clear tend_chan

% % dirty tricks only for Åknes data
% allind=1:DATA(1).nsamp;
% corind=allind>2.1/DATA(1).dt;
% 
% DATA(1).data(corind)=0.0;
% DATA(2).data(corind)=0.0;
% DATA(3).data(corind)=0.0;
% DATA(10).data(corind)=0.0;
% DATA(11).data(corind)=0.0;
% DATA(12).data(corind)=0.0;
%% dirty tricks


[rec_id,rec_chan]=unique([DATA.geonum]); % rec_chan gives the maximum channel index of each receiver
nrec=length(rec_id);

%optional=[]; % Optional choices

tp=[]; ts=[]; % tp and ts are empty by default
recx=reshape([DATA(rec_chan).xyz],3,nrec)';
% Compute arrival times
if ~isempty(model) 
    %[rec_id,rec_chan]=unique([DATA.geonum]); % rec_chan gives the maximum channel index of each receiver
    %nrec=length(rec_id);
    %recx=reshape([DATA(rec_chan).xyz],3,nrec)'; % reshape fill the elements column-wise   
        
%     if max(recx(:,3))<0
%         recx(:,3)=-recx(:,3);
%     end
%     fprintf(1,'WARNING: signs of receiver Z-cordinates are changed\n');
    
    % compute travel times and basis vectors only for live receivers
    [tp,ts,pbasis,shbasis,svbasis]=timebasis_vect(model,rec_id,recx,model.src(1:3));
    tp=tp+model.src(4); ts=ts+model.src(4); % Scaled to each seismogram origin time
    %optional.tp=tp;
    %optional.ts=ts;
end


% Separate components
% Code below should be changed for mixed component data. In that case
% intializing the data with equal number of receivers should do the job.
datax=repmat(struct('data',0.0,'chnum',0,'geonum',0, ...
    'dt',0.0,'t0',0.0,'nsamp',1,'xyz',[NaN NaN NaN],'comp','E'),1,nrec);
datay=repmat(struct('data',0.0,'chnum',0,'geonum',0, ...
    'dt',0.0,'t0',0.0,'nsamp',1,'xyz',[NaN NaN NaN],'comp','N'),1,nrec);
dataz=repmat(struct('data',0.0,'chnum',0,'geonum',0, ...
    'dt',0.0,'t0',0.0,'nsamp',1,'xyz',[NaN NaN NaN],'comp','Z'),1,nrec);

datax=orderfields(datax);
datay=orderfields(datay);
dataz=orderfields(dataz);

% set appropriate DATA structures
ichan=[DATA.comp]=='E';
irec=[DATA(ichan).geonum];
[tf,irec_process]=ismember(irec,rec_id); % For all 3 components data only datax(1:nrec)=DATA(ichan); is enough
datax(irec_process)=DATA(ichan);
rrx=rr_comb(ichan);


ichan=[DATA.comp]=='N';
irec=[DATA(ichan).geonum];
[tf,irec_process]=ismember(irec,rec_id); % For all 3 components data only datay(1:nrec)=DATA(ichan); is enough
datay(irec_process)=DATA(ichan);
rry=rr_comb(ichan);

ichan=[DATA.comp]=='Z';
irec=[DATA(ichan).geonum];
[tf,irec_process]=ismember(irec,rec_id); % For all 3 components data only dataz(1:nrec)=DATA(ichan); is enough
dataz(irec_process)=DATA(ichan);
rrz=rr_comb(ichan);

% set to appropriate geophone number if any
% this avoids the error during plotting some components of a geophone are missing
 for i_rec=1:nrec
     geonum_correct=max([datax(i_rec).geonum datay(i_rec).geonum dataz(i_rec).geonum]);
     datax(i_rec).geonum=geonum_correct;
     datay(i_rec).geonum=geonum_correct;
     dataz(i_rec).geonum=geonum_correct;
 end
 
% for i=1:nchan
%     xcrec(i)=DATA(i).xyz(1);
%     ycrec(i)=DATA(i).xyz(2);
% end
% figure;
% plot(ycrec,xcrec,'o');

% for NGI observed data
%datay(10).nsamp=1;
%datay(10).data=0;

clear DATA ichan irec irec_process tf rec_process geonum_correct
% Denosing
% fprintf(1,'denoising the signals...');
% for i_rec=1:nrec
%     fprintf(1,'%d\n',i_rec);
%     if norm(datax(i_rec).data)~=0
%         datax(i_rec).data=denoise(datax(i_rec).data,datax(i_rec).data(1:50),10);
%     end
%     if norm(datay(i_rec).data)~=0
%         datay(i_rec).data=denoise(datay(i_rec).data,datay(i_rec).data(1:50),10);
%     end
%     if norm(dataz(i_rec).data)~=0
%         dataz(i_rec).data=denoise(dataz(i_rec).data,dataz(i_rec).data(1:50),10);
%     end
% end
% fprintf(1,'complete!\n');

% Frequency filter
if max(process.ffreq)>0
    for i_rec=1:nrec
        if norm(datax(i_rec).data)~=0
            datax(i_rec).data=buttwo(process.forder,process.ffreq,datax(i_rec).dt,process.ftype,1,datax(i_rec).data);
        end
        if norm(datay(i_rec).data)~=0
            datay(i_rec).data=buttwo(process.forder,process.ffreq,datay(i_rec).dt,process.ftype,1,datay(i_rec).data);
        end
        if norm(dataz(i_rec).data)~=0
            dataz(i_rec).data=buttwo(process.forder,process.ffreq,dataz(i_rec).dt,process.ftype,1,dataz(i_rec).data);
        end
    end
end

% Error prediction filtering
if max(process.erpf) > 0
    for i_rec = 1:nrec
        if norm(datax(i_rec).data)~=0
            lwnoise=round(process.tnoise/datax(i_rec).dt);
            ncoef=min(process.ncoef, lwnoise-1);
            a=lpc(datax(i_rec).data(1:lwnoise),ncoef); %[a,g]=lpc(datax(i_rec).data(1:lwnoise),ncoef);
            prederr=filter([1 real(a(2:ncoef))],1,datax(i_rec).data);
            prederr(1:ncoef-1)=prederr(ncoef);
            datax(i_rec).data=prederr;
        end
        if norm(datay(i_rec).data)~=0
            lwnoise=round(process.tnoise/datay(i_rec).dt);
            ncoef=min(process.ncoef, lwnoise-1);
            a=lpc(datay(i_rec).data(1:lwnoise),ncoef); %[a,g]=lpc(datay(i_rec).data(1:lwnoise),ncoef);
            prederr=filter([1 real(a(2:ncoef))],1,datay(i_rec).data);
            prederr(1:ncoef-1)=prederr(ncoef);
            datay(i_rec).data=prederr;
        end
        if norm(dataz(i_rec).data)~=0
            lwnoise=round(process.tnoise/dataz(i_rec).dt);
            ncoef=min(process.ncoef, lwnoise-1);
            a=lpc(dataz(i_rec).data(1:lwnoise),ncoef); %[a,g]=lpc(dataz(i_rec).data(1:lwnoise),ncoef);
            prederr=filter([1 real(a(2:ncoef))],1,dataz(i_rec).data);
            prederr(1:ncoef-1)=prederr(ncoef);
            dataz(i_rec).data=prederr;
        end
    end
end

% Original data
datax_old=datax;
datay_old=datay;
dataz_old=dataz;

% % cross spectrum matrix
% for i_rec=1:nrec
%     fft_data(:,i_rec)=fft(datax(i_rec).data);
% end
% f_ind=100;
% for i=1:length(fft_data(:,1))
% for i_rec=1:nrec
%     for j_rec=1:nrec
%         c(i,i_rec,j_rec)=fft_data(i,i_rec)*conj(fft_data(i,j_rec));
%     end
% end
% end
% %figure;
% %imagesc(abs(c));

% % temporary error
% for i_rec=1:nrec
%     mean_data(i_rec)=mean(abs(dataz(i_rec).data));
%     max_data(i_rec)=max(abs(dataz(i_rec).data));
% end
% figure
% plot(1:nrec,mean_data,'-ok');
% hold on
% plot(1:nrec,max_data,'-or');
% set(gca,'xtick',1:nrec,'xticklabel',geolab);

nplot=length(seisplot);

nptype=zeros(6,1); %(currently we have type 0, 1, 2, 3, 4, 5)
% Start plotting
for i_plot=1:nplot
    % count plot types
    nptype(seisplot(i_plot).type+1)=nptype(seisplot(i_plot).type+1)+1;
    
    scount=strcat('_n',int2str(nptype(seisplot(i_plot).type+1))); % string of a counter
    if seisplot(i_plot).type<=4
        superpose=seisplot(i_plot).optional.superpose;
        cmap=seisplot(i_plot).par.cmap;
        norm_fact=seisplot(i_plot).par.norm;
        plot_ampl=seisplot(i_plot).par.amp;
        clip_ampl=seisplot(i_plot).par.clipamp;        
        
        %rr=(0:nrec-1)*dn+n0; % Only for equal spacing plot
        geonum_active=[datax.geonum];
        rr=(double(geonum_active)-1.0)*dn+n0; % Only for equal spacing plot
        %rrx=([dataz_old.chnum]-1)*dn+n0; % Plot E on top
        %rry=([datay_old.chnum]-1)*dn+n0;
        %rrz=([datax_old.chnum]-1)*dn+n0; % Plot Z on bottom
        ylim=[min(rr_comb)-clip_ampl max(rr_comb)+clip_ampl];        
        
        optional=seisplot(i_plot).optional;
        optional.tp=tp;
        optional.ts=ts;
        
        % seismogram colormap is plotted only for separate components
        if cmap==1
            optional.superpose=0;
            superpose=0;
        end
    end
    % Plot traces
    if seisplot(i_plot).type==0    
        fprintf(1,'plotting original traces...');
        geomax=zeros(1,nrec);
        if norm_fact==2            
            % Geophone maximum    
            for i_rec=1:nrec
                geomax(i_rec)=max([max(abs(datax_old(i_rec).data)) max(abs(datay_old(i_rec).data)) max(abs(dataz_old(i_rec).data))]);
            end
        end
        
%         % ---------Plot ordered data-------------
%         for i_rec=1:nrec
%             src2g(i_rec)=sqrt(sum((recx(i_rec,:)-process.sort_src').^2));
%         end
%         [src2g,ind]=sort(src2g);
%         %optional.onset=0;
%         if strfind(process.comp,'E')>0
%             figure(fig_prop1,fig_val1)
%             plot_data(datax_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%             %title(['East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%             if preinfo.save>0
%                 outf=strcat(preinfo.fighead,'_original_trace_ordered_east_%d',nplot);
%                 savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%             end
%         end    
% 
%         if strfind(process.comp,'N')>0
%             figure(fig_prop1,fig_val1)
%             plot_data(datay_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%             %title(['North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%             if preinfo.save>0
%                 outf=strcat(preinfo.fighead,'_original_trace_ordered_north_%d',nplot);
%                 savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%             end
%         end
% 
%         if strfind(process.comp,'Z')>0
%             figure(fig_prop1,fig_val1)
%             plot_data(dataz_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%             %title(['Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%             if preinfo.save>0
%                 outf=strcat(preinfo.fighead,'_original_trace_ordered_vertical_%d',nplot);
%                 savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%             end
%         end
%         
%         figure(fig_prop1,fig_val1)
%         hold on
%         if strfind(process.comp,'E')>0
%             plot_data(datax_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%         end
%         if strfind(process.comp,'N')>0
%             plot_data(datay_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)
%         end
%         if strfind(process.comp,'Z')>0
%             plot_data(dataz_old(ind),src2g,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)
%         end
%         title(['East (black), North (red) and Z (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%         if preinfo.save>0
%             outf=strcat(preinfo.fighead,'_ordered_original_traces');
%             savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%         end
%         % --------Plot ordered data-------------

%         % Master segment
%         %i1=ceil((2.13-dataz_old(3).t0)/dataz_old(3).dt);
%         %i2=ceil((2.203-dataz_old(3).t0)/dataz_old(3).dt);
%         i1=ceil((2.21-datax_old(3).t0)/datax_old(3).dt);
%         i2=ceil((2.309-datax_old(3).t0)/datax_old(3).dt);
%         
%         nsamp_cor=i2-i1+1;
%         seg_master=dataz_old(3).data(i1:i2);
%         
%         % Original data
%         datax_cor=datax_old;
%         datay_cor=datay_old;
%         dataz_cor=dataz_old;
%         seg_cor=zeros(nsamp_cor,1);
%         for i_rec=1:nrec
%             for i_samp=1:datax_old(i_rec).nsamp
%                 if i_samp+nsamp_cor-1>datax_old(i_rec).nsamp
%                     i2=datax_old(i_rec).nsamp;
%                 else
%                     i2=i_samp+nsamp_cor-1;
%                 end
%                 nsamp_use=i2-i_samp+1;
%                 seg_cor(1:nsamp_use,1)=datax_old(i_rec).data(i_samp:i2);
%                 datax_cor(i_rec).data(i_samp)=xcorr(seg_master,seg_cor,0,'coeff');
%                 seg_cor(1:nsamp_use,1)=datay_old(i_rec).data(i_samp:i2);
%                 datay_cor(i_rec).data(i_samp)=xcorr(seg_master,seg_cor,0,'coeff');
%                 seg_cor(1:nsamp_use,1)=dataz_old(i_rec).data(i_samp:i2);
%                 dataz_cor(i_rec).data(i_samp)=xcorr(seg_master,seg_cor,0,'coeff');
%             end
%             datax_cor(i_rec).data=abs(hilbert(datax_cor(i_rec).data));
%             datay_cor(i_rec).data=abs(hilbert(datay_cor(i_rec).data));
%             dataz_cor(i_rec).data=abs(hilbert(dataz_cor(i_rec).data));
%             %figure;plot(data_cor);
%         end
%         %if superpose==0
%             if strfind(process.comp,'E')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(datax_cor,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_correlated_trace_east_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end    
% 
%             if strfind(process.comp,'N')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(datay_cor,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_correlated_trace_north_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end
% 
%             if strfind(process.comp,'Z')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(dataz_cor,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_correlated_trace_vertical_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end
%         %end
                
            

%         if superpose==0
%             if strfind(process.comp,'E')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(datax_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_original_trace_east_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end    
% 
%             if strfind(process.comp,'N')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(datay_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_original_trace_north_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end
% 
%             if strfind(process.comp,'Z')>0
%                 figure(fig_prop1,fig_val1)
%                 plot_data(dataz_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
%                 title(['Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
%                 if preinfo.save>0
%                     outf=strcat(preinfo.fighead,'_original_trace_vertical_%d',nplot);
%                     savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
%                 end
%             end
%         end

        if cmap==1
            % plot colormap
            if strfind(process.comp,'E')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(datax_old,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['East component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_map_east',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end    

            if strfind(process.comp,'N')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(datay_old,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['North component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_map_north',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'Z')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(dataz_old,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['Z component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_map_vertical',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end
            fprintf(1,'complete!\n');
            continue;
        end


        if superpose==0
            if strfind(process.comp,'E')>0
                figure(fig_prop1,fig_val1)
                plot_data(datax_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                %title(['East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_east',scount);
                    %print -depsc -r300 ~/publication/draft/seg2011/figure/py_myz_noise40_east.eps
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end    

            if strfind(process.comp,'N')>0
                figure(fig_prop1,fig_val1)
                %datay_old(10).nsamp=1; datay_old(10).data=0.; %only for plot_py_2005_355_00_27_26_4.in
                plot_data(datay_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                %title(['North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_north',scount);
                    %print -depsc -r300 ~/publication/draft/seg2011/figure/py_myz_noise40_north.eps
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'Z')>0
                figure(fig_prop1,fig_val1)
                plot_data(dataz_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                %title(['Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_original_trace_vertical',scount);
                    %print -depsc -r300 ~/publication/draft/seg2011/figure/py_myz_noise40_vert.eps
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end
        elseif superpose==1
            figure(fig_prop1,fig_val1)
            hold on
            if strfind(process.comp,'E')>0
                plot_data(datax_old,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay_old,rr,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz_old,rr,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)
            end
            title(['East (black), North (red) and Z (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_original_traces',scount);
                %print -depsc -r300 aaknes_synt.eps
                %print -depsc -r300 ~/publication/draft/seg2011/figure/py_myz_noise40_all.eps
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        elseif superpose==2            
            figure(fig_prop1,fig_val1)
            hold all
            if strfind(process.comp,'E')>0
                plot_data(datax_old,rrx,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay_old,rry,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz_old,rrz,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            ylabel('Channels');
            set(gca,'box','on','ylim',ylim);
            set(gca,'YTick',rr_comb,'YTickLabel',geolab,'FontSize',14);
            %title(['Original traces, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_original_traces',scount);
                %set(gcf,'position',[1,1,540,850]);
                %set(gcf,'PaperPositionMode','auto');
                %print -depsc -r300 ~/publication/draft/seg2011/figure/test.eps
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        else
            error('wrong superpose option %d! possible values are (0, 1, 2)', superpose);
        end
        % ******************End of traces ********************************
        fprintf(1,'complete!\n');
    elseif seisplot(i_plot).type==1
        % Plot rotated signals
        fprintf(1,'plotting rotated signals...');
        
        datap=datax;
        datash=datax;
        datasv=datax;

%         for i_rec=1:nrec
%             datap(i_rec).data =[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*pbasis(i_rec,:)';
%             datash(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*shbasis(i_rec,:)';
%             datasv(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*svbasis(i_rec,:)';        
%         end
        
        % Rotation
        for i_rec=1:nrec
            %datap(i_rec).data =[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*pbasis(i_rec,:)';
            %datash(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*shbasis(i_rec,:)';
            %datasv(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*svbasis(i_rec,:)';
            datap(i_rec).data =(datax_old(i_rec).data*pbasis(i_rec,1))+(datay_old(i_rec).data*pbasis(i_rec,2))+(dataz_old(i_rec).data*pbasis(i_rec,3));
            datash(i_rec).data=(datax_old(i_rec).data*shbasis(i_rec,1))+(datay_old(i_rec).data*shbasis(i_rec,2))+(dataz_old(i_rec).data*shbasis(i_rec,3));
            datasv(i_rec).data=(datax_old(i_rec).data*svbasis(i_rec,1))+(datay_old(i_rec).data*svbasis(i_rec,2))+(dataz_old(i_rec).data*svbasis(i_rec,3));
        end
        
        geomax=zeros(1,nrec);
        if norm_fact==2
            % Geophone maximum
            for i_rec=1:nrec
                geomax(i_rec)=max([max(abs(datap(i_rec).data)) max(abs(datash(i_rec).data)) max(abs(datasv(i_rec).data))]);
            end
        end
        
        if cmap==1
            figure(fig_prop1,fig_val1)
            plot_seismap(datap,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated signals map -- P, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_map_p',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end

            figure(fig_prop1,fig_val1)
            plot_seismap(datash,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated signals map -- SH, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_map_sh',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end

            figure(fig_prop1,fig_val1)
            plot_seismap(datasv,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated signals map -- SV, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_map_sv',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            fprintf(1,'complete!\n');            
            continue;
        end
        
        if optional.area>0
            % Compute vectors for stacked regions
            [tvectp,fvectp] = vect4bound({datap.data},[datap.t0]',[datap.dt]',optional.tp,optional.tp+optional.twin(1));
            [tvects,fvects] = vect4bound({datash.data},[datash.t0]',[datash.dt]',optional.ts,optional.ts+optional.twin(2));
        
            % Set optional parameters
            optional.tvectp=tvectp;
            optional.tvects=tvects;
            optional.fvectp=fvectp;
            optional.fvects=fvects;
        end

        if superpose==0
            figure(fig_prop1,fig_val1)
            plot_data(datap,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated signals -- P, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_p',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end

            figure(fig_prop1,fig_val1)
            plot_data(datash,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated signals -- SH, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_sh',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end

            figure(fig_prop1,fig_val1)
            plot_data(datasv,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated signals -- SV, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces_sv',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        elseif superpose==1
            figure(fig_prop1,fig_val1)
            hold on
            plot_data(datap,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)

            plot_data(datash,rr,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)

            plot_data(datasv,rr,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)

            title(['Rotated signals -- P (black), SH (red) and SV (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        elseif superpose==2
            figure(fig_prop1,fig_val1)
            hold all
            if strfind(process.comp,'E')>0
                plot_data(datap,rrx,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datash,rry,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(datasv,rrz,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            
            set(gca,'box','on','ylim',ylim);
            set(gca,'YTick',rr_comb,'YTickLabel',geolab);
            title(['Rotated signals, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end            
        else
            error('wrong superpose option %d! possible values are (0, 1, 2)', superpose);
        end
        fprintf(1,'complete!\n');
    elseif seisplot(i_plot).type==2
        % Plot envelops
        fprintf(1,'plotting envelopes...');
        % Set to original data
        datax=datax_old;
        datay=datay_old;
        dataz=dataz_old;
        for i_rec=1:nrec
            datax(i_rec).data=abs(hilbert(datax_old(i_rec).data));
            datay(i_rec).data=abs(hilbert(datay_old(i_rec).data));
            dataz(i_rec).data=abs(hilbert(dataz_old(i_rec).data));
        end
    
        geomax=zeros(1,nrec);
        if norm_fact==2
            % Geophone maximum
            for i_rec=1:nrec
                geomax(i_rec)=max([max(abs(datax(i_rec).data)) max(abs(datay(i_rec).data)) max(abs(dataz(i_rec).data))]);
            end
        end
        
        if cmap==1
            if strfind(process.comp,'E')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(datax,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['Envelope -- East component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_map_east',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'N')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(datay,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['Envelope -- North component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_map_north',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'Z')>0
                figure(fig_prop1,fig_val1)
                plot_seismap(dataz,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
                title(['Envelope -- Z component map, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_map_vertical',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end
            fprintf(1,'complete!\n');
            continue;
        end
        
        if optional.area>0
            % Compute vectors for stacked regions
            [tvectp,fvectp] = vect4bound({datax.data},[datax.t0]',[datax.dt]',optional.tp,optional.tp+optional.twin(1));
            [tvects,fvects] = vect4bound({datax.data},[datax.t0]',[datax.dt]',optional.ts,optional.ts+optional.twin(2));
        
            % Set optional parameters
            optional.tvectp=tvectp;
            optional.tvects=tvects;
            optional.fvectp=fvectp;
            optional.fvects=fvects;
        end

        if superpose==0
            if strfind(process.comp,'E')>0
                figure(fig_prop1,fig_val1)
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['Envelope -- East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_east',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'N')>0
                figure(fig_prop1,fig_val1)
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['Envelope -- North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_north',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'Z')>0
                figure(fig_prop1,fig_val1)
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['Envelope -- Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_envelope_vertical',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end            
        elseif superpose==1
            figure(fig_prop1,fig_val1)
            hold on
            if strfind(process.comp,'E')>0
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)
            end
            title(['Envelopes -- East (black), North (red) and Z (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_envelopes',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end       
        elseif superpose==2
            figure(fig_prop1,fig_val1)
            hold all
            if strfind(process.comp,'E')>0
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            
            set(gca,'box','on','ylim',ylim);
            set(gca,'YTick',rr_comb,'YTickLabel',geolab);
            title(['Envelopes, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end            
        else
            error('wrong superpose option %d! possible values are (0, 1, 2)', superpose);
        end
        fprintf(1,'complete!\n');
    elseif seisplot(i_plot).type==3
        % Plot rotated envelops
        fprintf(1,'plotting rotated envelopes...');
        
        datap=datax_old;
        datash=datax_old;
        datasv=datax_old;
        % Rotation
        for i_rec=1:nrec
            %datap(i_rec).data =[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*pbasis(i_rec,:)';
            %datash(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*shbasis(i_rec,:)';
            %datasv(i_rec).data=[datax_old(i_rec).data datay_old(i_rec).data dataz_old(i_rec).data]*svbasis(i_rec,:)';
            datap(i_rec).data =(datax_old(i_rec).data*pbasis(i_rec,1))+(datay_old(i_rec).data*pbasis(i_rec,2))+(dataz_old(i_rec).data*pbasis(i_rec,3));
            datash(i_rec).data=(datax_old(i_rec).data*shbasis(i_rec,1))+(datay_old(i_rec).data*shbasis(i_rec,2))+(dataz_old(i_rec).data*shbasis(i_rec,3));
            datasv(i_rec).data=(datax_old(i_rec).data*svbasis(i_rec,1))+(datay_old(i_rec).data*svbasis(i_rec,2))+(dataz_old(i_rec).data*svbasis(i_rec,3));
        end
        % Envelopes
        for i_rec=1:nrec
            datap(i_rec).data=abs(hilbert(datap(i_rec).data));
            datash(i_rec).data=abs(hilbert(datash(i_rec).data));
            datasv(i_rec).data=abs(hilbert(datasv(i_rec).data));
        end
    	
        geomax=zeros(1,nrec);
        if norm_fact==2
            % Geophone maximum
            for i_rec=1:nrec
                geomax(i_rec)=max([max(abs(datap(i_rec).data)) max(abs(datash(i_rec).data)) max(abs(datasv(i_rec).data))]);
            end
        end
        
        if cmap==1
            figure(fig_prop1,fig_val1)
            plot_seismap(datap,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated envelopes map -- P, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_map_p',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            
            if seisplot(i_plot).optional.align==4
                optional.align=2;
            end
            if seisplot(i_plot).optional.area==4
                % Compute vectors for stacked regions
                [tvects,fvects] = vect4bound({datash.data},[datash.t0]',[datash.dt]',optional.ts,optional.ts+optional.twin(2));                
                % Set optional parameters
                optional.tvects=tvects; optional.fvects=fvects;
                optional.area=2; 
                optional.onset=1;
            end
            figure(fig_prop1,fig_val1)
            plot_seismap(datash,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated envelopes map -- SH, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_map_sh',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            
            if seisplot(i_plot).optional.align==4
                optional.align=2;
            end
            if seisplot(i_plot).optional.area==4
                % Compute vectors for stacked regions
                [tvects,fvects] = vect4bound({datasv.data},[datasv.t0]',[datasv.dt]',optional.ts,optional.ts+optional.twin(2));                
                % Set optional parameters
                optional.tvects=tvects; optional.fvects=fvects;
                optional.area=2;
                optional.onset=1;
            end
            figure(fig_prop1,fig_val1)
            plot_seismap(datasv,rr,norm_fact,plot_ampl,clip_ampl,geomax,optional)
            title(['Rotated envelopes map -- SV, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_map_sv',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            fprintf(1,'complete!\n');
            continue;
        end

        if superpose==0
            if seisplot(i_plot).optional.align==4
                optional.align=1;
            end
            if seisplot(i_plot).optional.area==4
                % Compute vectors for stacked regions
                [tvectp,fvectp] = vect4bound({datap.data},[datap.t0]',[datap.dt]',optional.tp,optional.tp+optional.twin(1));                
                % Set optional parameters
                optional.tvectp=tvectp; optional.fvectp=fvectp;
                optional.area=1;
                optional.onset=2;
            end
            figure(fig_prop1,fig_val1)
            plot_data(datap,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated envelopes -- P, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_p',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            
            if seisplot(i_plot).optional.align==4
                optional.align=2;
            end
            if seisplot(i_plot).optional.area==4
                % Compute vectors for stacked regions
                [tvects,fvects] = vect4bound({datash.data},[datash.t0]',[datash.dt]',optional.ts,optional.ts+optional.twin(2));                
                % Set optional parameters
                optional.tvects=tvects; optional.fvects=fvects;
                optional.area=2; 
                optional.onset=1;
            end
            figure(fig_prop1,fig_val1)
            plot_data(datash,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated envelopes -- SH, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_sh',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
            
            if seisplot(i_plot).optional.align==4
                optional.align=2;
            end
            if seisplot(i_plot).optional.area==4
                % Compute vectors for stacked regions
                [tvects,fvects] = vect4bound({datasv.data},[datasv.t0]',[datasv.dt]',optional.ts,optional.ts+optional.twin(2));                
                % Set optional parameters
                optional.tvects=tvects; optional.fvects=fvects;
                optional.area=2;
                optional.onset=1;
            end
            figure(fig_prop1,fig_val1)
            plot_data(datasv,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            title(['Rotated envelopes -- SV, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes_sv',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        elseif superpose==1
            figure(fig_prop1,fig_val1)
            hold on
            plot_data(datap,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)

            plot_data(datash,rr,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)

            plot_data(datasv,rr,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)

            title(['Rotated envelopes -- P (black), SH (red) and SV (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_envelopes',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end        
        else
            error('wrong superpose option %d! possible values are (0, 1, 2)', superpose);
        end
        fprintf(1,'complete!\n');
    elseif seisplot(i_plot).type==4
        % Plot SNR
        fprintf(1,'plotting SNR...');
        %computation of sta, lta and sn for the filtered traces/geophones
        %lwsta=round(seisplot(i_plot).par.sta/datax(1).dt);
        %lwlta=round(seisplot(i_plot).par.lta/datax(1).dt);
        %lwin=lwsta+lwlta;
        
        % Set to original data
        datax=datax_old;
        datay=datay_old;
        dataz=dataz_old;
        
        for i_rec = 1:nrec
            if datax(i_rec).nsamp>1 || datax(i_rec).dt~=0.0
                %computation of sta, lta and sn for the filtered traces/geophones
                lwsta=round(seisplot(i_plot).par.sta/datax(i_rec).dt);
                lwlta=round(seisplot(i_plot).par.lta/datax(i_rec).dt);
                lwin=lwsta+lwlta;
                % Compute SNR
                [datax(i_rec).data, datax(i_rec).lta, datax(i_rec).sta]=compute_SNR(datax(i_rec),lwin,lwlta,lwsta);
            end
            if datay(i_rec).nsamp>1 || datay(i_rec).dt~=0.0
                %computation of sta, lta and sn for the filtered traces/geophones
                lwsta=round(seisplot(i_plot).par.sta/datay(i_rec).dt);
                lwlta=round(seisplot(i_plot).par.lta/datay(i_rec).dt);
                lwin=lwsta+lwlta;
                % Compute SNR
                [datay(i_rec).data, datay(i_rec).lta, datay(i_rec).sta]=compute_SNR(datay(i_rec),lwin,lwlta,lwsta);
            end
            if dataz(i_rec).nsamp>1 || dataz(i_rec).dt~=0.0
                %computation of sta, lta and sn for the filtered traces/geophones
                lwsta=round(seisplot(i_plot).par.sta/dataz(i_rec).dt);
                lwlta=round(seisplot(i_plot).par.lta/dataz(i_rec).dt);
                lwin=lwsta+lwlta;
                % Compute SNR
                [dataz(i_rec).data, dataz(i_rec).lta, dataz(i_rec).sta]=compute_SNR(dataz(i_rec),lwin,lwlta,lwsta);
            end
        end
        
        geomax=zeros(1,nrec);
        if norm_fact==1
            % Geophone maximum
            for i_rec=1:nrec
                geomax(i_rec)=max([max(abs(datax(i_rec).data)) max(abs(datay(i_rec).data)) max(abs(dataz(i_rec).data))]);
            end
        end

        if superpose==0
            if strfind(process.comp,'E')>0
                figure(fig_prop1,fig_val1)
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['SNR -- East component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_snr_east',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'N')>0
                figure(fig_prop1,fig_val1)
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['SNR -- North component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_snr_north',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end

            if strfind(process.comp,'Z')>0
                figure(fig_prop1,fig_val1)
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
                title(['SNR -- Z component, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
                if preinfo.save>0
                    outf=strcat(preinfo.fighead,'_snr_vertical',scount);
                    savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                end
            end
        elseif superpose==1
            figure(fig_prop1,fig_val1)
            hold on
            if strfind(process.comp,'E')>0
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'r',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'g',geomax,optional)
            end
            title(['SNR -- East (black), North (red) and Z (green) components, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_snr',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end        
        elseif superpose==2
            figure(fig_prop1,fig_val1)
            hold all
            if strfind(process.comp,'E')>0
                plot_data(datax,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'N')>0
                plot_data(datay,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            if strfind(process.comp,'Z')>0
                plot_data(dataz,rr,norm_fact,plot_ampl,clip_ampl,'k',geomax,optional)
            end
            
            set(gca,'box','on','ylim',ylim);
            set(gca,'YTick',rr_comb,'YTickLabel',geolab);
            title(['SNR, filter: ' num2str(process.ffreq(1)) ' to ' num2str(process.ffreq(2)) ' Hz']);
            
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_rotated_traces',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end            
        else
            error('wrong superpose option %d! possible values are (0, 1, 2)', superpose);
        end
        fprintf(1,'complete!\n');
    elseif seisplot(i_plot).type==5
        fprintf(1,'plotting polarization...');
        
        % Set to original data
        datax=datax_old;
        datay=datay_old;
        dataz=dataz_old;
        
        geonum=[datax.geonum];
        alpha=cell(1,nrec);
        gamma=cell(1,nrec);
        Q1=cell(1,nrec);
        Q2=cell(1,nrec);
        
        move_pol=seisplot(i_plot).par.move;
        twin_pol=seisplot(i_plot).par.twin;
        if move_pol<=0
            fprintf(1,'WARNING: value of ''move_pol'' changed to 1 from wrong value %d\n',move_pol);
            move_pol=1;
        end
        % Segment below works for geophones having same origin time, sampling
        % rate etc., otherwise use the segement within the loop below.
        ntwin_pol=round(twin_pol/datax(1).dt)+1;
        if (mod(ntwin_pol,2)==0)
            ntwin_pol=ntwin_pol+1;
            twin_pol=(ntwin_pol-1)*datax(1).dt;
        end
        tstart=datax(1).t0+0.5*twin_pol; % Center of the time window
        tend=datax(1).t0+(datax(1).nsamp-1)*datax(1).dt-0.5*twin_pol;  % Center of the time window
        tpol=(tstart:datax(1).dt*move_pol:tend); % Time vector to plot polarization
        np=length(tpol);
        %ntend=datax(1).nsamp-(ntwin_pol-1);
        % Segment...
        
        if seisplot(i_plot).par.map>0
            figure(fig_prop1,fig_val1)
            hold on    
            for i_rec=1:nrec
                %ntwin_pol=round(twin_pol/datax(i_rec).dt)+1;
                %if (mod(ntwin_pol,2)==0)
                %    ntwin_pol=ntwin_pol+1;
                %    twin_pol=(ntwin_pol-1)*datax(i_rec).dt;
                %end        
                %tstart=datax(i_rec).t0+0.5*twin_pol; % Center of the time window
                %tend=datax(i_rec).t0+(datax(i_rec).nsamp-1)*datax(i_rec).dt-0.5*twin_pol;  % Center of the time window
                %tpol=(tstart:datax(i_rec).dt*move_pol:tend); % Time vector to plot polarization
                %np=length(tpol);        
                %ntend=datax(i_rec).nsamp-(ntwin_pol-1);

                loc(1:np)=rr(i_rec);
                ipol=0;
                nwin1=1; % Starts from first smaple
                for i_pol=1:np
                    ipol=ipol+1;                    
                    nwin2=nwin1+(ntwin_pol-1);
                    [alpha{i_rec}(ipol),gamma{i_rec}(ipol),Q1{i_rec}(ipol),Q2{i_rec}(ipol)]=polarize([datax(i_rec).data(nwin1:nwin2) datay(i_rec).data(nwin1:nwin2) dataz(i_rec).data(nwin1:nwin2)]);
                    nwin1=nwin1+move_pol;
                end

                scatter(tpol,loc,15,Q1{i_rec},'o','filled');       

            end
            xlim=[min(tpol) max(tpol)];
            ylim=[min(rr) max(rr)+10];
            set(gca,'box','on','xlim',xlim','ylim',ylim);
            set(gca,'YTick',rr,'YTickLabel',geonum(end:-1:1));
            xlabel('Time (s)');
            ylabel('Receivers');    
            colorbar;
            title('Polarization strength');
            if preinfo.save>0
                outf=strcat(preinfo.fighead,'_polarization_map',scount);
                savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
            end
        end

        if seisplot(i_plot).par.azimuth>0
            pr1=seisplot(i_plot).par.r1;
            prstep=seisplot(i_plot).par.rstep;
            pr2=seisplot(i_plot).par.r2;
            rp=pr1:prstep:pr2;
            nprec=length(rp);
            [tf,irp]=ismember(rp,geonum);
            for i_rec=1:nprec
                if irp(i_rec) ~= 0
                    figure(fig_prop1,fig_val1)
                    scatter(tpol,abs(alpha{irp(i_rec)}),15,Q1{irp(i_rec)},'o','filled');
                    set(gca,'box','on','xlim',xlim,'ylim',[0 360]);
                    xlabel('Time (s)');
                    ylabel('Azimuth (^o)');
                    colormap(jet);
                    tit=sprintf('Receiver: %d - azimuth and polarization strength',rp(i_rec));
                    title(tit);
                    colorbar;
                    if preinfo.save>0
                        outf=strcat(preinfo.fighead,'_polarization_',num2str(rp(i_rec)),scount);
                        savefigure(outf,'Dpi',preinfo.figres,'po',preinfo.figform);
                    end
                else
                    fprintf('WARNING: receiver %d was not processed!\n',rp(i_rec));
                end                
            end
        end
        fprintf(1,'complete!\n');    
    end    
end
%======================================================

% ****************** Other functions ***************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plot_data(data,chnum,geonum,t0,dt,nsamp,plotmask,norm_fact,plot_ampl,clip_ampl.norm,component)
% plots a seismogram montage and returns handle mont
% input : DATA structure with fields
%         data       : vector with the seismograms 
%         chnum      : scalar, channel numbers
%         geonum     : scalar, geophone numbers
%         t0         : scalar, time of first sample of individual seismogram
%         dt         : scalar, time increment for each individual seismogram
%         nsamp      : scalar, number of samples for each individual seismogram
%         plotmask   : scalar, zeroes, ones or twoes (masked, black colored or
%                                                         red colored)
%         norm_fact,plot_ampl,clip_ampl.norm  : normalization on trace (0), geophone (1), or total maximum (2)
%         component  : character 
%         horv       : 1 horizontal seismograms, 2 vertical seismograms
%
%   Michael Roth 21.03.2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HNG,SEP 18,2008
function plot_data(DATPAR,rr,norm_fact,plot_ampl,clip_ampl,plot_color,geomax,varargin)

data={DATPAR.data};
t0=[DATPAR.t0]';
dt=[DATPAR.dt]';
nsamp=[DATPAR.nsamp]';
nrec=length(DATPAR);
%comp=DATPAR(1).comp;
geonum=[DATPAR.geonum]';
optional.onset=0;
optional.align=0;
optional.area=0;
optional.rect=0;
if nargin>7
    optional=set_optional(varargin{1},optional);    
end

normf=ones(1,nrec);

% Normalize first to trace mean if necessary
if norm_fact==4
    for i_rec=1:nrec        
        data{i_rec}=data{i_rec}/mean(abs(data{i_rec}));
    end
    norm_fact=3;
end

%find the trace maxima
tracemax=zeros(1,nrec);
for i_rec=1:nrec
    tracemax(i_rec)=max(abs(data{i_rec}));    
end

%find the absolute maximum of all selected traces
allmax(1:nrec)=max(tracemax);

if (norm_fact == 1)
    normf=tracemax;
elseif (norm_fact == 2)
    normf=geomax;
elseif (norm_fact == 3)
    normf=allmax;
end

% Normalize
splot=cell(1,nrec);
for i_rec=1:nrec
    if normf(i_rec)==0
        splot{i_rec}=data{i_rec}*plot_ampl;
    else
        splot{i_rec}=data{i_rec}/normf(i_rec)*plot_ampl;
    end
end

%clipping
for i_rec=1:nrec;
 isamp=find(abs(splot{i_rec}) > clip_ampl);
 splot{i_rec}(isamp)=sign(splot{i_rec}(isamp))*clip_ampl;
end

tcor=zeros(nrec,1);
% Align
if optional.align>0    
    if optional.align == 1 % Alignment with respect to P
        tmin=min(optional.tp);
        tcor=optional.tp-tmin;
    elseif optional.align == 2 % Alignment with respect to S
        tmin=min(optional.ts);
        tcor=optional.ts-tmin;
    end
    if isfield(optional,'tp')
        optional.tp=optional.tp-tcor;
    end
    if isfield(optional,'ts')
        optional.ts=optional.ts-tcor;
    end
end

if optional.onset>0    
    onset_h = 0.5*clip_ampl;
    onset_w = 2.5;
    onset_wp = 2.0;
    onset_ws = 2.0;
    blue = [0.0 0.4 1.0];
end

hold on

% Plot background rectangle
if optional.rect==1
    if optional.align==1        
        rcol = [1.0 0.7 0.7]; % light red        
        rx=optional.tp(1);
        ry=rr(1)-1.0*clip_ampl;
        rw=optional.twin(1);
        rh=rr(nrec)-rr(1)+2.0*clip_ampl;
    elseif optional.align==2      
        rcol = [0.65 0.85 1.0]; % light blue
        rx=optional.ts(1);
        ry=rr(1)-1.0*clip_ampl;
        rw=optional.twin(2);
        rh=rr(nrec)-rr(1)+2.0*clip_ampl;
    end
    rectangle('Position',[rx,ry,rw,rh],'FaceColor',rcol,'EdgeColor','none','LineWidth',1);
end
maxt=0;

for i_rec=1:nrec                     
    t=(0:nsamp(i_rec)-1)'*dt(i_rec)+t0(i_rec)-tcor(i_rec);

    maxt=max([t; maxt]);
    pos=rr(i_rec);
                       
    plot(t,splot{i_rec}+pos,'Color',plot_color);    
    
    if optional.area==1 || optional.area==3
        if max(size(optional.tvectp{i_rec})) > 0 && max(size(optional.fvectp{i_rec})) > 0
            farea=area(optional.tvectp{i_rec}-tcor(i_rec),optional.fvectp{i_rec}*plot_ampl/normf(i_rec)+pos,pos);
            set(farea,'FaceColor','r','EdgeColor','r');
        end
    end
    if optional.area==2 || optional.area==3
        if max(size(optional.tvects{i_rec})) > 0 && max(size(optional.fvects{i_rec})) > 0
            farea=area(optional.tvects{i_rec}-tcor(i_rec),optional.fvects{i_rec}*plot_ampl/normf(i_rec)+pos,pos);
            set(farea,'FaceColor',blue,'EdgeColor',blue);
        end
    end

    % Plot onsets
    if optional.onset == 1 || optional.onset == 3
        plot([optional.tp(i_rec),optional.tp(i_rec)],[pos-onset_h,pos+onset_h],'k','LineWidth',onset_wp); % ,'r','LineWidth',onset_w);
    end
	if optional.onset == 2 || optional.onset == 3    
        plot([optional.ts(i_rec),optional.ts(i_rec)],[pos-onset_h,pos+onset_h],'color',blue,'LineWidth',onset_ws); % ,blue,'LineWidth',onset_w);
	end
end

irec_on=nsamp>1;

ylim=[min(rr)-clip_ampl max(rr)+clip_ampl];
set(gca,'YTick',rr,'YTickLabel',geonum,'FontSize',14);
set(gca,'box','on','XLIM',[min(t0(irec_on)-tcor(irec_on)) maxt+0.0000001],'ylim',ylim);
xlabel('Time (s)');
ylabel('Receivers');
%set(gca,'FontSize',12);
%xlhand=get(gca,'xlabel');
%set(xlhand,'FontSize',20);
%ylhand=get(gca,'ylabel');
%set(ylhand,'FontSize',20);


%hold off
%======================================================

function plot_seismap(DATPAR,rr,norm_fact,plot_ampl,clip_ampl,geomax,varargin)
% this function plots the seismgram map in colorscale

data={DATPAR.data};
t0=[DATPAR.t0]';
dt=[DATPAR.dt]';
nsamp=[DATPAR.nsamp]';
nrec=length(DATPAR);
%comp=DATPAR(1).comp;
geonum=[DATPAR.geonum]';

optional.onset=0;
optional.align=0;
optional.area=0;
optional.rect=0;
if nargin>6
    varargin{1}.area=0; % we don't want to plot area in seismogram map
    optional=set_optional(varargin{1},optional);    
end

normf=ones(1,nrec);
% Normalize first to trace mean if necessary
if norm_fact==4
    for i_rec=1:nrec        
        data{i_rec}=data{i_rec}/mean(abs(data{i_rec}));
    end
    norm_fact=3;
end

%find the trace maxima
tracemax=zeros(1,nrec);
for i_rec=1:nrec
    tracemax(i_rec)=max(abs(data{i_rec}));    
end

%find the absolute maximum of all selected traces
allmax(1:nrec)=max(tracemax);

if (norm_fact == 1)
    normf=tracemax;
elseif (norm_fact == 2)
    normf=geomax;
elseif (norm_fact == 3)
    normf=allmax;
end

% Normalize
splot=cell(1,nrec);
for i_rec=1:nrec
    if normf(i_rec)==0
        splot{i_rec}=data{i_rec}*plot_ampl;
    else
        splot{i_rec}=data{i_rec}/normf(i_rec)*plot_ampl;
    end
end

% t vector
tmin=min(t0);
tmax=max(t0+(nsamp-1).*dt);

dt_new=max([0 min(dt)]);
nsamp_new=round((tmax-tmin)/dt_new)+1;

tvect_new=(0:nsamp_new-1)'*dt_new+tmin;

% generate colormap values
cmap=zeros(nrec,nsamp_new);
for i_rec=1:nrec
    tvect=(0:nsamp(i_rec)-1)'*dt(i_rec)+t0(i_rec);
    if nsamp(i_rec) ~= nsamp_new || ~isequal(tvect,tvect_new)
        tvect=(0:nsamp(i_rec)-1)'*dt(i_rec)+t0(i_rec);
        cmap(i_rec,:)=interp1(tvect,splot{i_rec},tvect_new,'linear');
    else
        cmap(i_rec,:)=splot{i_rec};
    end
end

hold on
%imagesc(tvect_new,rr,cmap); % imagesc always plots top o botom therefore flipud is necessary
pcolor(tvect_new,rr,cmap); % imagesc always plots top o botom therefore flipud is necessary
% ---------------------------------------------------------------
% colormap code was used from http://meyavuz.blogspot.com/2011/02/meep-colormap-for-matlab.html
mincolor    = [1 0 0]; % red
mediancolor = [1 1 1]; % white   
maxcolor    = [0 0 1]; % blue      

ColorMapSize = 64; % lower value will supress smaller values
int1 = zeros(ColorMapSize,3); 
int2 = zeros(ColorMapSize,3);
for k=1:3
    int1(:,k) = linspace(mincolor(k), mediancolor(k), ColorMapSize);
    int2(:,k) = linspace(mediancolor(k), maxcolor(k), ColorMapSize);
end
meep = [int1(1:end-1,:); int2];

colormap(meep); 
shading flat;
colorbar
% ---------------------------------------------------------------

% plot onsets
if optional.onset>0   
    onset_h = 0.5*clip_ampl;
    onset_w = 2.5;
    blue = [0.0 0.4 1.0];

    for i_rec=1:nrec      

        pos=rr(i_rec);

        % Plot onsets
        if optional.onset == 1 || optional.onset == 3
            plot([optional.tp(i_rec),optional.tp(i_rec)],[pos-onset_h,pos+onset_h],'r','LineWidth',onset_w);
        end
        if optional.onset == 2 || optional.onset == 3    
            plot([optional.ts(i_rec),optional.ts(i_rec)],[pos-onset_h,pos+onset_h],'color',blue,'LineWidth',onset_w);
        end
    end
end

ylim=[min(rr) max(rr)];
set(gca,'box','on','XLIM',[tmin tmax],'ylim',ylim);
set(gca,'YTick',rr,'YTickLabel',geonum);

ylabel('Receivers');
xlabel('Time (s)','FontSize',20);
%hold off
%======================================================

% Funtion to set optiona cell variables
function optional=set_optional(cellvar, optional)

% Set main variables 'onset', 'align', 'area'
if isfield(cellvar,'onset')
    optional.onset=cellvar.onset;
end
if isfield(cellvar,'align')
    optional.align=cellvar.align;
end
if isfield(cellvar,'area')
    optional.area=cellvar.area;
end

% Check and set 'tp'
if (optional.onset==1 || optional.onset==3) || ...
   (optional.align==1 || optional.align==3) || ...
   (optional.area==1  || optional.area==3)
    if isfield(cellvar,'tp')
        optional.tp=cellvar.tp;
    else
        fprintf(1,'WARNING: field ''tp'' not found in optional cell variable!\n optional choices onset, align and area may be deactivated!');
        optional.onset=max([0 optional.onset]);
        optional.align=max([0 optional.align]);
        optional.area=max([0 optional.area]);
    end
end

% Check and set 'ts'
if (optional.onset==2 || optional.onset==3) || ...
   (optional.align==2 || optional.align==3) || ...
   (optional.area==2  || optional.area==3)
    if isfield(cellvar,'ts')
        optional.ts=cellvar.ts;
    else
        fprintf(1,'WARNING: field ''ts'' not found in optional cell variable!\n optional choices - onset, align and area may be deactivated!');
        optional.onset=max([0 optional.onset]);
        optional.align=max([0 optional.align]);
        optional.area=max([0 optional.area]);
    end
end

% Check and set 'tvectp' and 'fvectp'
if optional.area==1 || optional.area==3
    if isfield(cellvar,'tvectp')
        optional.tvectp=cellvar.tvectp;
    else
        error('field ''tvectp'' not found in optional cell variable to plot stacked area for P!');
    end
    if isfield(cellvar,'fvectp')
        optional.fvectp=cellvar.fvectp;
    else
        error('field ''fvectp'' not found in optional cell variable to plot stacked area for P!');
    end
end

% Check and set 'tvects' and 'fvects'
if optional.area==2 || optional.area==3
    if isfield(cellvar,'tvects')
        optional.tvects=cellvar.tvects;
    else
        error('field ''tvects'' not found in optional cell variable to plot stacked area for S!');
    end
    if isfield(cellvar,'fvects')
        optional.fvects=cellvar.fvects;
    else
        error('field ''fvects'' not found in optional cell variable to plot stacked area for S!');
    end
end

% Check and set 'rect'
if optional.align>0 && isfield(cellvar,'rect')
    optional.rect=cellvar.rect;        
end

% Check and set 'twin'
if (optional.rect==1)
    if isfield(cellvar,'twin')
        optional.twin=cellvar.twin;
    else
        fprintf(1,'WARNING: field ''twin'' not found in optional cell variable!\n optional choice ''rect'' will be deactivated!');
        optional.rect=0;
    end
end
% end function optional--------------------------------
%======================================================


function [SNR, lta, sta] = compute_SNR(data,lwin,lwlta,lwsta)
%computation of cumulative sum over all channels of a geophone
cuss=cumsum(sum(abs(data.data),2));
sta=(cuss(lwin+1:data.nsamp)-cuss(lwin+1-lwsta:data.nsamp-lwsta))/lwsta;
lta=(cuss(lwlta+1:data.nsamp-lwsta)-cuss(1:data.nsamp-lwin))/lwlta;
zz = find(lta==0);
lta(zz)=1;
sta(zz)=1;
sn=sta./lta;
% padding
sta=[zeros(lwin,1);sta];
lta=[zeros(lwin,1);lta];
SNR=[ones(lwin,1);sn];
%======================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [alpha,gamma,Q1,Q2]=polaris_func(t_start,t_win,nsamp,dt,data,t0,comp,freq,ord)
% computes polarization within the time specified time window 
% 
% input: 
%        t_start: start time [s] of analyse window for the geophone, scalar
%        t_win  : length of the time window, scalar 
%        nsamp  : number of the seismogram samples, scalar integer
%        dt     : sampling interval, scalar
%        data   : seismograms, structure with 3 vectors data{1:3}(1:nsamp) 
%                 sign convention for particle motion:
%                                  E (- in west direction, + in east direction)
%                                  N (- in south direction, + in north direction)
%                                  Z (- down, + up)
%        t0     : time of the first sample, scalar 
%        comp   : channel-component assignement
%                 sructure comp{1:3} containing the identifier E N or Z
% 
%                  
% output:
%        alpha  : azimuth, scalar
%                 +- 180 deg ambiguity !!!
%                 (north 0 deg -> east 90 -> south 180 -> west 270 -> north)
%        gamma  : incidence angle, scalar
%                 (-z direction 0 deg -> +z direction 180 deg)
%                 90+- 2*(gamma-90) ambiguity !!!
%        Q1     : strength of polarization (-1 - 1)
%                 (-1: no prefered polarization direction,  
%                   0: the largest eigenvalue equals the sum of intermediate and smallest eigenvalue 
%                   1: perfect linear polarization)
%        Q2     : degree of planar polarization (0 - 1)
%                 (0: intermediate and smallest components are comparable,
%                  1: intermediate component is much larger than smallest component)
%
% Michael Roth 04.03.03
%
% bug fix in length of the analysis window
% Volker Oye 07.03.07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function [alpha,gamma,Q1,Q2]=polarize(data_win)

if norm(data_win,inf)<=0
    fprintf('WARNING: matrix norm=0 for polarization!\n');
    alpha=0; gamma=0;
    Q1=0; Q2=0;
    return;
end

%normalize frequency for filtering 
% freq=freq*dt*2;
% freq(1)=max([0.001 freq(1)]);
% freq(2)=min([freq(2) 0.999]);

%computation of filter coefficients 
% [B,A]=butter(ord,freq);
% filt_len=10*ord;

%find order of traces
% n1=find(strcmp('E',comp));
% n2=find(strcmp('N',comp));
% n3=find(strcmp('Z',comp));
 
%sorting of the ENZ-components
%dat=[data{n1} data{n2} data{n3}];

%compute sample indices for start and end of time window 
% it_start=round((t_start-t0)/dt); 
% it_end=it_start+round(t_win/dt);

% %start index for filter window
% i_filwin1=max([it_start-filt_len 1]);
% %end index for filter window
% i_filwin2=min([it_end+filt_len nsamp]);
% 
% if i_filwin1 == 1
%   it_start=filt_len;
%   it_end=it_start+round(t_win/dt);
%   i_filwin2 = min([it_end+filt_len nsamp]);
% end
% if i_filwin2 == nsamp
%   i_filwin1 = i_filwin2-2*filt_len-round(t_win/dt);
%   it_start=i_filwin1+filt_len;
%   it_end=it_start+filt_len+round(t_win/dt);
% end
% 
% filt_win=i_filwin1:i_filwin2;
% 
% %filtering of data
% dat(filt_win,:)=filter(B,A,dat(filt_win,:));

%compute covarianz matrix
%cov_mat=cov(dat(it_start:it_end,:));
cov_mat=cov(data_win);
  
%solve eigenwert problem (V=eigenvectors, D=eigenvalues)
[V,D]=eig(cov_mat);
  
%sort eigenvalues (dum(1)<dum(2)<dum(3))
[dum,ii]=sort(diag(D));
  
%compute quality factors
Q1=1-(dum(1)+dum(2))/dum(3);
Q2=1-dum(1)/dum(2);
 
%get eigenvector for the largest eigenvalue 
eigenvec=V(:,ii(3));
  
%compute azimuth
alpha=atan2(eigenvec(1),eigenvec(2))*180/pi;
if alpha<0.;
  alpha=alpha+360.; 
end;
  
%compute incidence angle
gamma=acos(-eigenvec(3))*180/pi;
%======================================================

% function vect4bound
function [tvect,fvect] = vect4bound(data,data_t0,data_dt,data_t1,data_t2)
% PURPOSE
%   This function computes the definite integration of a equi-spaced discrete function
% INPUT
%   data: vector containing function values
%   data_t0: Initial time of recording
%   data_dt: spacing of function values
%   t1,t2: Bound for integration, in principle bound can be defined anywhere in the t-axis
% CHANGE HISTORY
%   HNG,OCT 02,2009;HNG,SEP 18,2008;HNG,Aug 15,2008
nrec=length(data);
tvect=cell(1,nrec);
fvect=cell(1,nrec);
for i_rec=1:nrec
    t1=data_t1(i_rec)-data_t0(i_rec);
    t2=data_t2(i_rec)-data_t0(i_rec);
    dt=data_dt(i_rec);
    
    nlen=max(size(data{i_rec}));
    t(1:nlen,1)=(0:nlen-1)*dt;   
    
	tvect{i_rec}=[];
	fvect{i_rec}=[];
   
    if nlen==0 || t2 < 0 || t1 > t(nlen) || t1 == t2
        return;
    end

    if t1 > t2
        error('Bad bound for integration: t1 > t2!');
    end

    if t1 < 0
        t1=0;
    end

    if t2 > t(end)
        t2=t(end); %tend;
    end

    % Indices of the first and the last values
    ind0=floor(t1/dt)+1;
    indn=ceil(t2/dt)+1;

    ind1=ind0+1; % Index just after first value f0
    ind2=indn-1; % Index just before last value fn

    % Find the first and the last values of the region by interpolation
    f0=interp_line(t(ind0:ind1),data{i_rec}(ind0:ind1),t1);
    fn=interp_line(t(ind2:indn),data{i_rec}(ind2:indn),t2);        

    n=ind2-ind1+1; % Total number of intermediate values

    tvect{i_rec}=zeros(n+2,1); fvect{i_rec}=zeros(n+2,1);

    tvect{i_rec}(1)=t1; fvect{i_rec}(1)=f0;
    tvect{i_rec}(n+2)=t2; fvect{i_rec}(n+2)=fn;

    % Compute area using trapezoidal rule
    if n == 1    
        tvect{i_rec}(2)=(ind1-1)*dt;
        fvect{i_rec}(2)=data{i_rec}(ind1);
    elseif n > 1    
        tvect{i_rec}(2:n+1)=(ind1-1:ind2-1)*dt;
        fvect{i_rec}(2:n+1)=data{i_rec}(ind1:ind2);
    end
    tvect{i_rec}=tvect{i_rec}+data_t0(i_rec);  
end
%---function vect4bound---
%======================================================

%---function interp_line---
% This function interpolates a point (xp,yp) on a line given by two points
% (x(1),y(1)) and (x(2),y(2)) to find yp.
function [yp] = interp_line(x,y,xp)
dx=x(2)-x(1); dy=y(2)-y(1);

if dy == 0
    yp=y(1);
    return;
end

if dx == 0
    error('ERROR: wrong values for the interpolation!');
end

yp=y(1)+(xp-x(1))*dy/dx;
return
%---function interp_line---
%======================================================
