function [preinfo, rdata, process, seisplot, deplot, model] = read_plot(fname)
%----------------------------------------------------------------
% This function reads the main input file for mig3d routines
% Revision
%   HNG, Dec 08,2009; HNG, Nov 11,2009
if nargin==0
    error('no input file name!')
end
if ~exist(fname,'file')
    error('file ''%s'' doesn''t exist!',fname);
end

MAXNLINE=1000; % maximum number of lines in the input file

% Initialize output structures
preinfo=[];
model=[];
rdata=[];
process=[];
deplot=[];

need_model=0;
count=0;

inpf=fopen(fname,'r');

% Default values
preinfo_stat=0; rdata_stat=0; process_stat=0; seisplot_stat=0; deplot_stat=0; model_stat=0; 
rdata.nrec=[]; rdata.rfile=[];
process.select=0; process.xrec=0; process.comp='ENZP'; process.tclip=[0 0];  process.noise=0; process.ffreq=[0 0]; process.norm=0;
process.erpf=0; process.ncoef=[]; process.tnoise=[]; process.sort_src=[];
seisplot=[];

% Start reading file
while ~feof(inpf)    
    str_line=fgetl(inpf);
    if isempty(str_line)
        continue;
    end
    [token, remain]=strtok(str_line);
    % Ignore comment lines
    if token(1)=='#' || isempty(remain)
        continue;
    end
    
    % Look for continuation line
    for i_line=1:MAXNLINE
        if remain(end) ~= '&'
            remain=regexprep(remain,'&',' '); % Replace '&' with space
            break;
        end        
        remain=strcat(remain,fgetl(inpf));        
    end 
    
    % Read receivers data information    
    if strcmp(token,'preinfo')
        split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
        list=[split_line{:}];
        nv=length(list);
        vnames=list(1:2:nv-1);
        values=list(2:2:nv);
        
        preinfo.type=get_nval('type',vnames,values);
        preinfo.plot=get_nval('plot',vnames,values);
        preinfo.save=get_nval('save',vnames,values);
        if preinfo.save==1
            preinfo.figform=get_sval('figform',vnames,values);
            preinfo.figres=get_sval('figres',vnames,values);
            preinfo.fighead=get_sval('fighead',vnames,values);        
        end        

        preinfo_stat=1;
        continue
    end
    % Read receivers data...
    
    % Read receivers data information    
    if strcmp(token,'rdata')
        split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
        list=[split_line{:}];
        nv=length(list);
        vnames=list(1:2:nv-1);
        values=list(2:2:nv);
        
        rdata.type=get_nval('type',vnames,values);
        switch rdata.type
            case 0
                rdata.fheader=get_sval('fheader',vnames,values);
                rdata.nrec=get_nval('nrec',vnames,values);
                rdata.rfile=get_sval('rfile',vnames,values);                
            case 1
                rdata.mfile=get_sval('mfile',vnames,values);
                tmp=look_sval('rfile',vnames,values);
                if ~isempty(tmp)
                    rdata.rfile=tmp;
                end
            case 2 % SEG2 data
                rdata.seg2file=get_sval('seg2file',vnames,values);                
                rdata.rfile=get_sval('rfile',vnames,values);
                rdata.rsetfile=get_sval('rsetfile',vnames,values);
                rdata.nrec=get_nval('nrec',vnames,values);                
            case 3 % SPECFEM3D data
                rdata.fheader=get_sval('fheader',vnames,values);
                rdata.sfile=get_sval('sfile',vnames,values);
                rdata.ext=get_sval('ext',vnames,values);
            case 32 % SPECFEM3D data difference
                rdata.fheader=get_sval('fheader',vnames,values);
                rdata.fheader1=get_sval('fheader1',vnames,values);
                rdata.sfile=get_sval('sfile',vnames,values);
                rdata.ext=get_sval('ext',vnames,values);
            case 4 %SEGYdata
              rdata.segyfile=get_sval('segyfile',vnames,values);
              rdata.rfile=get_sval('rfile',vnames,values);
            otherwise
                error('wrong ''rdata.source'' option: %d! valid options are (0,1,2,3,4)',rdata.source);
        end  

        rdata_stat=1;
        continue
    end
    % Read receivers data...
    
    % Read process information    
    if strcmp(token,'process')
        split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
        list=[split_line{:}];
        nv=length(list);
        vnames=list(1:2:nv-1);
        values=list(2:2:nv);
        
        tmp=look_nval('select',vnames,values);
        if ~isempty(tmp)
            process.select=tmp;
        end
        switch process.select
            case 0
                % do nothing
            case 1
                process.r1=get_nval('r1',vnames,values);
                process.rstep=get_nval('rstep',vnames,values);
                process.r2=get_nval('r2',vnames,values);
            case 2                       
                process.refx=get_nval('refx',vnames,values);
                process.rad=get_nval('rad',vnames,values);
            case 3        
                process.rnfile=get_sval('rnfile',vnames,values);
            case -4                       
                process.refrec=get_nval('refrec',vnames,values);
            case 4                       
                process.refrec=get_nval('refrec',vnames,values);                
            otherwise
                error('wrong ''process.select'' option: %d! valid options are (0, 1, 2, 3)',process.select);
        end
        tmp=look_nval('xrec',vnames,values);
        if ~isempty(tmp)
            process.xrec=tmp;
        end
        switch process.xrec
            case 0
                % do nothing
            case 1
                process.xr1=get_nval('xr1',vnames,values);
                process.xrstep=get_nval('xrstep',vnames,values);
                process.xr2=get_nval('xr2',vnames,values);
            case 2                       
                process.xrefx=get_nval('xrefx',vnames,values);
                process.xrad=get_nval('xrad',vnames,values);
            case 3        
                process.xrnfile=get_sval('xrnfile',vnames,values);
            case -4                       
                process.xrefrec=get_nval('xrefrec',vnames,values);
            case 4                       
                process.xrefrec=get_nval('xrefrec',vnames,values);
            otherwise
                error('wrong ''process.xrec'' option: %d! valid options are (0, 1, 2, 3)',process.select);
        end
        
        tmp=look_sval('comp',vnames,values);
        if ~isempty(tmp)
            process.comp=tmp;
        end
        tmp=look_nval('tclip',vnames,values);
        if ~isempty(tmp)
            process.tclip=tmp;
        end
        tmp=look_nval('noise',vnames,values);
        if ~isempty(tmp)
            process.noise=tmp;
        end
        tmp=look_nval('ffreq',vnames,values);
        if ~isempty(tmp)
            process.ffreq=tmp;
        end
        if max(process.ffreq)>0
            process.ftype=get_nval('ftype',vnames,values);
            process.forder=get_nval('forder',vnames,values);
        end        
        tmp=look_nval('erpf',vnames,values);
        if ~isempty(tmp)
            process.erpf=tmp;
        end
        if max(process.erpf)>0
            process.tnoise=get_nval('tnoise',vnames,values);
            process.ncoef=get_nval('ncoef',vnames,values);
        end
        tmp=look_nval('norm',vnames,values);
        if ~isempty(tmp)
            process.norm=tmp;
        end
        
        tmp=look_nval('sort_src',vnames,values);
        if ~isempty(tmp)
            process.sort_src=tmp;
        end

        process_stat=1;
        continue
    end
    % Read process...    
   
	if preinfo_stat==1
        switch preinfo.type
            case 0                                
                % Read seisplot field parameters
                if strcmp(token,'seisplot')                    
                    split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
                    list=[split_line{:}];
                    nv=length(list);
                    vnames=list(1:2:nv-1);
                    values=list(2:2:nv);
                    
                    count=count+1;
                    
                    % default values
                    seisplot(count).par.cmap=0;
                    
                    seisplot(count).type=get_nval('type',vnames,values);                 
                    
                    if seisplot(count).type==1 || seisplot(count).type==3 % rotation
                        need_model=2;
                    end
                    if seisplot(count).type==4 % SNR
                        seisplot(count).par.lta=get_nval('lta',vnames,values);
                        seisplot(count).par.sta=get_nval('sta',vnames,values);
                    end                  
                    

                    if seisplot(count).type<=4
                        % Parameters
                        tmp=look_nval('cmap',vnames,values);
                        if ~isempty(tmp) && tmp>0
                            seisplot(count).par.cmap=tmp;
                        end
                        
                        seisplot(count).par.norm=get_nval('norm',vnames,values);
                        seisplot(count).par.amp=get_nval('amp',vnames,values);
                        seisplot(count).par.clipamp=get_nval('clipamp',vnames,values);
                    end
                    % Optional parameters
                    % Default optional parameters
                    seisplot(count).optional.superpose=0;
                    seisplot(count).optional.onset=0;
                    seisplot(count).optional.align=0;
                    seisplot(count).optional.area=0;
                    seisplot(count).optional.rect=0;
                    
                    tmp=look_nval('superpose',vnames,values);
                    if ~isempty(tmp) && tmp>0
                        seisplot(count).optional.superpose=tmp;
                    end
                    tmp=look_nval('onset',vnames,values);
                    if ~isempty(tmp)  && tmp>0
                        seisplot(count).optional.onset=tmp;
                        need_model=1;
                    end
                    tmp=look_nval('align',vnames,values);
                    if ~isempty(tmp)  && tmp>0
                        seisplot(count).optional.align=tmp;
                        need_model=1;
                    end
                    tmp=look_nval('area',vnames,values);
                    if ~isempty(tmp)  && tmp>0
                        seisplot(count).optional.area=tmp;
                        need_model=1;
                    end
                    tmp=look_nval('rect',vnames,values);
                    if ~isempty(tmp)  && tmp>0
                        seisplot(count).optional.rect=tmp;
                        need_model=1;
                    end
                    if seisplot(count).optional.area>0
                        seisplot(count).optional.twin=get_nval('twin',vnames,values);
                    end
                    
                    if seisplot(count).type==5 % Polarization
                        seisplot(count).par.map=0;
                        seisplot(count).par.azimuth=0;
                        seisplot(count).par.twin=get_nval('twin',vnames,values);
                        seisplot(count).par.move=get_nval('move',vnames,values);
                        tmp=look_nval('map',vnames,values);
                        if ~isempty(tmp) && tmp>0
                            seisplot(count).par.map=tmp;
                        end
                        tmp=look_nval('azimuth',vnames,values);
                        if ~isempty(tmp) && tmp>0
                            seisplot(count).par.azimuth=tmp;
                        end
                        if seisplot(count).par.azimuth>0
                            seisplot(count).par.r1=get_nval('r1',vnames,values);
                            seisplot(count).par.rstep=get_nval('rstep',vnames,values);
                            seisplot(count).par.r2=get_nval('r2',vnames,values);
                        end
                        
                        if seisplot(count).par.map<1 && seisplot(count).par.azimuth<1
                            fprintf('WARNING: nothing to plot for type=5 (polarization)!\n');
                            seisplot(count)=[];
                            count=count-1;
                        end
                    end
                    
                    seisplot_stat=1;
                    
                    continue
                end
                % Read seisplot...
            case 1
                % Read differential evolution parameters    
                if strcmp(token,'deplot')
                    % Default values
                    deplot.lastpop=0; deplot.allpop=0; deplot.par_gen=0; deplot.par_of=0;
                    deplot.par_par=0; deplot.gen_of=0; deplot.gen_sd=0; deplot.anim=0;
                    
                    split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
                    list=[split_line{:}];
                    nv=length(list);
                    vnames=list(1:2:nv-1);
                    values=list(2:2:nv);
                    
                    deplot.range=get_nval('range',vnames,values);
                    deplot.mfile=get_sval('mfile',vnames,values);
                    deplot.lunit=get_sval('lunit',vnames,values);
                    deplot.rtick=get_nval('rtick',vnames,values);
                    deplot.ntick=get_nval('ntick',vnames,values);
                    tmp=look_nval('lastpop',vnames,values);
                    if ~isempty(tmp)
                        deplot.lastpop=tmp;
                    end
                    tmp=look_nval('allpop',vnames,values);
                    if ~isempty(tmp)
                        deplot.allpop=tmp;
                    end
                    tmp=look_nval('par_gen',vnames,values);
                    if ~isempty(tmp)
                        deplot.par_gen=tmp;
                    end
                    tmp=look_nval('par_of',vnames,values);
                    if ~isempty(tmp)
                        deplot.par_of=tmp;
                    end
                    tmp=look_nval('par_par',vnames,values);
                    if ~isempty(tmp)
                        deplot.par_par=tmp;
                    end
                    tmp=look_nval('gen_of',vnames,values);
                    if ~isempty(tmp)
                        deplot.gen_of=tmp;
                    end
                    tmp=look_nval('gen_sd',vnames,values);
                    if ~isempty(tmp)
                        deplot.gen_sd=tmp;
                    end
                    tmp=look_nval('anim',vnames,values);
                    if ~isempty(tmp)
                        deplot.anim=tmp;
                    end
                    
                    deplot_stat=1;
                    continue
                end            
                % Read differential...
            otherwise
                error('wrong ''compute.type'' option: %d! valid options are (0,1)',compute.type);
        end % switch
	end %switch compute.type
        
    if need_model>0
        % Read model information    
        if strcmp(token,'model')
            split_line=regexp(strtrim(regexp(remain,',','split')),'=','split');
            list=[split_line{:}];
            nv=length(list);
            vnames=list(1:2:nv-1);
            values=list(2:2:nv);
            
            model.lookup=get_nval('lookup',vnames,values);
            switch model.lookup
                case 0
                    model.vp=get_nval('vp',vnames,values);
                    model.vs=get_nval('vs',vnames,values);
                case {1,2}
                    model.pheader=get_sval('pheader',vnames,values);
                    model.sheader=get_sval('sheader',vnames,values);
                otherwise
                    error('wrong ''model.lookup'' option: %d! valid options are (0, 1, 2)',model.lookup);
            end
            model.src=[0.0; 0.0; 0.0; 0.0];
            if need_model==2
                model.lookup=-1;
            end
            if model.lookup==2 || model.lookup==-1
                % other model is optional (for origin time)
                tmp=look_nval('dh',vnames,values);
                if ~isempty(tmp)
                    model.dh=tmp;
                end
                tmp=look_nval('ox',vnames,values);
                if ~isempty(tmp)
                    model.ox=tmp;
                end
                tmp=look_nval('nx',vnames,values);
                if ~isempty(tmp)
                    model.nx=tmp;
                end
                tmp=look_nval('src',vnames,values);
                if ~isempty(tmp)
                    model.src=tmp;
                end                
            else 
                % other model information is compulsory
                model.dh=get_nval('dh',vnames,values);
                model.ox=get_nval('ox',vnames,values);
                model.nx=get_nval('nx',vnames,values);
                model.src=get_nval('src',vnames,values);
            end
            
            model_stat=1;
            continue
        end
        % Read model...
    end % if need_model>0 
end % while ~feof(inpf)

% Check status
if preinfo_stat ~= 1
    error('Preliminary information not found! use ''preinfo'' line with legitimate information.');
end
if preinfo.type==0 && rdata_stat ~= 1
    error('receivers data information not found! use ''rdata'' line with legitimate information.');
end
if preinfo.type==0 && process_stat ~= 1
    error('process information not found! use ''process'' line with legitimate information.');
end
if preinfo.type==0 && seisplot_stat ~= 1
    error('seismogram plotting information not found! use ''seisplot'' line with legitimate information.');
end
if preinfo.type==1 && deplot_stat ~= 1
    error('DE plotting information not found! use ''deplot'' line with legitimate information.');
end
if need_model>0 && model_stat ~= 1
    error('model information not found! use ''model'' line with legitimate information.');
end
% function read_input
%----------------------------------------------------------------

%----------------------------------------------------------------
% Auxilliary routines
% Get numeric value/s
function [val] = get_nval(v,vnames, values)
vind=strcmp(v,vnames);
if any(vind)
    [val, stat]=str2num(char(values(vind)));
    if ~stat
        error('illegal value for ''%s''!\n',v)
    end	
else
	error('variable ''%s'' not found!',v);
end

% Get string value/s
function [val] = get_sval(v,vnames, values)
vind=strcmp(v,vnames);
if any(vind)
    val=char(strtrim(values(vind)));	
else
	error('variable ''%s'' not found!',v);
end

% Look for numeric value/s. Return empty if not found.
function [val] = look_nval(v,vnames, values)
val=[];
vind=strcmp(v,vnames);
if any(vind)
	[val,stat]=str2num(char(values(vind)));
    if ~stat
        error('illegal value for %s!\n',v)
    end
end

function [val] = look_sval(v,vnames, values)
val=[];
vind=strcmp(v,vnames);
if any(vind)
	val=char(strtrim(values(vind)));
end

