function [DATA_NEW] = select_receivers(process,DATA)
% Revision
%   HNG, Nov 11,2009
% Geophone numbers in DATA
[rec_id,rec_chan]=unique([DATA.geonum]); % rec_chan gives the maximum channel index of each receiver
nrec=length(rec_id);

% Determine which geophones to process
if process.select==0 % All receivers
    rec_process=rec_id;
elseif process.select==1 % Indicial selection
    rec_process=process.r1:process.rstep:process.r2;
elseif process.select==2 % Circular selection
    % Compute distances
    recx=reshape([DATA(rec_chan).xyz],3,nrec)'; % reshape fill the elements column-wise
    dist=sqrt((recx(:,1)-process.refx(1)).^2+(recx(:,2)-process.refx(2)).^2+(recx(:,3)-process.refx(3)).^2);
    rec_ind=unique(find(dist<=process.rad));
    rec_process=rec_id(rec_ind);
elseif process.select==3 % List from file
    recn_inpf=fopen(process.rnfile,'r');
    fgetl(recn_inpf); % skip 1 line
    rec_process=fscanf(recn_inpf,'%d',[1 inf]);
    fclose(recn_inpf);
elseif process.select==-4 % less than the reference excluding the reference
    rec_process=find(rec_id<process.refrec);
elseif process.select==4 % larger than the reference excluding the reference
    rec_process=find(rec_id>process.refrec);
else
    error('Wrong process.select option!');
end

% List of receivers to be discarded
% Determine which geophones to process
if process.xrec==0 % All receivers
    xrec_list=[];
elseif process.xrec==1 % Indicial selection
    xrec_list=process.xr1:process.xrstep:process.xr2;
elseif process.xrec==2 % Circular selection
    % Compute distances
    recx=reshape([DATA(rec_chan).xyz],3,nrec)'; % reshape fill the elements column-wise
    dist=sqrt((recx(:,1)-process.xrefx(1)).^2+(recx(:,2)-process.xrefx(2)).^2+(recx(:,3)-process.xrefx(3)).^2);
    rec_ind=unique(find(dist<=process.xrad));
    xrec_list=rec_id(rec_ind);
elseif process.xrec==3 % List from file
    recn_inpf=fopen(process.xrnfile,'r');
    fgetl(recn_inpf); % skip 1 line
    xrec_list=fscanf(recn_inpf,'%d',[1 inf]);
    fclose(recn_inpf);
elseif process.xrec==-4 % less than the reference excluding the reference
    xrec_list=find(rec_id<process.xrefrec);
elseif process.xrec==4 % larger than the reference excluding the reference
    xrec_list=find(rec_id>process.xrefrec);
else
    error('Wrong process.select option!');
end
clear dist rec_ind

% Throw other geophones' data
rec_throw=[setdiff(rec_id,rec_process) xrec_list];
for i_throw=1:nrec % maximum recursion is nrec
	[~,ichan_throw]=ismember(rec_throw,[DATA.geonum]); % This gives only the maximum index therefore I need to do recursively    
	if isempty(ichan_throw) || max(ichan_throw)==0
        break;
	end
    chan_throw=nonzeros(ichan_throw)';    
	
    DATA(chan_throw)=[];    
end
DATA_NEW=DATA;
