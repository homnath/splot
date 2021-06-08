function [DATA_NEW] = throw_data(process,DATA)
% Geophone numbers in DATA
[rec_id,rec_chan]=unique([DATA.geonum]); % rec_chan gives the maximum channel index of each receiver
nrec=length(rec_id);

% Determine which geophones to process
if process.select==0 % All receivers
    rec_process=1:nrec;
elseif process.select==1 % Indicial selection
    rec_process=process.r1:process.rstep:process.r2;
elseif process.select==2 % Circular selection
    % Compute distances
    recx=reshape([DATA(rec_chan).xyz],3,nrec)'; % reshape fill the elements column-wise
    dist=sqrt((recx(:,1)-process.refx(1)).^2+(recx(:,2)-process.refx(2)).^2+(recx(:,3)-process.refx(3)).^2);
    rec_ind=unique(find(dist<=process.rad));
    rec_process=rec_id(rec_ind);
elseif process.select==3 % List from file
    rec_select_inpf=fopen(rec_select_filename,'r');
    rec_process=fscanf(rec_select_inpf,'%d',[1 inf]);
    fclose(rec_select_inpf);
else
    error('Wrong process.select option!');
end
clear dist rec_ind

% Throw other geophones' data
rec_throw=setdiff(rec_id,rec_process);
for i_throw=1:inf
	[tf,ichan_throw]=ismember(rec_throw,[DATA.geonum]); % This gives only the maximum index therefoere I need to do recursively
    
	if isempty(ichan_throw) || max(ichan_throw)==0
        break;
	end
    chan_throw=nonzeros(ichan_throw)';    
	
    DATA(chan_throw)=[];
end
DATA_NEW=DATA;
