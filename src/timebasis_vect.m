function [tp_src,ts_src,pbasis,shbasis,svbasis] = timebasis_vect(model,rec_id,recx,srcx)
% This fucntion computes the P and S arrival times for homogeneous medium 
% or extract the arrival times from lookup ables. This function also computes 
% the basis vectors along P,SH and SV (L, T and Q) in 3-dimensional homogeneous 
% domain. For the simplicty, SH (T) axis is always taken on the horizontal plane.
% P is along the ray path, SV is perpendicular to P and lies on the ZP plane, 
% and SH is perpendiculer to both P and SV, and direction being congruent with
% the cross product of P and SV all units of dimensional variables are considered
% to be consistent as called by main routine. The basis vectors returned will be
% in the Z axis up convention.
% INPUT:
%   grid: this structures contains origin (ox, oy, z), sampling interval
%           (dh), and number of grid points (nx, ny, nz)
%   rec_id: live receivers ID (number)
%   srcx: Coordinates of trial source in z down convention
%   recx: Coordinates of receivers in z down convention
%   model.vp: P wave velocity for homogeneous model
%   model.vs: S wave velocity for homogeneous model
%   lookup: this structure contains pheader, sheader, option, and optinally tp and ts
% OUTPUT:
%   tp_src: P wave travel time
%   ts_src: S wave travel time
%   pbasis: P basis vector
%   shbasis: SH basis vector
%   svbasis: SV basis vector
% HISTORY: HNG,Oct 15,2009;HNG,Oct08,2008; HNG,May26,2009
%TODO:
% - make pbasis,shbasis,svbasis optional (required for no rotation)
% - Need to check why dot(pbasis,svbasis) is not always 0.0

nrec=length(rec_id); % Number of live receivers

pbasis=zeros(nrec,3);
shbasis=zeros(nrec,3);
svbasis=zeros(nrec,3);

zvect=zeros(1,3); zvect(3)=1.0; % Normal to horizontal plane and upward

% If source and receiver coincide P, SH and SV directions are taken as just
% the original axis
pbasis(:,1)=1.0;
shbasis(:,2)=1.0;
svbasis(:,3)=1.0;

tp_src=zeros(nrec,1); ts_src=zeros(nrec,1);

% P vector not yet basis vector though
pbasis(:,1)=recx(:,1)-srcx(1);
pbasis(:,2)=recx(:,2)-srcx(2);
pbasis(:,3)=-recx(:,3)+srcx(3); % Z axis is converted to positive upward only here

pnorm=sqrt(dot(pbasis,pbasis,2));

irec_comp=find(pnorm~=0.0); % Indices of receivers for which following computation is necessary
nrec_comp=length(irec_comp); % Number of receivers for which following computation is necessary
% Basis vectors are computed assuming homogeneous isotropic model
for i_rec=1:nrec_comp
    irec=irec_comp(i_rec);       
    
    shbasis(irec,:)=cross(pbasis(irec,:),zvect);
    if dot(shbasis(irec,:),shbasis(irec,:))==0
        shbasis(irec,:)=[1 0 0]; % Take x-axes as shbasis
    end
    svbasis(irec,:)=cross(pbasis(irec,:),shbasis(irec,:));    
    % We need unit vectors   
    pbasis(irec,:)=pbasis(irec,:)/pnorm(irec);    
    shbasis(irec,:)=shbasis(irec,:)/sqrt(dot(shbasis(irec,:),shbasis(irec,:)));    
    svbasis(irec,:)=svbasis(irec,:)/sqrt(dot(svbasis(irec,:),svbasis(irec,:)));    
end
if model.lookup<0
    % no need to compute arrival times
    return;
end
if model.lookup==0
    tp_src(irec_comp)=pnorm(irec_comp)/model.vp;
    ts_src(irec_comp)=pnorm(irec_comp)/model.vs;
elseif model.lookup==1
    % Use reciprocity lookup table for arrival times
    % compute indices of 8 corners
    frac_ix = (srcx(1)-model.ox(1))/model.dh; 
    frac_iy = (srcx(2)-model.ox(2))/model.dh; 
    frac_iz = (srcx(3)-model.ox(3))/model.dh; 
    
    ix(1) = floor(frac_ix)+1; 
    iy(1) = floor(frac_iy)+1; 
    iz(1) = floor(frac_iz)+1; 
    
    if frac_ix == floor(frac_ix) && frac_iy == floor(frac_iy) && frac_iz == floor(frac_iz)
        % point exactly coincides with the grid point
        ind = (iz(1)-1)*model.nx(2)*model.nx(1)+(iy(1)-1)*model.nx(1)+ix(1);
        
        if model.lookup == 1
            for i_rec=1:nrec_comp
                irec = irec_comp(i_rec); % receiver number to compute
                % P times
                file_p = sprintf('%s%d.vti',model.pheader,rec_id(irec)-1); % rec_id(irec) is actual receiver number
                tp_src(irec) = value_from_vti(file_p,ind,1);
                % S times
                file_s = sprintf('%s%d.vti',model.sheader,rec_id(irec)-1);                                
                ts_src(irec) = value_from_vti(file_s,ind,1);
            end
        else
            % We have already loaded the lookup table
            tp_src(irec_comp) = model.tp(irec_comp,ind);
            ts_src(irec_comp) = model.ts(irec_comp,ind);
        end
%         if isnan(ptime) | isnan(stime) 
%             error('Arrival time is NaN!');
%         end
        return
    end % if frac_ix
    
    tp=zeros(nrec_comp,8);
    ts=zeros(nrec_comp,8);
    % use forward points except for furthest boundary nodes
    if ix(1) == model.nx(1)
        ix(2)=ix(1)-1;
        [ix(1), ix(2)]=deal(ix(2),ix(1)); % Swap
    else
        ix(2)=ix(1)+1;
    end
    if iy(1) == model.nx(2)
        iy(2)=iy(1)-1;
        [iy(1), iy(2)]=deal(iy(2),iy(1)); % Swap
    else
        iy(2)=iy(1)+1;
    end
    if iz(1) == model.nx(3)
        iz(2)=iz(1)-1;
        [iz(1), iz(2)]=deal(iz(2),iz(1)); % Swap
    else
        iz(2)=iz(1)+1;
    end
    
    % ix,iy and iz are indexed from 1 not 0    
    x=zeros(2,1); % 2 points in a first side
    y=zeros(2,1); % 2 points in a first side
    z=zeros(2,1); % 2 points in a first side
    
    for i_x=1:2 %i_x = i_y = i_z
        x(i_x)=model.ox(1)+(ix(i_x)-1)*model.dh;
        y(i_x)=model.ox(2)+(iy(i_x)-1)*model.dh;
        z(i_x)=model.ox(3)+(iz(i_x)-1)*model.dh;
    end
    [X,Y,Z]=meshgrid(x,y,z);
    
    inum=0;
    ind=zeros(8,1); % indices for 8 corners of a cube
    for i_z=1:2 
        for i_y=1:2 
            for i_x=1:2                
                inum=inum+1;                
                ind(inum) = (iz(i_z)-1)*model.nx(2)*model.nx(1)+(iy(i_y)-1)*model.nx(1)+ix(i_x);   
            end
        end
    end
    
    if max(ind)>model.nx(1)*model.nx(2)*model.nx(3)
        error('Offset is ouside the range!');
    end
    if model.lookup == 1
        for i_rec=1:nrec_comp
            irec = irec_comp(i_rec); % Actual receiver number
            % P times
            file_p = sprintf('%s%d.vti',model.pheader,rec_id(irec)-1);
            tp(i_rec,:) = value_from_vti(file_p,ind,8);
            % S times
            file_s = sprintf('%s%d.vti',model.sheader,rec_id(irec)-1);            
            ts(i_rec,:) = value_from_vti(file_s,ind,8);
        end
    else
        % We have already loaded the lookup table
        tp(1:nrec_comp,:) = model.tp(irec_comp,ind);
        ts(1:nrec_comp,:) = model.ts(irec_comp,ind);
    end        
    for i_rec=1:nrec_comp
        % Change the shape of the array for interp3, it needs 3D array 
        tp_cube=reshape(tp(i_rec,:),[2 2 2]);
        ts_cube=reshape(ts(i_rec,:),[2 2 2]);
        
        irec=irec_comp(i_rec);
        % Interpolation
        tp_src(irec) = interp3(X,Y,Z,tp_cube,srcx(1),srcx(2),srcx(3));
        ts_src(irec) = interp3(X,Y,Z,ts_cube,srcx(1),srcx(2),srcx(3));
%         if isnan(ptime) | isnan(stime)
%             error('Arrival time is NaN!');        
%         end
    end
elseif model.lookup==2
    % read arrival times directly from the respective files
    % p arrival time
    inpf=fopen(model.pheader,'r');
    fgetl(inpf);
    fgetl(inpf);
    tp_tmp=fscanf(inpf,'%f',inf);
    tp_src(irec_comp)=tp_tmp(irec_comp);
    fclose(inpf);
    % s arrival time
    inpf=fopen(model.sheader,'r');
    fgetl(inpf);
    fgetl(inpf);
    ts_tmp=fscanf(inpf,'%f',inf);
    ts_src(irec_comp)=ts_tmp(irec_comp);
    fclose(inpf);
else
    error('wrong value for model.lookup: %d',model.lookup)
            
        
end % if model.lookup==0
%return