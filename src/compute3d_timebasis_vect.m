% This fucntion computes the P and S arrival times, and basis vectors along P,
% SH and SV taking SH and SV directions as mentioned below in 3-dimensional
% homogeneous domain
function [ptime,stime,pbasis,shbasis,svbasis] = timebasis_vect(gnum,xs,ys,zs,src_t0,xr,yr,zr,pvel,svel)
global lookup_header_p lookup_header_s use_lookup
global ox oy oz dh nx ny nz
global tp_lookup ts_lookup
% Need to check why dot(pbasis,svbasis) is not always 0.0
% P is along the ray path
% SV is perpendicular to P and lies on the ZP plane
% SH is perpendiculer to both P and SV, and direction being congruent with
% the cross product of P and SV
% all units of dimensional variables are considered to be consistent as called by main routine 
% HISTORY: HNG,Oct08,2008; HNG,May26,2009
%TODO: make pbasis,shbasis,svbasis optional (required for no rotation)

pbasis=zeros(1,3);
shbasis=zeros(1,3);
svbasis=zeros(1,3);
zvect=zeros(1,3); zvect(3)=1.0; % Normal to horzontal plane

% If source and receiver coincide P, SH and SV directions are taken as just
% the original axis
pbasis(1)=1.0;
shbasis(2)=1.0;
svbasis(3)=1.0;

ptime=0; stime=0;

% P vector not yet basis vector though
pbasis(1)=(xr-xs);
pbasis(2)=(yr-ys);
pbasis(3)=(zr-zs);

pnorm=norm(pbasis);

if pnorm ~= 0.0
    ptime=pnorm/pvel+src_t0;
    stime=pnorm/svel+src_t0;    
    
    shbasis=cross(pbasis,zvect);    
    svbasis=cross(pbasis,shbasis);    
    % We need unit vectors   
    pbasis=pbasis/pnorm;
    shbasis=shbasis/norm(shbasis); 
    svbasis=svbasis/norm(svbasis);       
end

% Use reciprocity lookup table for arrival times
% Note: we approximate the ray direction assuming homogeneous model
if use_lookup>0
    file_p = sprintf('%s%d.vti',lookup_header_p,gnum-1);
    file_s = sprintf('%s%d.vti',lookup_header_s,gnum-1);
    % compute indices of 8 corners
    frac_ix = (xs-ox)/dh; 
    frac_iy = (ys-oy)/dh; 
    frac_iz = (-zs-oz)/dh; %sign of zs was changed when passed into this function 
    
    ix(1) = floor(frac_ix)+1; 
    iy(1) = floor(frac_iy)+1; 
    iz(1) = floor(frac_iz)+1; %sign of zs was changed when passed into this function 
    
    if frac_ix == floor(frac_ix) && frac_iy == floor(frac_iy) && frac_iz == floor(frac_iz)
        % point exactly coincides with the grid point
        ind = (iz(1)-1)*ny*nx+(iy(1)-1)*nx+ix(1);
        
        if use_lookup == 1
            % Read particular value/s from the file/s
            tp_src = value_from_vti(file_p,ind,1);
            ts_src = value_from_vti(file_s,ind,1);
        else
            % We have already loaded the lookup table
            tp_src = tp_lookup(gnum,ind);
            ts_src = ts_lookup(gnum,ind);
        end
        
        ptime = tp_src+src_t0;
        stime = ts_src+src_t0;
        if isnan(ptime) || isnan(stime)
            disp(ptime);
            disp(stime);
            disp(tp_src);
            disp(ts_src);
            disp(ind);
            disp(ix(1));
            disp(iy(1));
            disp(iz(1));            
            error('Arrival time is NaN!');
        end
        return
    end
    
    % use forward points except for furthest boundary nodes
    if ix(1) == nx
        ix(2)=ix(1)-1;
        [ix(1), ix(2)]=deal(ix(2),ix(1)); % Swap
    else
        ix(2)=ix(1)+1;
    end
    if iy(1) == ny
        iy(2)=iy(1)-1;
        [iy(1), iy(2)]=deal(iy(2),iy(1)); % Swap
    else
        iy(2)=iy(1)+1;
    end
    if iz(1) == nz
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
        x(i_x)=ox+(ix(i_x)-1)*dh;
        y(i_x)=oy+(iy(i_x)-1)*dh;
        z(i_x)=oz+(iz(i_x)-1)*dh;
    end
    [X,Y,Z]=meshgrid(x,y,z);
    
    inum=0;
    ind=zeros(8,1); % indices for 8 corners of a cube
    for i_z=1:2 
        for i_y=1:2 
            for i_x=1:2                
                inum=inum+1;                
                ind(inum) = (iz(i_z)-1)*ny*nx+(iy(i_y)-1)*nx+ix(i_x);   
            end
        end
    end
    
    if max(ind)>nx*ny*nz
        error('Offset is ouside the range!');
    end
    if use_lookup == 1
        % Read particular value/s from the file/s
        tp = value_from_vti(file_p,ind,8);
        ts = value_from_vti(file_s,ind,8);
    else
        % We have already loaded the lookup table
        tp = tp_lookup(gnum,ind);
        ts = ts_lookup(gnum,ind);
    end        
    
    % Change the shape of the array for interp3, it needs 3D array 
    tp=reshape(tp,[2 2 2]);
    ts=reshape(ts,[2 2 2]);
    
    % Interpolation
    tp_src = interp3(X,Y,Z,tp,xs,ys,-zs);
    ts_src = interp3(X,Y,Z,ts,xs,ys,-zs);
    
    ptime = tp_src+src_t0;
    stime = ts_src+src_t0;
    if isnan(ptime) || isnan(stime)
        disp(tp);
        disp(ts);
        disp(tp_src);
        disp(ts_src);
        disp(x);
        disp(y);
        disp(z);
        disp(xs);
        disp(ys);
        disp(zs);
        error('Arrival time is NaN!');        
    end
end % if use_lookup>0

return