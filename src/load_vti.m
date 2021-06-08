function [ox,dh,ns,array,data_type] = load_vti(infname)
% ox: origin vector (x,y,z)
% dh: sampling interval vector (dx,dy,dz)
% ns: number of samples vector (nx,ny,nz)
% array: values stored in 1D array of size nx*ny*nz in the order of x index
% changes first
% array_name: name of the array
inpf=fopen(infname,'r');
%headerline1 = fgetl(inpf);
%headerline2 = fgetl(inpf);
%source_location = fgetl(inpf);
%headerline4 = fgetl(inpf);
headerline5 = fgetl(inpf);
headerline6 = fgetl(inpf);
Endian = headerline6(findstr(headerline6,'byte_order="')+12:end-2);
if strcmp(lower(Endian),'littleendian');
  Endian='l';
else
  Endian='b';
end
bulk=textscan(inpf,'<%*s WholeExtent="%d %d %d %d %d %d" Origin="%f %f %f" Spacing="%f %f %f">');
fgetl(inpf);
data_type = textscan(inpf,'<%*s Scalars=%q');
data_type = char(data_type{1}); % variable name        

fclose(inpf); 

wext=[bulk{1} bulk{2} bulk{3} bulk{4} bulk{5} bulk{6}]; % Whole extent
ox=[bulk{7} bulk{8} bulk{9}]; % Origin
dh=[bulk{10} bulk{11} bulk{12}]; % Spacing
ns(1)=wext(2)-wext(1)+1;
ns(2)=wext(4)-wext(3)+1;
ns(3)=wext(6)-wext(5)+1;
ns=double(ns);

% Read lookup table values only for indices in 'ind', here for all points
% (prod(ns)).
array = value_from_file(infname,prod(ns),Endian);

return

function [value] = value_from_file(file_name,n,Endian)
% This function read a value specified by
% file_name: file name [string]
% skip_nline: number of lines to skip [int]
% skip_nchar: number of characters to skip [int]
% ind: number of values to skip (float32 values) [int]
% Use
% Endian: checks whether big of little endian data

%   value_from_file(file_name,skip_nline,skip_nchar,ind)
% History
%   HNG,May 26,2009
inpf = fopen(file_name,'r',Endian);

% Read header and determine pos   
% Skip lines
% for i_skip = 1:skip_nline
%     fgetl(fin);
% end

% Skip characters
% for i_skip = 1:skip_nchar
%     c=fread(fin,1,'uchar=>char');
% end
for i=1:inf
    sline=fgetl(inpf);  
    if size(strfind(sline, '<AppendedData'))>0
        break;
    end
end

c=fread(inpf,1,'uchar=>char'); % There is '_' in the beginning of data
fread(inpf,1,'*int32'); % offset =

%fread(fin,1,'*int32');
value = fread(inpf,n,'*float32');
fclose(inpf);

return;
