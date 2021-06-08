function [value] = value_from_file(file_name,skip_nline,skip_nchar,ind,n)
global pos
% This function read a value specified by
% file_name: file name [string]
% skip_nline: number of lines to skip [int]
% skip_nchar: number of characters to skip [int]
% ind: number of values to skip (float32 values) [int]
% Use
%   value_from_file(file_name,skip_nline,skip_nchar,ind)
% History
%   HNG,May 26,2009
fin = fopen(file_name,'r');

% Read header and determine pos only once
if isnan(pos)    
    % Skip lines
    for i_skip = 1:skip_nline
        fgetl(fin);
    end

    % Skip characters
    for i_skip = 1:skip_nchar
        fread(fin,1,'uchar=>char'); %dum =
        %disp(dum);
    end
    fread(fin,1,'*int32'); % offset =
    pos = ftell(fin); % file position
end
nind=length(ind);

% Read entire data
if nind == 0
    fseek(fin,pos,'bof'); %float32 = 4 bytes    
    value = fread(fin,n,'*float32');
    fclose(fin);
    return;
end

if nind ~= n
    error('length of the indices /= size specified!');
end

% Read data only specified by the ind
value=zeros(nind,1);
for i_ind = 1:nind
    fseek(fin,pos+4*(ind(i_ind)-1),'bof'); %float32 = 4 bytes    
    value(i_ind) = fread(fin,1,'*float32');
end

fclose(fin);
return;

