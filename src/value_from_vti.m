function [value] = value_from_vti(file_name,ind,n)
% This function read value/s from a .vti binary file
% INPUT:
%   file_name   : file name [string]
%   ind         : indice/s to read specific value/s(if empty read all values)
%   n           : total number of values to read 
% History
%   HNG,Oct 15,2009;HNG,Oct 09,2009;HNG,May 26,2009

persistent pos % Position in the lookup table from where the data begins

inpf = fopen(file_name,'r');
% Read header and determine pos only once
if isempty(pos)    
    for i=1:inf
        sline=fgetl(inpf);  
        if size(strfind(sline, '<AppendedData'))>0
            break;
        end
    end

    fread(inpf,1,'uchar=>char'); % There is '_' in the beginning of data
    fread(inpf,1,'*int32'); % offset =
    pos = ftell(inpf); % file position
end
nind=length(ind);

% Read entire data
if nind == 0
    fseek(inpf,pos,'bof'); %float32 = 4 bytes    
    value = fread(inpf,n,'*float32');
    fclose(inpf);
    return;
end

if nind ~= n
    error('length of the indices /= specified size!');
end

% Read data only specified by the ind
value=zeros(nind,1);
for i_ind = 1:nind
    fseek(inpf,pos+4*(ind(i_ind)-1),'bof'); %float32 = 4 bytes    
    value(i_ind) = fread(inpf,1,'*float32');
end

fclose(inpf);
return;

