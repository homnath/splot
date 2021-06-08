%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [stat,fdb,tdb,data]=read_seg2 ('filename')
% reads a SEG-2 (standard SEG-2) file from disk.
% input  :
%          filename (character, absolute filename)
% output : 
%           stat (integer, status 1 means ok, 0 means not ok
%           fdb (character array, file descriptor block)
%           tdb (cell array with character arrays, trace descriptor blocks)
%           data (cell array with vectors, seismograms)
%   Michael Roth 04.02.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [istat,fdb,tdb,data] = read_seg2(filename)

istat=0;
fdb=[];tdb=[];data=[];



%open data file, SEG2 should be in IEEE little endien format 
fid=fopen (filename,'rb','ieee-le');
   

% check for SEG-2 file type
% first 2 bytes equal '3a55h' (14933) for PC/Windows
fileType=fread(fid,1,'ushort');
if (fileType ~= 14933)
  fclose (fid);
  str=strcat('routine read_seg2: not a standard SEG2 file:',filename);
  disp(str)
  return
end;

try
  % read file descriptor block
  revNumber          = fread (fid,1,'short');
  sizeOfTracePointer = fread (fid,1,'ushort');
  ntrace             = fread (fid,1,'ushort');
  sizeOfST           = fread (fid,1,'uchar');
  firstST            = fread (fid,1,'char');
  secondST           = fread (fid,1,'char');
  sizeOfLT           = fread (fid,1,'uchar');
  firstLT            = fread (fid,1,'char');
  secondLT           = fread (fid,1,'char');
  reserved           = fread (fid,18,'uchar');
  tracePointers      = fread (fid,ntrace,'ulong');
  
  % read free strings 
  iflag=0;
  fseek (fid,32+sizeOfTracePointer,'bof');
  offset = fread(fid,1,'ushort');
  while (offset > 0)
    %freeString = setstr(fread (fid,offset-2,'char'))';
    freeString = fread (fid,offset-2,'*char')';
    if (iflag==0);   
      fdb=freeString;
      iflag=1;
    else
      fdb=char(fdb,freeString);
    end; 
    offset = fread(fid,1,'ushort');
  end;

  % read data 
  % vo pre-allocate RAM-space for data
  data = cell(1,ntrace);
  for i=1:ntrace;
  % read trace descriptor block
    tracePointers(i);
    fseek (fid,tracePointers(i),'bof');
    traceId     = fread (fid,1,'ushort');
    sizeOfBlock = fread (fid,1,'ushort');
    sizeOfData  = fread (fid,1,'ulong');
    nsamp       = fread (fid,1,'ulong');
    dataCode    = fread (fid,1,'uchar');
    reserved    = fread (fid,19,'uchar');
    iflag=0;
    offset = fread (fid,1,'ushort');
    while (offset > 0);
      %freeString = setstr(fread (fid,offset-2,'char'))';
      freeString = fread (fid,offset-2,'*char')';
      if (iflag==0);   
        temp=freeString;
        iflag=1;
      else
        temp=char(temp,freeString);
      end; 
      offset = fread(fid,1,'ushort');
    end
    tdb{i}=temp;
    fseek (fid,tracePointers(i)+sizeOfBlock,'bof');
    % read data block    
    if(dataCode==1)
      data{i}=fread(fid,nsamp,'int16');
    elseif(dataCode==2) 
      data{i}=fread(fid,nsamp,'int32');
    elseif(dataCode==4)
      data{i}=fread(fid,nsamp,'float32');
    elseif(dataCode==5)
      data{i}=fread(fid,nsamp,'float64');
    end 
  end;
  istat=1;
  fclose (fid);  
catch
  str=strcat('routine read_seg2: error',lasterr);    
  disp(str);
  fclose (fid);
end
