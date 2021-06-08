%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ASC_HEAD,BIN_HEAD,NUM_TRACES,NR_PTS,SAMPLE_INTV,DATA_FORMAT,FORM_CODE,TRACE_LENGTH,EXTRA_BYTES]=segy_header(filename)
% This function will read the EBCDIC and BINARY header of a segy file
% Automatic detection of blocked format
% input : 
%         filename : Segy file name [complete]
%          
% output:
%         ASC_HEAD    : EBCDIC HEADER 3200 bytes [ Ascii format assumed , no conversions]
%         BIN_HEAD    : The Binary header 400 bytes
%         NUM_TRACES  : Total number of traces on the file
%         NR_PTS      : Trace length in words NOT including 60 word header
%         SAMPLE_INTV : Sample rate in seconds
%         DATA_FORMAT : Format code indicating IBMIEEE for instance
%         FORM_CODE   : int32 , int16 , uint or float32
%         TRACE_LENGTH: Trace length in bytes 
%         EXTRA_BYTES : In case of fortran unformatted ther will be 8 extra bytes per record
%
%  Michael Roth 10.02.03 (modified)
%  Leif Jahren 25.09.01 (original)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
function [ASC_HEAD,BIN_HEAD,NUM_TRACES,NR_PTS,SAMPLE_INTV,DATA_FORMAT,FORM_CODE,TRACE_LENGTH,EXTRA_BYTES]=segy_header(filename)
 
 
% In the SEGY standard all binary values are defined as using big-endian byte ordering  
[fid,message] = fopen(filename,'r','ieee-be');

% the following check is for the Read Well (RW) blocked format
fortran_asctest = fread(fid,1 ,'int32');
if fortran_asctest==3200
  EXTRA_BYTES = 8;
  ASC_HEAD        = fread(fid,3200  ,'uchar');
  fortran_bintest = fread(fid,2     ,'int32');
else 
  EXTRA_BYTES = 0;
  frewind(fid);
  ASC_HEAD = fread(fid,3200  ,'uchar');
end
 

%read the binary header
BIN_HEAD(1:3)   = fread(fid,  3,'uint32');
BIN_HEAD(4:197) = fread(fid,194,'uint16');


%get the important parameter from the binary header 
SAMPLE_INTV = BIN_HEAD(6)/1000/1000;
NR_PTS      = BIN_HEAD(8);
DATA_FORMAT = BIN_HEAD(10);
  
%get the word size and format of the data  
fact = 4;

if DATA_FORMAT == 1
  %data in 4-byte IBM floating-point format 
  %in order to convert the floating point format into ieee-floating point
  %format the data must be read as unsigned integers and forwarded to the 
  %routine ibm2ieee
  form_code = 'uint';       
elseif DATA_FORMAT == 2
  %data in 4-byte two's complement integer (not tested)  
  form_code = 'int32';
elseif DATA_FORMAT == 3
   %data in 2-byte two's complement integer (not tested)   
  form_code = 'int16'; 
  fact=2;
elseif DATA_FORMAT == 5
  %data in 4-byte IEEE floating-point format  
  form_code = 'float32';
end 
  

% lj 27 SEP  discrepancy between DATA_FORMAT and actual storage format
% in case of RW blocked format the data are stored as 4-byte IEEE floating-point 
if fortran_asctest==3200
  form_code = 'float32';
  % RWS blocked data are formatted as IEEE floating-point, but the DATA_FORMAT is set to 1 in their files.
  % Craig Woerpel 15.07.30
  DATA_FORMAT = 5;
end
  
 
%get the file size in bytes
fseek(fid,0,'eof'); FILE_SIZE=ftell(fid);
 
 
  
%  FILE_SIZE = getfield(dir(filename),'bytes');
%compute the number of traces 
NUM_TRACES=(FILE_SIZE-(3200 +400 +2*EXTRA_BYTES))/(240+NR_PTS*fact+EXTRA_BYTES);
FORM_CODE = form_code;
TRACE_LENGTH = 240+EXTRA_BYTES+NR_PTS*fact;
fclose(fid);

 