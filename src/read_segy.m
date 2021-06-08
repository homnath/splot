
function [TRACE_HEAD,TRACES] = read_segy(filename,trace_num,TRACE_LENGTH,EXTRA_BYTES,NR_PTS,FORM_CODE,DATA_FORMAT)


% Lj 25SEP 2001 implementing fortran format 
% Scan the header

%[ASC_HEAD,BIN_HEAD,NUM_TRACES,NR_PTS,SAMPLE_RATE,DATA_FORMAT,FORM_CODE,TRACE_LENGTH,EXTRA_BYTES] = segy_header(filename);

[fid,message] = fopen(filename,'r','ieee-be');

%Preallocation of RAM, speeds up factor more than 10 plus avoids OUT OF
%MEMORY error message with big files.
%t=cputime;
TRACE_HEAD=zeros(94,length(trace_num));
TRACES=zeros(NR_PTS,length(trace_num));
for n=1:length(trace_num)
    % potential extra bytes are also included in variable TRACE_LENGTH
	position = 3200+400 + (5*EXTRA_BYTES/2) + (trace_num(n)-1)*TRACE_LENGTH;
	status=fseek(fid,position,'bof');
	
	TRACE_HEAD(1:7,n)  = fread(fid,7, 'int32');
	TRACE_HEAD(8:11,n) = fread(fid,4, 'int16');
	TRACE_HEAD(12:19,n)= fread(fid,8, 'int32');
	TRACE_HEAD(20:21,n)= fread(fid,2, 'int16');
	TRACE_HEAD(22:25,n)= fread(fid,4, 'int32');
	TRACE_HEAD(26:71,n)= fread(fid,46,'int16');
	TRACE_HEAD(72:77,n)= fread(fid,6, 'float32');
	TRACE_HEAD(78:78,n)= fread(fid,1, 'int32');
	TRACE_HEAD(79:94,n)= fread(fid,16,'int16');
 %   TRACES{n} = fread(fid,NR_PTS,FORM_CODE);
	TRACES(1:NR_PTS,n) = fread(fid,NR_PTS,FORM_CODE);
end

fclose(fid);
	
% if the data are in 4-byte IBM-floating point format they have to be converted into IEEE format  
if  DATA_FORMAT == 1
  for n=1:length(trace_num)
    TRACES(1:NR_PTS,n) = ibm2ieee(TRACES(1:NR_PTS,n));
  end
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = ibm2ieee (ibmf)
% Name:         ibm2ieee
% Abstract:     convert a matrix of IBM/360 32-bit floats
%               to IEEE doubles.
%
%               IBMF is the matrix of IBM/360 32-bit floats each
%               stored as a 32 bit unsigned big-endian integer
%               in a MATLAB double.
%
%               The format of a IBM/360 32-bit float is:
%
%                  sign 7-bit exponent  24 bit fraction
%                  The base is 16. The decimal point is to
%                  the left of the fraction. The exponent is
%                  biased by +64.
%
%               The basic idea is that you use floating point on
%               the various fields.
%
%               ieee = sign * 16 ^ (exponent - 64) * fraction / 16 ^ 6
%
% By:           Martin Knapp-Cordes
%               The MathWorks, Inc.
%
% Date(s):      Jun 95 - 28, 29

% $Revision: 1.2 $  $Date: 1995/06/29 14:50:03 $
%----------------------------------------------------------------------------
%
    if (nargin ~= 1)
        error ('Wrong number of arguments.');
    elseif (isempty(ibmf))
        error ('Argument is an empty matrix');
    end
%
    aibmf = sprintf('%08x',ibmf);
%
% hexd(1) - 1st hex digit - first bit is sign, next 3 bits high order exponent
% hexd(2) - 2nd hex digit - bits of low order exponent
% hexd(3) - 3rd-8th hex digit - bits of fraction
%
    hexd = sscanf(aibmf,'%1x%1x%6x',[3,inf]);
    d = (1 - (hexd(1,:) >= 8) .* 2) .* ...
        16 .^ ((hexd(1,:) - (hexd(1,:) >= 8) .* 8) .* 16 + hexd(2,:) ...
        - 70).* hexd(3,:);
    d = reshape(d,size(ibmf));



 
  function revword = byte_reversal( word )
  % return a byte - reversed  word
   
  for i=1:length(word)
    
    ib1 =  bitand( word(i)               ,   2^8-1);
    ib2 =  bitand( bitshift(word(i), -8) ,   2^8-1);
    ib3 =  bitand( bitshift(word(i),-16) ,   2^8-1);
    ib4 =  bitand( bitshift(word(i),-24) ,   2^8-1);
    
    revword(i) = bitor(  bitor(ib4 , bitshift(ib3,8) )  , bitor ( bitshift(ib2,16) , bitshift(ib1,24) )  );
  end
  






