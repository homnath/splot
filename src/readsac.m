function [data,varargout] = readsac(sacfile)
%sac2mat  loads a SAC2000 file into Matlab, along with header information.
%  [data] = readsac(sacfile);
%  [data,header] = readsac(sacfile);
%  [data,header,kstnm] = readsac(sacfile);
% 3rd argument is optional
%
%  sacfile must be a STRING, eg.  [fff,hhh,name]=sac2mat('dir1/file1.sac')
%
%Output:
%data=	the data from the file
%header=the first 105 header variables, as described on the SAC Home Page
%	http://www.llnl.gov/sac/
%	http://www.llnl.gov/sac/SAC_Manuals/FileFormatPt1.html
%kstnm=	KSTNM, the station name (string)
%	(header and kstnm are optional)
%
%Some common/useful headers:
%header(1:3) DELTA,DEPMIN,DEPMAX
%header(57) DEPMEN
%header(6:7) B,E
%header(8) Event origin time (seconds relative to reference time)
%header(32:39) STLA,STLO,STEL,STDP,EVLA,EVLO,EVEL,EVDP
%header(51:54) DIST,AZ,BAZ,GCARC        
%header(71:76)  NZYEAR,NZDAY,NZHOUR,NZMIN,NZSEC,NZMSEC
%header(80) NPTS
%						Case Bradford, 8 June 2004
% Change history: Around 2004, Volker Oye

sacfid = fopen(sacfile,'r','n');
if sacfid==-1;
    %msg=sprintf('Invalid path: %s',sacfile);
    error('Invalid path: %s',sacfile);
    return;
end

%The first 70 header variables are floating point
header1 = fread(sacfid,70,'float32');

%The next 35 are integers (after that they're a mix of strings)
header2 = fread(sacfid,35,'int32');

%Outputs the number of points being read to the screen
header=[header1;header2];
%sprintf('%6f',header(1))
npts = header(80);
%sprintf('Hello everybody %d\n',npts)
%kstnm, the station name, is a string stored in the 110th header variable
%kstnm = KSTNM (First 3 letters of the station name)
fseek(sacfid,4*110,-1);
if nargout > 1
    varargout(1)={header};
end
if nargout > 2
    varargout(2) = {char(fread(sacfid,4)')}; % Station name is optional
end

%Now read the data...
fseek(sacfid,158*4,-1);
data(1:npts) = fread(sacfid,npts,'float32');
fclose(sacfid);
return