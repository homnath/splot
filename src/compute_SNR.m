function [SNR, lta, sta] = compute_SNR(data,nsamp,lwin,lwlta,lwsta)
%computation of cumulative sum over all channels of a geophone
cuss=cumsum(sum(abs(data),2));
sta=(cuss(lwin+1:nsamp)-cuss(lwin+1-lwsta:nsamp-lwsta))/lwsta;
lta=(cuss(lwlta+1:nsamp-lwsta)-cuss(1:nsamp-lwin))/lwlta;
zz = find(lta==0);
lta(zz)=1;
sta(zz)=1;
sn=sta./lta;
% padding
sta=[zeros(lwin,1);sta];
lta=[zeros(lwin,1);lta];
SNR=[ones(lwin,1);sn];