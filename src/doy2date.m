function [dd,mm] = doy2date(doy,year)
% VO function to calculate from day of year to month and day, including
% lead years.
% works for leap years from about 1582...
if ((mod(year,4)==0 && mod(year,100)~=0) || mod(year,400)==0)
  %leap year!
   dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
else
  % no leap year!
  dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
end

%if ismember(year,[1988 1992 1996 2000 2004 2008 2012])
%  dayspermonth = [0 31 29 31 30 31 30 31 31 30 31 30 31];
%else
%  dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
%end
mm = ones(size(doy))*NaN;
dd = ones(size(doy))*NaN;

for im = 1:12
   I = find(doy > sum(dayspermonth(1:im)) & doy <= sum(dayspermonth(1:im+1)));
   mm(I) = ones(size(I)).*im;
   dd(I) = doy(I) - ones(size(I))*sum(dayspermonth(1:im));
end

