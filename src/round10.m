function rnum=round10(num,varargin)
% This function round the number to nearest pth power of 10
% This is similar to the Mapping Toolbox function roundn
% Hom Nath Gharti, Sep 02,2009

if nargin == 1
    p=1; %Default power of 10
elseif nargin == 2
    p=varargin{1};
else
    error 'wrong number of arguments!'
end

fact=10^fix(-p);
rnum=round(fact*num)/fact;
% nsign=sign(num); % Sign of the number
% 
% num=abs(round(num)); % Round the number to nearest integer
% bnum=10^round(p); % Base number
% r=mod(num,bnum); % Reminder
% if r < bnum/2
%     rnum=num-r; % Rounded number
% else
%     rnum=num+(bnum-r); % Rounded number
% end
% rnum=nsign*rnum;
