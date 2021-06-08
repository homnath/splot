function fig_hand = plot_color(x,y,c,clim,mark)
% Modified from plot3k (Ken Garrard, North Carolina State University, 2005, 2008
% based on plot3c by Uli Theune, University of Alberta)
% PLOT_COLOR color coded 3D scatterplot
%    fig_hand = plot3_color(x,y,c,clim,mark)
% Input
%    x, y    location vectors%              
%    c          color vector, default = 3rd column of pmat
%    clim       range of color indices     [min max], [max min] or []
%                                         default = [min(c) max(c)]
%    mark      marker type and size, default = {'.',6}
% Output
%    fig_hand         figure handle for the plot
%
% Generate a 2D scatter plot of pmat=(x,y) using the values in the vector c to
% determine the color of each point.  If c is empty, then z (column 3 of pmat)
% is used to color the plot.  The data points are sorted so that plot3 is
% only called once for each group of points that map to the same color.  The
% upper and lower limits of the color range (and the z axis) can be defined
% with clim.  This is useful for creating a series of plots with the same
% coloring.  The figure handle is returned if an output argument is given.  
% HISTORY
%   Dec 03,2009: Hom Nath Gharti

if length(x) ~= length(c)
   error('Location vector and color vector must be the same length');
end

% validate marker character
if ~iscell(mark), mark_type = mark;  mark_size = 6;
else
    mark_type = mark{1};
    if length(mark)>1, mark_size = mark{2};
    else               mark_size = 6;
    end
end
mark_str = '+o*.xsd^v><ph';
if (length(mark_type) ~= 1) || isempty(strfind(mark_str,mark_type))
   error('Invalid marker character, select one of ''%cind_first''', mark_str);
end
mark_size = max(mark_size,0.5);

% find color limits and range
if nargin < 3 || isempty(clim) % auto
   cmin   = min(c);
   cmax   = max(c);
else % user specified
   cmin = clim(1);
   cmax = clim(2);
end
crange = cmax - cmin;

pmat = [x(:) y(:)];

% Set colormap
cmap = colormap(jet(256));
if crange < 0
    cmap = flipud(cmap);
end
clen = length(cmap);

% calculate color map index for each point
pmat(:,3) = min(max(round((c-min(cmin,cmax))*(clen-1)/abs(crange)),1),clen);

% sort by color map index
pmat = sortrows(pmat,3);

% build vector of last indices corresponding to each unique color index (last point for each color)
cind_last = [find(diff(pmat(:,3))>0); length(pmat)];

% plot data points in groups by color map index
hold on;                           % add points to one set of axes
cind_first = 1; % index of 1st point in a  color group
for i_cind = 1:length(cind_last) % loop over each non-empty color group
   plot(pmat(cind_first:cind_last(i_cind),1), ... % call plot3 once for each color group
        pmat(cind_first:cind_last(i_cind),2), ...
        mark_type,           ...
        'MarkerSize',mark_size, ...
        'MarkerEdgeColor',cmap(pmat(cind_first,3),:), ... % same marker color from cmap
        'MarkerFaceColor',cmap(pmat(cind_first,3),:) ); % for all points in group
   cind_first = cind_last(i_cind)+1;                  % next group starts at next point
end
hold off;

if nargout > 0
    fig_hand = gcf; % return figure handle
end
