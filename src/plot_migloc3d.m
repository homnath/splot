function plot_migloc3d(infname,varargin)
% This program plots seismograms, envelopes and SNR profiles with and without
% onsets, and polarization analysis. See fo example, self-explanatory input 
% file: plot_clair.in for detail.
%HISTORY
%   Jan 16,2015: added option difference plot, fixed some warnings.
%   Dec 01,2009: added option to normalize first to trace mean and then to absolute maximum
%   Nov 13,2009: added option for rotation with respect to given source location
%   Nov 04,2009: added option for onset plot given the source location
%   Nov 03,2009: added option for polarization plot
% TODO
%   Add option to plot all components in a single figure, e.g. in the order
%   1E, 1N, 1Z, 2E, 2N, 2Z .....
fprintf(1,'WARNING: all length units must be in m!\n');
iscache=0;
if nargin>1
    iscache=varargin{1};
end
if iscache
    fprintf(1,'WARNING: plotting from the catch file!\n');
end
if ~exist(infname,'file')
    error('input file ''%s'' not found!',infname);
end
% Read plot input file
[preinfo, rdata, process, seisplot, deplot, model] = read_plot(infname);

if preinfo.save>0
    figform=lower(strtrim(preinfo.figform));    
    if strcmp(figform,'png')
        preinfo.figform='-dpng';
    elseif strcmp(figform,'eps')
        preinfo.figform='-depsc';
    elseif strcmp(figform,'jpeg') || strcmp(figform,'jpg')
        preinfo.figform='-djpeg';
    elseif strcmp(figform,'ps')
        preinfo.figform='-dpsc';
    elseif strcmp(figform,'pdf')
        preinfo.figform='-dpdf';
    else
        error('%s format not supported! supported formats: png, eps, jpeg, ps, pdf',preinfo.figform);
    end
    addpath('./savefigure');
    [path]=fileparts(preinfo.fighead);
    if ~exist(path,'dir')
        fprintf('WARNING: directory ''%s'' doesn''t exist! thus created\n',fighead);
        mkdir(path);
    end    
end

% Plot seismograms
if preinfo.type==0 
    % Load data
    fprintf(1,'loading data...');
    [DATA]=load_data(rdata,iscache);
    fprintf(1,'complete!\n');
    
    % Select receivers to process
    if process.select ~= 0 || process.xrec ~= 0
        fprintf(1,'selecting receivers...');
        [DATA]=select_receivers(process,DATA);
        fprintf(1,'complete!\n');
    end
    
    % Plot
    plot_seismo(preinfo,model,process,seisplot,DATA);
elseif preinfo.type==1
    plot_deresult(preinfo,deplot);
else
    error('Wrong value %d for preinfo.type! valid values are (0 and 1)',preinfo.type);
end

if preinfo.save>0
    rmpath('./savefigure');
end