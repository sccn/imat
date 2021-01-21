% plot original IC time-frequency decomposition, PC time-frequency backprojection, IM time-frequency backprojection and IM
% weights
% 
% 
% pop_plotIMtimecourse(EEG,varargin)
%
% 
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from functions written by Julie Onton 
% 
% 
% Example: plot components 1, 5 and 8 with frequency limit 6-40Hz and smoothing 0.1
%           - plot only IC time-frequency decomposition and PC time-frequency backprojections 
%
% >>  pop_plotIMtimecourse(EEG, 'comps', [1 5 8],...
%      'frqlim', [6 40], 'smoothing', 0.1, 'plotICtf', 'on', 'plotPCtf', 'on')
%
%
% Example: plot components 3, 5 and 6 and IMs 2, 6 and 9 with frequency limit 6-40Hz and smoothing 0.1
%          plot only IM time-frequency backprojections on ICs and IM weight
%           timecourses
%
% >>  pop_plotIMtimecourse(EEG, 'comps', [3 5 6], 'factors', [2 6 9],...
%      'frqlim', [6 40], 'smoothing', 0.1, 'plotIMtf', 'on', 'plotIMtime', 'on')
%
%
% INPUTS:
% required Inputs:
% EEG - EEG structure pointing to precomputed IMA information and EEG
%         datasets
%
% optional Inputs
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject 
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject% frqlim - [min max] frequency limits for plotting
% frqlim - [min max] frequency limits for plotting
% smoothing - [integer] how many timepoints to use for moving average for
%                        smoothing. A higher number means more smoothing. default 40
% plotICtf  -  [string]  {'on' 'off'} plot IC time-frequency decomposition, default 'off'                    
% plotPCtf  -   [string] {'on' 'off'}  plot PC time-frequency backprojection, default 'off'                      
% plotIMtf  -    [string] {'on' 'off'}  plot IM time-frequency backprojection, default 'off'                      
% plotIMtime  -  [string] {'on' 'off'}   plot IM weights timecourse, default 'off'                     


function pop_plotIMtimecourse(EEG,varargin)

if ~isfield(EEG.etc, 'IMA'), return; end
try
    tmpIMA = load([EEG.etc.IMA.filepath '/' EEG.etc.IMA.filename], '-mat');
catch
    return;
end

g = finputcheck(varargin, {'comps'          'integer'   []             []; ...
                           'factors'        'integer'   []             []; ...
                           'frqlim'         'real'      []             []; ...
                           'smoothing'      'real'   []             [40];...
                           'plotICtf'       'string'    {'on' 'off'}   'on';... 
                           'plotPCtf'       'string'    {'on' 'off'}   'on';...
                           'plotIMtf'       'string'    {'on' 'off'}   'on';...
                           'plotIMtime'     'string'    {'on' 'off'}   'on';...
                            }, 'inputgui');
if isstr(g), error(g); end;

if nargin ==1
    plotTypes2funtc = {'plotICtf', 'plotPCtf', 'plotIMtfims', 'plotIMtime'};
     plotTypes = {'IC spectrogram', 'Summed IM backprojection', 'Combined IC-IM spectrogram', 'IM timecourse'};
    freqLim = tmpIMA.IMA.freqlim;
    ic_list = sprintfc('%d',tmpIMA.IMA.complist);
    im_list = sprintfc('%d',[1:tmpIMA.IMA.npcs]);
    plotopt = repmat({'off'}, 1,4);
    
    uilist = {{'style' 'text' 'string' 'Plot type'} {'style' 'popupmenu'  'string' plotTypes 'tag' 'plottype' 'value' 1}...
              {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit' 'string' num2str(freqLim) 'tag' 'freqlimits'}...
              {'style' 'text' 'string' 'Select ICs', }     {'style' 'text' 'string' 'Select IMs'} ...
              {'style'  'list'  'string' ic_list 'max',200,'min',1,'Tag','icindx'} {'style'  'list'  'string' im_list 'max',200,'min',1,'Tag','imindx' } {}};
    
    ht = 6; wt = 2 ;
    geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [1 1]}...
            {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}...
            {wt ht [0 2]  [1 1]} {wt ht [1 2]  [1 1]}...
            {wt ht [0 3]  [1 3]} {wt ht [1 3]  [1 3]}...
            {wt ht [0 5]  [1 3]}};
   
    [result, ~, ~, resstruct, ~] = inputgui('title','Plot IM decomposition -- pop_plotspecdecomp', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_plotspecdecomp'');');
    if isempty(result), return; end;
    g.comps = tmpIMA.IMA.complist(resstruct.icindx);
    g.factors = resstruct.imindx;
    g.frqlim = str2num(resstruct.freqlimits);
    plotopt{resstruct.plottype} = 'on';
end

if isempty(g.frqlim)
    g.frqlim = tmpIMA.IMA.freqlim;
end

if isempty(g.comps)
    g.comps = tmpIMA.IMA.complist;
end

if isempty(g.factors)
    g.factors = 1:tmpIMA.IMA.npcs;
end

plotIMtimecourse(tmpIMA.IMA, EEG,'comps',g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
                 'smoothing', g.smoothing,...
                 'plotICtf',plotopt{1}, 'plotPCtf', plotopt{2},...
                 'plotIMtf', plotopt{3}, 'plotIMtime', plotopt{4});