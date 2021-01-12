%
% Plots IMA results: plots IMA templates in a matrix for all IMs and
% ICs separately, or all ICs for each IM superimposed on a single axis, or
% each IC scalp map with super imposed IM templates on a single axis
%
% pop_plotspecdecomp(EEG,varargin)
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
%
% Example: plot all ICs and IMs separately for IMA decomposition of a subject
% >> pop_plotspecdecomp(EEG, 'plottype', 'comb')
%
%
% Example: plot all ICs with a selection of superimposed IMs for IMA decomposition of a subject
% >> pop_plotspecdecomp(EEG, 'plottype', 'ims', 'factors', [1 3 4 6 7], 'maps', 'on')
%
%
% INPUTS
% EEG - EEG structure with precomputed IMA information included
% plotcomps - independent components to plot
% factors - IMs to plot
% frqlim - frequency limits for plotting
% plottype -- ['ims', 'ics' 'comb'] 'ics' will plot a single axis with all comps
%             superimposed.
%             'ims' will show each IC scalp map with super imposed IM templates.
%             'comb' will plot all ICs and IMs plotted separately.
% maps -- ['on','off'] if 'on', then will plot scalp maps.

function [IMA] = pop_plotspecdecomp(EEG,varargin);
IMA = [];

if ~isfield(EEG.etc, 'IMA'), return; end
try
    tmpIMA = load([EEG.etc.IMA.filepath '/' EEG.etc.IMA.filename], '-mat');
catch
    return;
end

% checking inputs
% ---------------
g = finputcheck(varargin, { 'comps'         'integer'    []                         []; ...
                            'factors'       'integer'    []                         []; ...
                            'frqlim'        'real'       []                         []; ...
                            'freqscale'     'string'     {'log' 'linear'}           'log';...
                            'plottype'      'string'     {'ics' 'ims' 'comb'}       'comb'; ...
                            'maps'          'string'     {'on', 'off'}              'on';...
                            }, 'inputgui');
if isstr(g), error(g); end;

if nargin ==1
    plotTypes2funtc = {'comb', 'ics', 'ims'};
     plotTypes = {'IM mode decomposition', 'Superimposed IC modes', 'Superimposed IM modes'};
    freqLim = tmpIMA.IMA.freqlim;
    ic_list = sprintfc('%d',tmpIMA.IMA.complist);
    im_list = sprintfc('%d',[1:tmpIMA.IMA.npcs]);
    
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
    g.plottype = plotTypes2funtc{resstruct.plottype};
    g.frqlim = str2num(resstruct.freqlimits);
    
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


plotspecdecomp(tmpIMA.IMA, EEG, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
    'freqscale', g.freqscale, 'plottype', g.plottype, 'maps', g.maps);



