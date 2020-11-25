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

% checking inputs
% ---------------
g = finputcheck(varargin, { 'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'freqscale'     'string'     {'log' 'linear'}         'log';...
    'plottype'      'string'     {'ics' 'ims' 'comb'}           'comb'; ...
    'maps'          'string'     {'on', 'off'}  'on';...
    }, 'inputgui');
if isstr(g), error(g); end;


load([EEG.etc.IMA.filepath '/' EEG.etc.IMA.filename], '-mat');

if isempty(g.frqlim)
    g.frqlim = IMA.freqlim;
end

if isempty(g.comps)
    g.plotcomps = IMA.complist;
end

if isempty(g.factors)
    g.factors = 1:IMA.npcs;
end


plotspecdecomp(IMA, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
    'freqscale', g.freqscale, 'plottype', g.plottype, 'maps', g.maps);


