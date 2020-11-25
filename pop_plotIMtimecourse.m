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
% smoothing - [number between 0 and 1] factor to smooth IM timecourses. A
%              number more close to 1 means less smoothing
% plotICtf  -  [string]  {'on' 'off'} plot IC time-frequency decomposition, default 'off'                    
% plotPCtf  -   [string] {'on' 'off'}  plot PC time-frequency backprojection, default 'off'                      
% plotIMtf  -    [string] {'on' 'off'}  plot IM time-frequency backprojection, default 'off'                      
% plotIMtime  -  [string] {'on' 'off'}   plot IM weights timecourse, default 'off'                     


function pop_plotIMtimecourse(EEG,varargin)

g = finputcheck(varargin, {'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'smoothing'      'integer'   []             [1];...
    'plotICtf'       'string'    {'on' 'off'}                       'on';...
    'plotPCtf'       'string'    {'on' 'off'}                       'on';...
    'plotIMtf'       'string'    {'on' 'off'}                       'on';...
    'plotIMtime'     'string'    {'on' 'off'}                       'on';...
    }, 'inputgui');
if isstr(g), error(g); end;


load([EEG.etc.IMA.filepath '/' EEG.etc.IMA.filename], '-mat');



if isempty(g.frqlim)
    g.frqlim = IMA.freqlim;
end

if isempty(g.comps)
    g.comps = IMA.complist;
end

if isempty(g.factors)
    g.factors = 1:IMA.npcs;
end

plotIMtimecourse(IMA,'comps',g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
                 'smoothing', g.smoothing,...
                 'plotICtf', g.plotICtf, 'plotPCtf', g.plotPCtf,...
                 'plotIMtf', g.plotIMtf, 'plotIMtime', g.plotIMtime)



