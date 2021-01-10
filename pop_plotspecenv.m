%
% plot envelopes of IC spectra, PC backprojection and IM templates on the same axis 
%
% pop_plotspecenv(EEG,varargin);
% 
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
% 
% Example plotting the spectra for a subject for ICs 1 3 5 6 and IMs 1 2 3
% 5 7 with frequency limits 3 to 40
% 
% >> pop_plotspecenv(EEG, 'comps', [1 3 5 6], 'factors', [1 2 3 5 7], 'frqlim', [3 40]);
%
%
% Example plotting the spectra for a subject for IC 1 and IMs 1 and 2 
%  with frequency limits 3 to 60, plotting the upper 95th percentile 
%
% >> pop_plotspecenv(EEG, 'comps', [1], 'factors', [1 2], 'frqlim', [3 60], 'plotenv', 'upper');%, 'setminmax', [30 80]);
%
%
%
% INPUTS:
% required Inputs:
% EEG - EEG structure with information on stored IMA file
%
% optional Inputs:
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject 
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject% frqlim - [min max] frequency limits for plotting
% plotenv - [string] possible values: 'env', 'upper', 'lower'. Either plot envelope of IMs, plot upper 95th percentile, 
%           plot lower 95th percentile. (default 'env')
% plotperc - [number between 0 and 1] define a percentile of trials to use to exclude outliers...
%            (default 95th percentile - 0.95)
% plotcond - off or on (default off). If IMA over different conditions was computed allows to plot conditions separately ;
%            the function is using the condition specific IC mean spectrum and the condition specific backprojection of IMs
%            the condition specific IC mean spectrum is added to the IC spectral envelope, to the PC spectral envelope and to the IM
%            backprojections
% setmimax - [min max] set minimum / maximum of spectrum yscale

function [IMA] = pop_plotspecenv(EEG,varargin);
IMA = [];

if ~isfield(EEG.etc, 'IMA'), return; end
try
    tmpIMA = load([EEG.etc.IMA.filepath '/' EEG.etc.IMA.filename], '-mat');
catch
    return;
end

g = finputcheck(varargin, { 'comps'        'integer'    []                           []; ...
                            'factors'      'integer'    []                           []; ...
                            'frqlim'       'real'       []                           []; ...
                            'plotenv'      'string'     { 'env' 'upper' 'lower'}     'env'; ...
                            'plotperc'     'real'       []                           [0.95]; ...
                            'setminmax'    'real'       []                           [];...
                           }, 'inputgui');
if isstr(g), error(g); end;

if nargin ==1
    EnvtTypes2funtc = {'upper' 'lower' 'env'};
    envTypes = {'Upper envelope      ', 'Lower envelope     ', 'Full Envelope                 '};
               
    freqLim = tmpIMA.IMA.freqlim;
    ic_list = sprintfc('%d',tmpIMA.IMA.complist);
    im_list = sprintfc('%d',[1:tmpIMA.IMA.npcs]);
    
    uilist = {{'style' 'text' 'string' 'Envelope type'} {'style' 'popupmenu'  'string' envTypes 'tag' 'plottype' 'value' 3}...
              {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit' 'string' num2str(freqLim) 'tag' 'freqlimits'}...
              {'style' 'text' 'string' 'Select ICs', }     {'style' 'text' 'string' 'Select IMs'} ...
              {'style'  'list'  'string' ic_list 'max',200,'min',1,'Tag','icindx'} {'style'  'list'  'string' im_list 'max',200,'min',1,'Tag','imindx' } {}};
    
    ht = 6; wt = 2 ;
    geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [1 1]}...
            {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}...
            {wt ht [0 2]  [1 1]} {wt ht [1 2]  [1 1]}...
            {wt ht [0 3]  [1 3]} {wt ht [1 3]  [1 3]}...
            {wt ht [0 5]  [1 3]}};
   
    [result, ~, ~, resstruct, ~] = inputgui('title','Plot spectral envelope -- pop_plotspecenv', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_plotspecenv'');');
    if isempty(result), return; end;
    g.comps = tmpIMA.IMA.complist(resstruct.icindx);
    g.factors = resstruct.imindx;
    g.plotenv = EnvtTypes2funtc{resstruct.plottype};
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


fprintf('\n Plotting overall data and meanspec...');
dataport = tmpIMA.IMA.timepntCond;
meanspec = tmpIMA.IMA.meanpwr;
plotspecenv(tmpIMA.IMA, EEG, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim, 'plotenv', g.plotenv, 'plotperc', g.plotperc, 'setminmax', g.setminmax, 'dataport', dataport, 'meanspec', meanspec);