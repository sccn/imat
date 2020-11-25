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


g = finputcheck(varargin, { 'comps'     'integer'   []             []; ...
                           'factors'       'integer'    []             []; ...
                           'frqlim'        'real'       []             []; ...
                           'plotenv'      'string'     { 'env' 'upper' 'lower'}           'env'; ...
                           'plotperc'      'real'       []           [0.95]; ...
                           'setminmax'        'real'      []             [];...
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


    fprintf('\n Plotting overall data and meanspec...');
    dataport = IMA.timepntCond;
    meanspec = IMA.meanpwr;
    plotspecenv(IMA, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim, 'plotenv', g.plotenv, 'plotperc', g.plotperc, 'setminmax', g.setminmax, 'dataport', dataport, 'meanspec', meanspec);



   
%