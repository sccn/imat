%
% plot envelopes of IC spectra, PC backprojection and IM templates on the same axis
%
% pop_plotspecenv_study(STUDY,varargin);
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
% Example plotting the spectra for subject 3 for ICs 1 3 5 6 and IMs 1 2 3
% 5 7 with frequency limits 3 to 40, plotting the condtions separately
%
% >> pop_plotspecenv_study(STUDY, 'subject', 'S3', 'comps', [1 3 5 6], 'factors', [1 2 3 5 7], 'frqlim', [3 40],'plotcond', 'on');
%
%
% Example plotting the spectra for subject 5 for IC 1 and IMs 1 and 2
%  with frequency limits 3 to 60, plotting the upper 95th percentile of the overall spectra over conditions
%
% >> pop_plotspecenv_study(STUDY, 'subject', 'S5', 'comps', [1], 'factors', [1 2], 'frqlim', [3 60], 'plotenv', 'upper');%, 'setminmax', [30 80]);
%
%
%
% INPUTS:
% required Inputs:
% STUDY - STUDY structure pointing to precomputed IMA information and EEG
%         datasets
%
% optional Inputs:
% subject - subject to plot - provide subject code/number as string i.e. 'S3'
%           if empty plots results for all subjects in the study
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject
% frqlim - [min max] frequency limits for plotting
% plotenv - [string] possible values: 'env', 'upper', 'lower'. Either plot envelope of IMs (upper and lower 95th percentile),
%            or plot upper 95th percentile, or plot lower 95th percentile. (default 'env')
% plotperc - [number between 0 and 1] define a percentile of trials to use to exclude outliers...
%            (default 95th percentile - 0.95)
% plotcond - [string] 'off' or 'on' (default off). If IMA over different conditions was computed allows to plot conditions separately;
%            the function is using the condition specific IC mean spectrum and the condition specific backprojection of IMs
%            the condition specific IC mean spectrum is added to the IC spectral envelope, to the PC spectral envelope and to the IM
%            backprojections
% setmimax - [min max] set minimum / maximum of spectrum yscale

function [IMA] = pop_plotspecenv_study(STUDY,varargin);


g = finputcheck(varargin, {'subject'     'string'   {}             ''; ...
    'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'plotenv'      'string'     { 'env' 'upper' 'lower'}           'env'; ...
    'plotperc'      'real'       []           [0.95]; ...
    'plotcond'       'string'    {'on' 'off'}                       'on';...
    'setminmax'        'real'      []             [];...
    }, 'inputgui');
if isstr(g), error(g); end;

if isempty(g.subject)
    subjcode = STUDY.subject;
else
    subjcode = {g.subject};
end

for iko = 1:length(subjcode)
    indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));
    
    %% load IMA file
    load([STUDY.datasetinfo(indsj(1)).filepath '/' STUDY.etc.IMA.imafilename(iko,:)], '-mat' );
    
    if isempty(g.frqlim)
        g.frqlim = IMA.freqlim;
    end
    
    if isempty(g.comps)
        g.comps = IMA.complist;
    end
    
    if isempty(g.factors)
        g.factors = 1:IMA.npcs;
    end
    
    
    if strcmp(g.plotcond, 'on') && length(IMA.ntrials) > 1
        fprintf('\n Plotting conditionwise data and meanspec...');
        for ib = 1:length(IMA.ntrials)
            dataport = IMA.timepntCond{ib}; % select trials for condition
            meanspec = squeeze(IMA.meanpwrCond(ib,:,:));
            plotspecenv(IMA, 'comps' ,g.comps, 'factors',g.factors, 'frqlim',g.frqlim, 'plotenv',g.plotenv, 'plotperc', g.plotperc, 'setminmax' ,g.setminmax, 'dataport', dataport, 'meanspec' ,meanspec);
            ph1 = gcf;
            ph1 = textsc([IMA.subj{1}  ' ' IMA.condition{ib}],'title');
            set(ph1,'fontsize',20);
        end
    else
        fprintf('\n Plotting overall data and meanspec...');
        dataport = 1:(IMA.ntrials*IMA.ntw_trials);
        meanspec = IMA.meanpwr;
        plotspecenv(IMA, 'comps' ,g.comps, 'factors',g.factors, 'frqlim',g.frqlim, 'plotenv',g.plotenv, 'plotperc', g.plotperc, 'setminmax' ,g.setminmax, 'dataport', dataport, 'meanspec' ,meanspec);
        ph1 = gcf;
        ph1 = textsc([IMA.subj{1}],'title');
        set(ph1,'fontsize',20);
    end
    
end
%