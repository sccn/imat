%
% Plots IMA results for a single subject in a study: plots IMA templates in a matrix for all IMs and
% ICs separately, or all ICs for each IM superimposed on a single axis, or
% each IC scalp map with super imposed IM templates on a single axis
%
%
%
% pop_plotspecdecomp_study(STUDY,varargin)
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
% Example: plot all ICs and IMs separately for IMA decomposition of subject 3
% >> pop_plotspecdecomp_study(STUDY, 'subject', 3, 'plottype', 'comb')
%
%
% Example: plot all ICs with a selection of superimposed IMs for IMA decomposition of subject 3
% >> pop_plotspecdecomp_study(STUDY, 'subject', 3, 'plottype', 'ims', 'factors', [1 3 4 6 7], 'maps', 'on')
%
%
%
% INPUTS
% STUDY - STUDY structure pointing to precomputed IMA information and EEG
%         datasets
% subject - subject to plot - provide subject code/number as string i.e. 'S3'
%           if empty plots results for all subjects in the study
% comps - independent components to plot
% factors - IMs to plot
% frqlim - frequency limits for plotting
% plottype -- ['ims', 'ics' 'comb'] 'ics' will plot a single axis with all comps
%             superimposed.
%             'ims' will show each IC scalp map with super imposed IM templates.
%             'comb' will plot all ICs and IMs plotted separately.
% maps -- ['on','off'] if 'on', then will plot scalp maps.

function pop_plotspecdecomp_study(STUDY,varargin)

% checking inputs
% ---------------
g = finputcheck(varargin, {'subject'     'string'   {}             ''; ...
    'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'freqscale'     'string'     {'log' 'linear'}         'log';...
    'plottype'      'string'     {'ics' 'ims' 'comb'}           'comb'; ...
    'maps'          'string'     {'on', 'off'}  'on';...
    }, 'inputgui');
if isstr(g), error(g); end;


%% getting activations
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
    
    
    plotspecdecomp(IMA, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
        'freqscale', g.freqscale, 'plottype', g.plottype, 'maps', g.maps);
    
end

