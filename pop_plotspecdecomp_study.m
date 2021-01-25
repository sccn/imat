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

function pop_plotspecdecomp_study(STUDY, varargin)

% checking inputs
% ---------------
g = finputcheck(varargin, { 'subject'       'string'     {}             ''; ...
                            'comps'         'integer'    []             []; ...
                            'factors'       'integer'    []             []; ...
                            'frqlim'        'real'       []             []; ...
                            'freqscale'     'string'     {'log' 'linear'}         'log';...
                            'plottype'      'string'     {'ics' 'ims' 'comb'}           'comb'; ...
                            'maps'          'string'     {'on', 'off'}  'on';...
                           }, 'inputgui');
                       
if isstr(g), error(g); end;

if ~isfield(STUDY.etc, 'IMA'), return; end
try
    tmpIMA = load([STUDY.etc.IMA.imafilepath{1} filesep STUDY.etc.IMA.imafilename{1}], '-mat');
catch
    return;
end

if nargin ==1
    plotTypes2funtc = {'comb', 'ics', 'ims'};
    plotTypes = {'IM mode decomposition', 'Superimposed IC modes', 'Superimposed IM modes'};
    subj_list = STUDY.etc.IMA.subject;
    
    % To be updated in GUI
    freqLim = tmpIMA.IMA.freqlim;
    ic_list = sprintfc('%d',tmpIMA.IMA.complist);
    im_list = sprintfc('%d',[1:tmpIMA.IMA.npcs]);
    
    subjcallback = ['subjindx = get(findobj(''Tag'', ''subjname''),''value'');'...
                   'tmpIMA = load([STUDY.etc.IMA.imafilepath{subjindx} filesep STUDY.etc.IMA.imafilename{subjindx}], ''-mat'');'...
                   'ic_list = sprintfc(''%d'',tmpIMA.IMA.complist);'...
                   'im_list = sprintfc(''%d'',[1:tmpIMA.IMA.npcs]);'...
                   'set(findobj(''Tag'', ''icindx''), ''string'', ic_list);'...
                   'set(findobj(''Tag'', ''imindx''), ''string'', im_list);'];
    
    uilist = {{'style' 'text' 'string' 'Subject'} {'style' 'popupmenu'  'string' subj_list 'tag' 'subjname' 'value' 1 'callback' subjcallback}...
              {'style' 'text' 'string' 'Plot type'} {'style' 'popupmenu'  'string' plotTypes 'tag' 'plottype' 'value' 1}...
              {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit' 'string' num2str(freqLim) 'tag' 'freqlimits'}...
              {'style' 'text' 'string' 'Select ICs', }     {'style' 'text' 'string' 'Select IMs'} ...
              {'style'  'list'  'string' ic_list 'max',200,'min',1,'Tag','icindx'} {'style'  'list'  'string' im_list 'max',200,'min',1,'Tag','imindx' } {}};
    
    ht = 6; wt = 2 ;
    geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [1 1]}...
            {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}...
            {wt ht [0 2]  [1 1]} {wt ht [1 2]  [1 1]}...
            {wt ht [0 3]  [1 1]} {wt ht [1 3]  [1 1]}...
            {wt ht [0 4]  [1 3]} {wt ht [1 4]  [1 3]}...
            {wt ht [0 5]  [1 3]}};
   
    [result, ~, ~, resstruct, ~] = inputgui('title','Plot IM decomposition -- pop_plotspecdecomp', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_plotspecdecomp'');');
    
    
    if isempty(result), return; end;
    g.subject = subj_list{resstruct.subjname};
    g.plottype = plotTypes2funtc{resstruct.plottype};
    g.comps = tmpIMA.IMA.complist(resstruct.icindx);
    g.factors = resstruct.imindx;
    
    tmpIMA = load([STUDY.etc.IMA.subjfilepath{1}{1} filesep STUDY.etc.IMA.subjfilename{1}{1}(1:end-4) '.ima'], '-mat'); % Fox to load correct subject file given info from GUI
    g.frqlim = str2num(resstruct.freqlimits);
end

%% getting activations
if isempty(g.subject)
    subjcode = STUDY.subject;
else
    subjcode = {g.subject};
end


for iko = 1:length(subjcode)
    indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));
    
    %% load IMA file
   % load([STUDY.datasetinfo(indsj(1)).filepath filesep STUDY.etc.IMA.imafilename{iko,:}], '-mat' );
   load([STUDY.etc.IMA.imafilepath{iko} filesep STUDY.etc.IMA.imafilename{iko}], '-mat' );
   
     
    EEG = pop_loadset('filename',IMA.subjfilename{1},'filepath',IMA.subjfilepath{1});

    if isempty(g.frqlim)
        g.frqlim = IMA.freqlim;
    end
    
    if isempty(g.comps)
        g.comps = IMA.complist;
    end
    
    if isempty(g.factors)
        g.factors = 1:IMA.npcs;
    end
    
    
    plotspecdecomp(IMA, EEG, 'comps', g.comps, 'factors', g.factors, 'frqlim', g.frqlim,...
        'freqscale', g.freqscale, 'plottype', g.plottype, 'maps', g.maps);
    
end

