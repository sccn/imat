% plot original IC time-frequency decomposition, PC time-frequency backprojection, IM time-frequency backprojection and IM
% weights
% 
% 
% pop_plotIMtimecourse_study(STUDY,varargin)
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from functions written by Julie Onton
%
%
%
% Example: for subject 4 and all conditions plot components 1, 5 and 8 with frequency limit 6-40Hz and smoothing 0.1
%           - plot only IC time-frequency decomposition and PC time-frequency backprojections 
%
% >>  pop_plotIMtimecourse_study(STUDY, 'subject', 'S4','comps', [1 5 8],...
%      'frqlim', [6 40], 'plotcond', 'on', 'plotICtf', 'on', 'plotPCtf', 'on')
%
%
% Example: for subject 5 and all conditions plot components 3, 5 and 6 and IMs 2, 6 and 9 with frequency limit 6-40Hz and smoothing 0.1
%          plot only IM time-frequency backprojections on ICs and IM weight
%           timecourses
%
% >>  pop_plotIMtimecourse_study(STUDY, 'subject', 'S5','comps', [3 5 6], 'factors', [2 6 9],...
%      'frqlim', [6 40], 'plotcond', 'on','smoothing', 0.1, 'plotIMtf', 'on', 'plotIMtime', 'on')
%
% INPUTS:
% STUDY - STUDY structure pointing to precomputed IMA information and EEG
%         datasets
% subject - subject to plot - provide subject code/number as string i.e. 'S3'
%           if empty plots results for all subjects in the study
% comps - [vector] independent components to plot - if empty plots all
%                      the independent components for a subject 
% factors - [vector] IMs to plot - if empty plots all the IMs for a subject% frqlim - [min max] frequency limits for plotting
% frqlim - [min max] frequency limits for plotting
% plotcond - [string] 'off' or 'on' (default 'off'). If IMA over different
%                      conditions was computed plots a vertical line in between conditions and
%                      plots conditions in different colors for IM weights plotting
% smoothing - [number between 0 and 1] factor to smooth IM timecourses. A
%              number more close to 1 means less smoothing
% plotICtf  -  [string]  {'on' 'off'} plot IC time-frequency decomposition, default 'off'                    
% plotPCtf  -   [string] {'on' 'off'}  plot PC time-frequency backprojection, default 'off'                      
% plotIMtf  -    [string] {'on' 'off'}  plot IM time-frequency backprojection, default 'off'                      
% plotIMtime  -  [string] {'on' 'off'}   plot IM weights timecourse, default 'off'                     


function pop_plotIMtimecourse_study(STUDY,varargin)

g = finputcheck(varargin, { 'subject'     'string'   {}             ''; ...
    'comps'     'integer'   []             []; ...
    'factors'       'integer'    []             []; ...
    'frqlim'        'real'       []             []; ...
    'plotcond'       'string'    {'on' 'off'}                       'off';...
    'smoothing'      'integer'   []             [1];...
    'plotICtf'       'string'    {'on' 'off'}                       'off';...
    'plotPCtf'       'string'    {'on' 'off'}                       'off';...
    'plotIMtf'       'string'    {'on' 'off'}                       'off';...
    'plotIMtime'     'string'    {'on' 'off'}                       'off';...
    }, 'inputgui');
if isstr(g), error(g); end;


%% check which subject/s to plot
if isempty(g.subject)
    subjcode = STUDY.subject;
else
    subjcode = {g.subject};
end


%% run over specified subjects
for iko = 1:length(subjcode)
indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));

%% load IMA file
load([STUDY.datasetinfo(indsj(1)).filepath '/' STUDY.etc.IMA.imafilename(iko,:)], '-mat' );


%% check if frequency limits, comps and IMs are empty 
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
                 'plotcond', g.plotcond, 'smoothing', g.smoothing,...
                 'plotICtf', g.plotICtf, 'plotPCtf', g.plotPCtf,...
                 'plotIMtf', g.plotIMtf, 'plotIMtime', g.plotIMtime)

end


