%
% plots spectral templates, dipole density and scalpmaps for IM clusters
%
%
% pop_plotIMAcluster(STUDY,varargin);
%
% 
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC
% San Diego, Feb 2020
% adapted from a function written by Julie Onton
%
% Example plotting cluster numbers 1 and 2 and restricting spectra to 3 60
% Hz and log frequency scale 
%  >>  pop_plotIMAcluster(STUDY, 'clust', [1 2], 'freqlim', [3 60], 'freqscale', 'log')
%
% Example plotting subclusters and scalpmaps
%  >>  pop_plotIMAcluster(STUDY, 'plotsubclusters', 'on', 'plotscalpmaps', 'on')
%
%
% INPUTS:
% required Inputs:
% STUDY - STUDY structure with information on stored IMA files
%
% optional Inputs:
% plotclust - 'string' 'on' or 'off' - default 'on' - may be turned off
%               when wanting to plot only subclusters
% clust - [vector] 'vector of integers' vector of cluster indices to plot
%                  if empty plot all clusters 
% freqscale -  'string' frequency scale to plot spectra 'linear' or 'log'
%                       default 'log'
% freqlim - [min max] 'integers'  provide min - max frequency range; ...
%
% plotsubclusters - 'string' 'on' or 'off' specifies whether to plot
%                    subclusters - plots only the subclusters previously
%                    defined in pop_subclusterIMAtemplates - default 'off'
% plottemplates - 'string' 'on' or 'off' specifies whether to plot
%                  spectral templates - default 'off'
% plotdipsources - 'string' 'on' or 'off' specifies whether to plot
%                  dipole density plots - default 'off'
% plotscalpmaps - 'string' 'on' or 'off' specifies whether to plot
%                  scalpmaps - default 'off'
%
% OUTPUTS: plots of diploedensity, spectra and scalpmaps of clusters and subclusters 



function pop_plotIMAcluster(STUDY,varargin);

g = finputcheck(varargin, {'freqlim'        'integer'       []             []; ...
    'clust'      'integer'     []         [];...
    'plotclust'      'string'   {'on' 'off'}       'on';...
    'freqscale'        'string'    {'linear' 'log'}       'log';...
    'plotsubclusters'      'string'     {'on' 'off'}         'off';...
    'plottemplates'     'string'    {'on' 'off'}      'off';...
    'plotdipsources'     'string'    {'on' 'off'}      'off';...
    'plotscalpmaps'     'string'    {'on' 'off'}      'off';...
    }, 'inputgui');
if isstr(g), error(g); end;

if nargin == 1
    onofflabel = {'off' 'on'};
    freqscalelabel = {'linear' 'log'};
    uilist = {{'style' 'text' 'string' 'IM clusters'} {'style' 'edit'  'string' '' 'tag' 'clustnum'}...
        {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit' 'string' '' 'tag' 'freqlimits'}...
        {'style' 'text' 'string' 'Freq. scale'} {'style' 'popupmenu' 'string' {'linear', 'log'} 'tag' 'frescale' 'value' 2}...
        {'style' 'text' 'string' 'Plotting options'} ...
        {'style' 'checkbox' 'string' 'Templates' 'Tag' 'cbox_plottemplates'} ...
        {'style' 'checkbox' 'string' 'Dipoles' 'Tag' 'cbox_plotdipoles'} ...
        {'style' 'checkbox' 'string' 'Scalp maps' 'Tag' 'cbox_plotscalpmaps'} ...
        };
    
%             {'style' 'checkbox' 'string' 'Clusters' 'Tag' 'cbox_plotcluster'} ...
%         {'style' 'checkbox' 'string' 'Subclusters' 'Tag' 'cbox_pltsbcluster'} ...
%     
    ht = 5; wt = 4 ;
    geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [3 1]}...
            {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}  {wt ht [2.2 1]  [1 1]} {wt ht [3 1]  [1 1]}...
            {wt ht [0 2]  [1 1]}...
            {wt ht [1 3]  [1 1]}...
            {wt ht [1 4]  [1 1]}...
            {wt ht [1 5]  [1 1]}...
            };
    
    [result, ~, ~, resstruct, ~] = inputgui('title','Plot IM clusters -- pop_plotIMAcluster', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_plotIMAcluster'');');
    if isempty(result), return; end;
    g.freqlim = str2num(resstruct.freqlimits);
    g.clust = str2num(resstruct.clustnum);
    g.freqscale = freqscalelabel{resstruct.frescale};
    g.plottemplates = onofflabel{resstruct.cbox_plottemplates+1};
    g.plotdipsources = onofflabel{resstruct.cbox_plotdipoles+1};
    g.plotscalpmaps = onofflabel{resstruct.cbox_plotscalpmaps+1};
end

subjcode = STUDY.subject;

templates = [];
IMICindex = [];
dipsources = [];
scalpmaps = [];

for iko = 1:length(subjcode)
    indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));
    
    %% load IMA file for curent subject
    load([STUDY.etc.IMA.imafilepath{iko} filesep STUDY.etc.IMA.imafilename{iko}], '-mat' );
  
    str = string(STUDY(iko).subject);
    sujnum = sscanf(str,'S%d');
    
    if isempty(g.freqlim)
        g.freqlim = IMA.freqlim;
    end
    
    [val, freqind1] = min(abs(IMA.freqvec-g.freqlim(1)));
    [val, freqind2] = min(abs(IMA.freqvec-g.freqlim(2)));
    
    templates = [templates IMA.precluster.templates(:,freqind1:freqind2)];
    dipsources = [dipsources IMA.precluster.dipsources];
    scalpmaps = [scalpmaps IMA.precluster.scalpmaps];
    freqvec = IMA.freqvec(freqind1:freqind2);
end

plotIMAcluster(STUDY.etc.IMA.clustidx, 'plotclust', g.plotclust, 'clust', g.clust,...
    'freqscale', IMA.freqscale,'plottemplates', g.plottemplates, 'templates', templates,...
    'freqvec', freqvec, 'plotdipsources', g.plotdipsources, 'dipsources', dipsources,...
    'plotscalpmaps', g.plotscalpmaps, 'scalpmaps', scalpmaps,'plotsubclusters', g.plotsubclusters)
