%  Takes WT comodulation spectral templates as collected in IMA.precluster.templates 
%  and clusters to specified number of clusters, uses k-means clustering
%
% [STUDY] = pop_clusterIMAtemplates(STUDY,ALLEEG, varargin);
%
% 
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
% Example:
% >> [STUDY] = pop_clusterIMAtemplates(STUDY, ALLEEG, 'freqlim', [3 120] ,'nclust', 3, 'pcs', 20);
%
% Example:
% >> [STUDY] = pop_clusterIMAtemplates(STUDY, ALLEEG, 'freqlim', [3 30] ,'nclust', 6);
%
%
% INPUTS:
% STUDY - STUDY structure with information on stored IMA files
% ALLEEG - ALLEEG structure
% freqlim - 'integer'   [min max] frequency limits to restrict clustering
%            to i.e. alpha frequency range [8 14] default whole frequency
%            range in IMA.freqlim
% nclust -- [integer] if not empty, will divide templates into nclust clusters
%                     by 'pdist','linkage','dendrogram' Matlab functions
% pcs -- number of principal IM spectral template dimensions to retain
% method -- ['k-means']  k-means clustering
% dipole_locs -- 'string' {'on' 'off'}  use dipole locations in addition to spectral templates for clustering, default 'off' 
% weightSP -- 'integer' A number between 1 and 20. Weighting to use for spectral templates when clustering on spectral templates and dipole locations default is 1.
%               A larger number will give more weight to spectral templates compared to dipole locations
% weightDP -- 'integer' A number between 1 and 20. Weighting to use for dipoles when clustering on spectral templates and dipole locations default is 1. a larger 
%               number will give more weight to dipole locations compared to spectral templates
 
%
% OUTPUT
% STUDY - STUDY structure with indexes of clusters (column 4 in matrix) and distance to cluster
% centroid saved in STUDY.etc.IMA.clustidx in the form [SJ IM ICs clusindx] and STUDY.etc.IMA.distance
%

function [STUDY] = pop_clusterIMAtemplates(STUDY, ALLEEG, varargin);

%clustmethods = {'corr' 'euc' 'kmeans'};
g = finputcheck(varargin, {'freqlim'        'integer'       []             []; ...
    'nclust'     'integer'     []         [];...
    'pcs'     'integer'     []         [10];...
    'method'      'string'       {'kmeans'}             'kmeans'; ...  
    'dipole_locs' 'string' {'on' 'off'}  'off';...
    'weightSP' 'integer' [] [1];...
    'weightDP' 'integer' [] [1];...
    }, 'inputgui');
if isstr(g), error(g); end;

if nargin == 2
    opt_offon = {'off', 'on'};
    clsutmethodslist = {'Kmeans'}; % for display in GUI only 
    clustmethods     = {'kmeans'};
    
     cbk_dipole = ['valchbx = get(findobj(''Tag'', ''chbx_enabledipole''), ''Value'');' ...
                  'if valchbx,'...'
                  'set(findobj(''Tag'', ''ed_weightSP''), ''enable'', ''on'');'...
                  'set(findobj(''Tag'', ''ed_weightDP''), ''enable'', ''on'');'...
                  'else,'...
                  'set(findobj(''Tag'', ''ed_weightSP''), ''enable'', ''off'');'...
                  'set(findobj(''Tag'', ''ed_weightDP''), ''enable'', ''off'');'...   
                  'end;'];
              
      uilist = {{'style' 'text' 'string' 'Method'} {'style' 'popupmenu'  'string' clsutmethodslist 'tag' 'method' 'value' 1 }...
              {'style' 'text' 'string' 'Number of clusters'} {'style' 'edit'  'string' ' ' 'tag' 'nclust'}...
              {'style' 'text' 'string' 'Number of PCs'} {'style' 'edit'  'string' ' ' 'tag' 'npcs'} {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit'  'string' ' ' 'tag' 'freqlim'}...
              { 'style' 'Checkbox'   'string' 'Complement clustering with dipole location' 'Tag', 'chbx_enabledipole', 'Value' 0,  'enable', 'on' 'callback' cbk_dipole} ...
              {'style' 'text' 'string' 'Template weight'} {'style' 'edit'  'string' num2str(g.weightSP) 'tag' 'ed_weightSP'  'enable', 'off'} {'style' 'text' 'string' 'Dipole weight'} {'style' 'edit'  'string' num2str(g.weightDP) 'tag' 'ed_weightDP'  'enable', 'off'}};
     
     ht = 4; wt = 2 ;
     geom = {{wt ht [0 0]  [1 1]} {wt ht [0.6 0]  [1 1]}...
             {wt ht [0 1]  [1 1]} {wt ht [0.6 1]  [0.4 1]}...
             {wt ht [0 2]  [1 1]} {wt ht [0.6 2]  [0.4 1]} {wt ht [1.2 2]  [1 1]} {wt ht [1.6 2]  [0.4 1]}...
             {wt ht [0 3]  [1 1]}...
             {wt ht [0 4]  [1 1]} {wt ht [0.6 4]  [0.4 1]} {wt ht [1.2 4]  [1 1]} {wt ht [1.6 4]  [0.4 1]}};
     
    [result, ~, ~, resstruct, ~] = inputgui('title','Cluster IM templates -- pop_clusterIMAtemplates', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_clusterIMAtemplates'');');
    
    if isempty(result), return; end;
    g.freqlim = str2num(resstruct.freqlim);
    g.nclust = str2num(resstruct.nclust);
    g.pcs = str2num(resstruct.npcs);
    g.method = clustmethods{resstruct.method};
    g.dipole_locs = opt_offon(resstruct.chbx_enabledipole+1);
    if resstruct.chbx_enabledipole
        g.weightSP = str2num(resstruct.ed_weightSP);
        g.weightDP = str2num(resstruct.ed_weightDP);
    else
        g.weightSP = 1;
        g.weightDP = 0;
    end
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

if ~isfield(IMA, 'precluster'), disp('pop_clusterIMAtemplates: Preclustering must be performed'); return; end;
templates = [templates IMA.precluster.templates(:,freqind1:freqind2)];
IMICindex = [IMICindex [repmat(sujnum,1,size(IMA.precluster.IMICindex,1))' IMA.precluster.IMICindex]]; %% add subject index
if strcmp(g.dipole_locs, 'on')
    dipsources = [dipsources IMA.precluster.dipsources];
else
    dipsources = struct;
end

[clustidx, distance]...
    = clusterIMAtemplates(templates, IMICindex, g.nclust, 'pcs', g.pcs, 'method', g.method, 'dipsources', dipsources, 'weightSP',g.weightSP ,'weightDP',g.weightDP);

STUDY.etc.IMA.clustidx = clustidx;
STUDY.etc.IMA.distance = distance;
[STUDY] = pop_savestudy( STUDY, ALLEEG, 'savemode','resave');
end


