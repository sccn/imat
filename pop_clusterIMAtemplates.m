%  Takes WT comodulation spectral templates as collected in IMA.precluster.templates 
%  and clusters to specified number of clusters
%
% [STUDY] = pop_clusterIMAtemplates(STUDY,varargin);
%
% 
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
%
% Example:
% >> [STUDY] = pop_clusterIMAtemplates(STUDY, 'freqlim', [3 120] ,'nclust', 3, 'pcs', 20);
%
% Example:
% >> [STUDY] = pop_clusterIMAtemplates(STUDY, 'freqlim', [3 30] ,'nclust', 6);
%
%
% INPUTS:
% STUDY - STUDY structure with information on stored IMA files
% freqlim - 'integer'   [min max] frequency limits to restrict clustering
%            to i.e. alpha frequency range [8 14] default whole frequency
%            range in IMA.freqlim
% nclust -- [integer] if not empty, will divide templates into nclust clusters
%                     by 'pdist','linkage','dendrogram' Matlab functions
% method -- ['corr', 'euc' or 'k-means'] for correlation or euclidean distance for pdist, or k-means clustering
% pcs -- number of principal IM spectral template dimensions to retain
%
%
% OUTPUT
% STUDY - STUDY structure with indexes of clusters (column 4 in matrix) and distance to cluster
% centroid saved in STUDY.etc.IMA.clustidx in the form [SJ IM ICs clusindx] and STUDY.etc.IMA.distance
%

function [STUDY] = pop_clusterIMAtemplates(STUDY, varargin);

clustmethods = {'corr' 'euc' 'kmeans'};
g = finputcheck(varargin, {'freqlim'        'integer'       []             []; ...
    'nclust'     'integer'     []         [];...
    'pcs'     'integer'     []         [10];...
    'method'      'string'       clustmethods             'kmeans'; ...
    'dipole_locs' 'struc' [] [];...
    'weightSP' 'integer' [] [1];...
    'weightDP' 'integer' [] [1];...
    }, 'inputgui');
if isstr(g), error(g); end;

if nargin == 2
    clsutmethodslist = {'Correlation' 'Euclidean' 'Kmeans'}; 
    
      uilist = {{'style' 'text' 'string' 'Method'} {'style' 'popupmenu'  'string' clsutmethodslist 'tag' 'method' 'value' 1 }...
              {'style' 'text' 'string' 'Number of clusters'} {'style' 'edit'  'string' ' ' 'tag' 'nclust'}...
              {'style' 'text' 'string' 'Number of PCs'} {'style' 'edit'  'string' ' ' 'tag' 'npcs'} {'style' 'text' 'string' 'Freq. limits (Hz)'} {'style' 'edit'  'string' ' ' 'tag' 'freqlim'}};
     
     ht = 3; wt = 4 ;
     geom = {{wt ht [0 0]  [1 1]} {wt ht [1 0]  [1 1]}...
             {wt ht [0 1]  [1 1]} {wt ht [1 1]  [1 1]}...
             {wt ht [0 2]  [1 1]} {wt ht [1 2]  [1 1]} {wt ht [2 2]  [1 1]} {wt ht [2.8 2]  [1 1]}};
     
    [result, ~, ~, resstruct, ~] = inputgui('title','Cluster IM templates -- pop_clusterIMAtemplates', 'geom', geom, 'uilist',uilist, 'helpcom','pophelp(''pop_clusterIMAtemplates'');');

    g.freqlim = str2num(resstruct.freqlim);
    g.nclust = str2num(resstruct.nclust);
    g.pcs = str2num(resstruct.nclust);
    g.method = clustmethods{resstruct.method};
end

subjcode = STUDY.subject; 
    
templates = [];
IMICindex = [];
dipsources = [];
scalpmaps = [];
    

for iko = 1:length(subjcode)
indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));

%% load IMA file for curent subject
load([STUDY.datasetinfo(indsj(1)).filepath filesep STUDY.etc.IMA.imafilename{iko,:}], '-mat' );


str = string(STUDY(iko).subject);
sujnum = sscanf(str,'S%d');


if isempty(g.freqlim)
   g.freqlim = IMA.freqlim;
end

[val, freqind1] = min(abs(IMA.freqvec-g.freqlim(1)));
[val, freqind2] = min(abs(IMA.freqvec-g.freqlim(2)));


templates = [templates IMA.precluster.templates(:,freqind1:freqind2)];
IMICindex = [IMICindex [repmat(sujnum,1,size(IMA.precluster.IMICindex,1))' IMA.precluster.IMICindex]]; %% add subject index
if dipole_locs = 'on'
dipsources = [dipsources IMA.precluster.dipsources];
end

[clustidx, distance]...
    = clusterIMAtemplates(templates, IMICindex, g.nclust, 'pcs', g.pcs, 'dipole_locs', g.dipsources, 'weightSP',g.weightSP ,'weightDP',g.weightDP);

STUDY.etc.IMA.clustidx = clustidx;
STUDY.etc.IMA.distance = distance;
[STUDY] = pop_savestudy( STUDY, ALLEEG, 'savemode','resave');
end


