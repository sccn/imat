%  Takes WT comodulation spectral templates as collected in IMA.precluster.templates 
%  and clusters to specified number of cluster
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
% pcs -- number of dimensions to reduce the spectra to
%
%
% OUTPUT
% STUDY - STUDY structure with indexes of clusters (column 4 in matrix) and distance to cluster
% centroid saved in STUDY.etc.IMA.clustidx in the form [SJ IM ICs clusindx] and STUDY.etc.IMA.distance
%

function [STUDY] = pop_clusterIMAtemplates(STUDY,varargin);


g = finputcheck(varargin, {'freqlim'        'integer'       []             []; ...
    'nclust'     'integer'     []         [];...
    'pcs'     'integer'     []         [10];...
    'method'      'string'       {'corr' 'euc' 'kmeans'}             'kmeans'; ...  
    }, 'inputgui');
if isstr(g), error(g); end;

subjcode = STUDY.subject; 


    
templates = [];
IMICindex = [];
dipsources = [];
scalpmaps = [];
    

for iko = 1:length(subjcode)
indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));

%% load IMA file for curent subject
load([STUDY.datasetinfo(indsj(1)).filepath '/' STUDY.etc.IMA.imafilename(iko,:)], '-mat' );


str = string(STUDY(iko).subject);
sujnum = sscanf(str,'S%d');


if isempty(g.freqlim)
   g.freqlim = IMA.freqlim;
end

[val, freqind1] = min(abs(IMA.freqvec-g.freqlim(1)));
[val, freqind2] = min(abs(IMA.freqvec-g.freqlim(2)));


templates = [templates IMA.precluster.templates(:,freqind1:freqind2)];
IMICindex = [IMICindex [repmat(sujnum,1,size(IMA.precluster.IMICindex,1))' IMA.precluster.IMICindex]]; %% add subject index

end

[clustidx, distance]...
    = clusterIMAtemplates(templates, IMICindex, g.nclust, 'pcs', g.pcs, 'method', g.method);

STUDY.etc.IMA.clustidx = clustidx;
STUDY.etc.IMA.distance = distance;
[STUDY] = pop_savestudy( STUDY, 'savemode','resave');


end


