% subclusters IM template cluster based on location of dipoles
%
% [STUDY] = pop_subclusterIMAtemplates(STUDY,varargin);
%
% %
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2019
% adapted from a function written by Julie Onton
%
% Example:  create 2 subclusters of cluster 1 and 2
% >> [STUDY] = pop_subclusterIMAtemplates(STUDY, 'clust', [1 2], 'nclust', 2);
%
%
%
% INPUTS:
% STUDY - STUDY structure with information on stored IMA files
% clust -- [integer] culster indices to subcluster - if empty subcluster all clusters 
% nclust -- [integer] will divide templates into nclust clusters
%                     by 'pdist','linkage','dendrogram' Matlab functions -
%                     default 2
% method -- ['corr', 'euc' or 'k-means'] for correlation or euclidean distance for pdist, or k-means clustering
%            default 'k-means'
%
% OUTPUT
% STUDY - STUDY structure with indexes of subclusters (column 5 in matrix) and distance to subcluster
%         centroid saved in STUDY.etc.IMA.clustidx in the form [SJ IM ICs clusindx subclusindx] 
%         and STUDY.etc.IMA.subclustdistance [ICs x Clusters]
%

function [STUDY] = pop_subclusterIMAtemplates(STUDY,varargin);


g = finputcheck(varargin, {'clust'      'integer'     []         [unique(STUDY.etc.IMA.clustidx(:,4))'];...
    'nclust'     'integer'     []         [2];...
    'method'      'string'       {'corr' 'euc' 'kmeans'}             'kmeans'; ...  
    }, 'inputgui');
if isstr(g), error(g); end;


subjcode = STUDY.subject; 

dipsources = [];
scalpmaps = [];
    

for iko = 1:length(subjcode)
indsj = find(ismember({STUDY.datasetinfo.subject}, STUDY(iko).subject));

%% load IMA file for curent subject
load([STUDY.datasetinfo(indsj(1)).filepath '/' STUDY.etc.IMA.imafilename(iko,:)], '-mat' );

str = string(STUDY(iko).subject);
sujnum = sscanf(str,'S%d');

dipsources = [dipsources IMA.precluster.dipsources];

end



[clustidx, subclustdistance]...
        = subclusterIMAtemplates(STUDY.etc.IMA.clustidx, dipsources, 'clust', g.clust, 'nclust', g.nclust, 'method', g.method);

    
STUDY.etc.IMA.clustidx = clustidx;
STUDY.etc.IMA.subclustdistance = subclustdistance;
[STUDY] = pop_savestudy( STUDY, 'savemode','resave');

end
