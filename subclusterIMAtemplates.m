%
% Takes clusters of templates and subclusters them based on diploe locations to specified number of subclusters
%
% [subclustidxs, distance]...
%        = subclusterIMAtemplates(clustidx, dipsources, varargin);
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC San Diego, 2020 
% adapted from a function written by Julie Onton
%
%
% Example:
% [clustidx, subclustdistance]...
%        = subclusterIMAtemplates(STUDY.etc.IMA.clustidx, dipsources,...
%                                 'clust', [1 3], 'nclust', [2], 'method', 'kmeans');
%
%
% 
% INPUT:
% required Inputs:
% clustidx -- matrix of indexes of [subject IM IC]
% dipsources -- [matrix] of (dipoles) as collected over subjects from IMA.precluster.dipsources
%        
%
% optional Inputs:
% clust -- [integer] culster indices to subcluster - if empty subcluster all clusters 
% nclust -- [integer] will divide templates into nclust clusters
%                     by 'pdist','linkage','dendrogram' Matlab functions -
%                     default 2
% method -- ['corr', 'euc' or 'k-means'] for correlation or euclidean distance for pdist, or k-means clustering
%            default 'k-means'
%
% OUTPUT:
% subclustidx -- matrix of indexes of subclusters (column 5 in matrix) [subject IM IC cluster subcluster] 
% subclusdistance -- distance of each IC to each subcluster centroid [ICs x Clusters]


function [clustidx, subclusdistance]...
        = subclusterIMAtemplates(clustidx, dipsources, varargin);

g = finputcheck(varargin, {'clust'      'integer'     []         [unique(clustidx(:,4))'];...
    'nclust'     'integer'     []         [2];...
    'method'      'string'       {'corr' 'euc' 'kmeans'}             'kmeans'; ...  
    }, 'inputgui');
if isstr(g), error(g); end;


% create variables with zeros    
subclusdistance = [zeros(size(clustidx,1),g.nclust)];
subclusidx = [clustidx(:,1:4) zeros(size(clustidx,1),1)];

for dind = 1:length(g.clust)

indclus = find(clustidx(:,4) == g.clust(dind));
dipole_locs = reshape([dipsources.posxyz]',3,size(clustidx,1))';

featurevec = zscore(dipole_locs(indclus',:));

if size(featurevec,1) > g.nclust

if strcmp(g.method,'corr')
    alldist = pdist(featurevec, 'correlation'); % euc better than seuc
    links = linkage(alldist,'complete');
    figttl = 'Template clusters -- correlation method';
    figure;[hnd,idx,perm]=  dendrogram(links,g.nclust);%close
elseif strcmp(g.method,'euc')
    alldist = pdist(featurevec, 'euclidean'); % euc better than seuc
    links = linkage(alldist,'ward');
    figttl = 'Template clusters -- euclidean method';
    figure;[hnd,idx,perm]=  dendrogram(links,g.nclust);%close
elseif strcmp(g.method,'kmeans')
    [idx,C,sumd,distancetmp] = kmeans(featurevec,g.nclust,'replicates',10,'emptyaction','drop');
    figttl = 'Template clusters -- kmeans method';
end;              

subclusidx(indclus,5) = idx;
subclusdistance(indclus,:) = distancetmp;

else
    fprintf(['\n \n cannot subcluster cluster ' num2str(g.clust(dind)) ' since number of ICs contained in cluster ' num2str(g.clust(dind)) ' \n smaller than number of subclusters \n \n'])
end
end
clustidx = subclusidx;
   
    
    
