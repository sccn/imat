%  Takes WT comodulation spectral templates and clusters to specified
%  number of cluster
%  
%
% [clustidxs, distance] = clusterIMAtemplates(templates, IMICindex, nclust, varargin);
%
%
% Author: Johanna Wagner, Swartz Center for Computational Neuroscience, UC
% San Diego, 2020
% adapted from a function written by Julie Onton
%
% Example:
% >> [clustidx, distance] = clusterIMAtemplates(templates, IMICindex, 3, 'pcs', 15, 'method', 'kmeans');
%
%
% INPUT
% required Inputs:
% templates -- [matrix] of (templates x freqs). Each row is a single modulator spectral template
%                       as collected over subjects from IMA.precluster.templates
% IMICindex -- matrix of indexes of templates matrix [subject IM IC] as
%              collected in pop_collecttemplates and pop_clusterIMAtemplates
% nclust -- [integer] if not empty, will divide templates into nclust clusters
%                     by 'pdist','linkage','dendrogram' Matlab functions
%
% optional Inputs:
% pcs -- number of principal IM spectral template dimensions to retain
% dipsources -- dipole locations as collected in pop_clusterIMAtemplates
% weightDP -- weight to assign to dipole locations for clustering
% weightSP -- weight to assign to spectra for clustering
% method -- ['k-means']  k-means clustering

% OUTPUT
% clustidx -- matrix of indexes of [cluster subject IM IC] 
% distance -- distance of each IC to each cluster centroid [ICs x Clusters]

function [clustidx, distance] = clusterIMAtemplates(templates, IMICindex, nclust, varargin);

g = finputcheck(varargin, {'pcs'        'integer'     []         [10];...
                          'method'      'string'       {'kmeans'}             'kmeans'; ...  
                          'dipsources'  'struct' [] struct;...
                          'weightSP'    'integer' [] [1];...
                          'weightDP'    'integer' [] [1];...
                           }, 'inputgui');
if isstr(g), error(g); end;

%'method'      'string'       {'corr' 'euc' 'kmeans'}             'kmeans'; ... 
  
distance = [];

%% reduce dimensions
fprintf('\nPCA''ing to %s dimensions.\n',int2str(g.pcs));%%%%%%%%%%%%%%%%%%%%%%
[U,S,V] = svds(templates,g.pcs);% if you scale 'acts', you can't scale the 'eigvec'
pc = (U*S); % scale 'activations' for appropriate weighting in decomp of pc
eigvec = V;



if isfield(g.dipsources,'posxyz')
    dipole_locs = reshape([g.dipsources.posxyz]',3,size(pc,1))';
    featurevec = [g.weightSP*zscore(pc) g.weightDP*zscore(dipole_locs)];
else
    featurevec = zscore(pc);
end


if isempty(nclust)
    nclust = (size(IMICindex,1)/length(unique(IMICindex(:,1))))/3;
end;
% if strcmp(g.method,'corr')
%     alldist = pdist(featurevec, 'correlation'); % euc better than seuc
%     links = linkage(alldist,'complete');
%     figttl = 'Template clusters -- correlation method';
%     figure;[hnd,idx,perm]=  dendrogram(links,nclust);%close
% elseif strcmp(g.method,'euc')
%     alldist = pdist(featurevec, 'euclidean'); % euc better than seuc
%     links = linkage(alldist,'ward');
%     figttl = 'Template clusters -- euclidean method';
%     figure;[hnd,idx,perm]=  dendrogram(links,nclust);%close
% elseif strcmp(g.method,'kmeans')
    [idx,C,sumd,distance] = kmeans(featurevec,nclust,'replicates',10,'emptyaction','drop');
    figttl = 'Template clusters -- kmeans method';
%end;              

clustidx = [IMICindex idx];



end

