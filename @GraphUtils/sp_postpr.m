function [T] = sp_postpr(Y, postpr, alg, nc, adj)
%
% Cluster the rows of the eigenvector matrix Y obtained from spectral 
% clustering (or effectively of any real matrix...) into nc clusters with 
% a postprocessing method postpr, for a "point representation" alg.
% If points are to be presented as a graph in the spectral embedding Y,
% then the graph connectivity structure in the form of adjacency matrix adj 
% is required.
% 
% Author: Ilya Tyuryukanov
% Date of first version: 18 August 2017
% Last revision: 5 November 2017

hlink = {'single', 'weighted', 'average', 'complete', 'ward', 'centroid',...
    'median'};   % "hierarchical" linkages
slink = {'pic'};   % "similarity" linkages
plink = {'kmeans', 'kmedoids'};   % point clustering in Euclidian space 
[nm, pts_repr] = parsealg(alg);

switch pts_repr
  case 'eucld'
    if ismember({postpr},hlink)
      T = GraphUtils.hclust(Y, [], postpr, nc);
    elseif ismember({postpr},plink)
      f = str2func(postpr);
      s = rng;
      rng_back = onCleanup(@() rng(s));
      rng(5);  % for reproducibility
      T = f(Y, nc, 'Replicates', 20);  % it's possible to put some other/additional options as well
    end  
  case 'embgr'
    if ismember({postpr},hlink)
      [~, dm] = GraphUtils.buildDistanceGraph(Y, adj, nm);
      T = GraphUtils.hclust([], dm, postpr, nc);
    elseif strcmp(postpr, 'pic')
      pts.Y = Y; pts.nm = nm;
      T = gacCluster(pts, adj, nc, 'path');
    end    
  otherwise
    error([mfilename,':WrongPublicInput'],...
    '[%s]. The 2nd input can only take values ''eucld'' and ''embgr''.',...
    mfilename);
end
end

function [nm, pts_repr] = parsealg(alg)
% Parse string with algorithm type and return its meaningful parts
% (normalization mode and points representation).
%
% Author: Ilya Tyuryukanov
% Last revision: 14 January 2016

assert(numel(alg)==4, [ mfilename, ':WrongPublicInput'],...
  '[%s] The alg string can be only four characters long', mfilename);

nm = alg(1:3);
if strcmp(alg(4), '+')
  pts_repr = 'eucld';  % just as points in Euclidian space
elseif strcmp(alg(4), '-')
  pts_repr = 'embgr';  % as 'embedded graph'
else
  error('[parsealg] The 4th character of the alg string can be either "+" or "-"');
end
end
