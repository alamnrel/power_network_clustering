function T = hclust(pts, dm, link, nc, varargin)
%
% Given a distance matrix, or bare point coordinates in an Euclidian space
% perform hierarchical clustering into a given number of clusters.
%
% Author: Ilya Tyuryukanov
% Date of first version: 18 November 2015 
% Last revision: 21 November 2015

if isempty(pts)&&~isempty(dm)
  dst = squareform(dm, 'tovector');  % transform distance/similarity matrix to pdist() format
  Z = linkage(dst, link);
  % dendrogram(Z,40)
  T = cluster(Z,'maxclust',nc);
elseif ~isempty(pts)&&isempty(dm)
  T = clusterdata(pts, 'maxclust', nc, 'linkage', link, varargin{:});
else
  error('[hclust] Either 1st or 2nd argument should be empty.')
end

end