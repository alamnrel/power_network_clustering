function allSP = allShortestPaths(adj, mode)
% Syntax: allSP = allShortestPaths(adj, mode)
%
% Purpose:
% Find all shortest paths in a graph (using its adjacency matrix)
%
% Input:
% adj - graph (weighted) adjacency matrix. adj is sparse symmetric real matrix.
% mode - the way function works
%  'inv': first convert weighted graph to an 'inverse' weghts  graph (i.e.
%    represent a tight coupling with a SMALL weight) and then find all
%    shortest paths. This option is useful if edge weights represent a
%    similarity measure between the graph nodes.
%  'orig': find all shortest paths in the original graph. This option is
%    useful if edge weights represent a dissimilarity measure (i.e. a
%    "distance") between the graph nodes.
%
% Output:
% allSP - matrix of all shortest paths between vertices
%
% Author: Ilya Tyuryukanov
% Date: 30 May 2015
% 25 November 2015: ??
% 29 June 2016: shortest paths are calculated with dijkstra.m, as the input
%   is not converted to a triangular matrix inside of simadj2distadj.m, and
%   it also makes the overall implementation more consistent.
%   

adj_dist = GraphUtils.simadj2distadj( adj, mode );
% adj_bin = double(logical(adj));

allSP = graphallshortestpaths(tril(adj_dist), 'Directed', false);
% allSP = dijkstra(adj_bin, adj_dist);
end