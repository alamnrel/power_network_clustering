function nSP = nShortestPaths(adj, gcg_vx, mode)
% Syntax: allSP = nShortestPaths(adj, gcg_vx, mode)
%
% Purpose:
% Find shortest paths in a graph (using its adjacency matrix) from 
%   its gcg_vx vertices to all other vertices. 
%
% Input:
% adj - graph (weighted) adjacency matrix. adj is a sparse symmetric real 
%   matrix.
% gcg_vx - sink vertices to calculate the paths to (they usually represent
%   generator coherent groups hence the name).
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
% nSP - matrix of shortest paths from sink vertices to all other vertices
%   It should have the following size: size(adj,1) x |gcg_vx|
%
% Author: Ilya Tyuryukanov
% Date: 30 May 2015
% 28 June 2016: Now the 2nd input represents an arbitrary set of vertex
%   indices, not exclusively the last ones. BioInformatics Toolbox is not
%   required anymore as now shortest paths are computed with dijkstra.m
%
old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath( fullfile(fileparts(fileparts(mfilename('fullpath'))), 'third_party') )

adj_dist = GraphUtils.simadj2distadj( adj, mode );
adj_bin = double(logical(adj));

m = size(adj_bin, 1);
nSP = dijkstra(adj_bin, adj_dist, 1:1:m, gcg_vx, false);
end
