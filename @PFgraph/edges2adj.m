function [idx_out, val] = edges2adj( g, linidx, NDEBUG)
% 
% edges2adj This method takes line indices (e.g. line1, line3, line5 would 
%  correspond to [1 0 1 0 1], the ordering is according to the incidence
%  matrix) and uses the correspondence between  inc and adj properties of a 
%  PFgraph object to return the corresponding pairs of bus indices for each 
%  line.
% 
%  Note that no "duplicated" index pairs with just swapped order of indices  
%  are returned (i.e. [i,j] is returned, but no [j,i] in addition to that),   
%  so if all the entries of a symmetric adjacency matrix need to be accessed,
%  these additional swapped pairs should be appended to the output of the 
%  method.
%
% NDEBUG is to switch off assertions in order to increase the speed (e.g. 
%   for many repeated calls). Switch off assertions by setting NDEBUG to
%   true.
%
% Author: Ilya Tyuryukanov
% Date of first version: 01 February 2016
% Last revision: 03 February 2016

if nargin<3, NDEBUG = false; end

adj = g.adj;
inc = g.inc;

if nargout > 1
  [idx_out, val] = GraphUtils.edges2adj( adj, inc, linidx, NDEBUG);
else
  idx_out = GraphUtils.edges2adj( adj, inc, linidx, NDEBUG);
end

% % % Convert to nominal bus numbers from bus numbers induced by adj. matr.
% % idx_out(:,1) = g.bus(idx_out(:,1));
% % idx_out(:,2) = g.bus(idx_out(:,2));
