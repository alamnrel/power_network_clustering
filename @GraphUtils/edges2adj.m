function [idx_out, val] = edges2adj( adj, inc, linidx, NDEBUG)
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

if nargin<4, NDEBUG = false; end

assert(NDEBUG || islogical(linidx), [mfilename,':WrongPublicInput'],...
  '[%s] The only non-PFgraph input argument (i.e. linidx) must be logical', mfilename)
assert(NDEBUG || ~isempty(adj)&&~isempty(inc), [mfilename,':WrongPublicInput'],...
  '[%s] Both adjacency and incidence matrices of g should not be empty', mfilename);

ij_binary = inc(:, linidx);
[row, ~] = find(ij_binary);
n_rows = numel(row)/2;
assert(NDEBUG || Utils.isint(n_rows), [ mfilename, ':InconsistentValue'],...
  ['[%s] The selected columns of the incidence matrix supplied with the ',...
  'graph object contain the odd total number of nonzero elements (but it ',...
  'should be even for any proper incidence matrix)'], mfilename)

idx_out = reshape(row, 2, n_rows)';
if nargout > 1
  lin_idx = sub2ind(size(adj), idx_out(:,1), idx_out(:,2));
  val = adj(lin_idx);
  val = full(val);  
end

