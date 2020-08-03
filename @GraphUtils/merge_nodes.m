function [adj_out, vw_out, map, ixsc] = merge_nodes(adj_in, ixs, vw_in)
%
% This method compresses the original nodes information and returns new
% adjacency matrix adj_out and new node weights vw_out. The second input
% can be set to [], which will result in each input vertex having the unity
% weight.
% 
% Original authors: (C)2012 Syama Sundar Rangapuram and Matthias Hein
% Date of first version: 28 January 2015
% Last revision: 01 February 2015

if isempty(ixs) || isempty(ixs{1})
  adj_out = adj_in;
  vw_out = vw_in;
  return;
end

n = size(adj_in, 1);
if nargin<3 || isempty(vw_in)
  vw_in = ones(n, 1);
end

n_deleted = 0;  % number of vertices will be reduced by this amount
ixs_all = [];  % ALL must-link constrained vertex indices
for i=1:length(ixs)
  n_deleted = n_deleted + length(ixs{i}) - 1;
  ixs_pls = ixs{i};
  ixs_all = [ixs_all; ixs_pls(:)];
end

ixsc = setdiff(1:n, ixs_all);  % "unconstrained" vertices
n_out = n - n_deleted;  % final number of vertices
n_ixsc = length(ixsc);  % number of unconstrained vertices
adj_out = sparse(n_out, n_out);  % final adjacency matrix

vw_out = zeros(n_out, 1);  % final vertex weights
vw_out(1:n_ixsc) = vw_in(ixsc);  % set the first n_ixsc vertex weights entries to those of unconstrained vertices

map = zeros(n,1);
map(ixsc) = 1:n_ixsc;  % set the mapping of unconstrained vertices to the first n_ixcs vertices
for i=1:length(ixs)
  map(ixs{i}) = n_ixsc + i;  % set the mapping of constrained vertices to the last length(ixs) vertices
end

all_idx = 1:n;  % old indices
all_idx_out = map(all_idx);  % to which new index an old index belongs..

for i=1:length(ixs)
  subset = ixs{i};
  tot_weight = 0;
  for index=1:length(subset)
    tot_weight = tot_weight + vw_in(subset(index)); % weight/degree of merged index
  end
  vw_out(n_ixsc+i) = tot_weight;  % weight/degree of merged index
  
  % Sum the weights from ALL vertices in 'subset' (merged as 1 vertex)
  % to each vertex which is the neighbour of AT LEAST ONE 'subset' vertex
  adj_out(:, n_ixsc+i) = sparse( all_idx_out, 1, sum(adj_in(:,subset), 2), n_out, 1);  % this line is naturally slow because it also performs accumarray() with all_indices_p
end

adj_out = [adj_in(ixsc, ixsc) adj_out(1:n_ixsc, n_ixsc+1:end); adj_out(:, n_ixsc+1:end)'];
adj_out(1:n_out+1:n_out*n_out) = 0;  % set the diagonal entries to zero
adj_out = (adj_out+adj_out.')/2;  % ensure symmetry
end