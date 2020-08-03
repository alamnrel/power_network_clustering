function [adj_sep, C, ixs, eXp, pcut] = cc_info(adj, inc, vw, cutind)
% 
% Syntax: [adj_sep, C, ixs, eXp, pcut] = cc_info(g, cutind)
% 
% Purpose: Separate graph into a set of connected components according to 
%   the index vectors cutind into the ordered set of graph edges (defined   
%   by the incidence matrix of the input graph g) and return some 
%   characteristics of each connected component.
%
% Input:
%   adj: graph adjacency matrix
%   inc: graph incidence matrix
%   cutind: ith row is an index vector into the network branch list (or 
%     into the incidence matrix), for showing the branches that should be 
%     cut to disconnect the ith identified island.                   
% 
% Output:
%   adj_sep: disconnected adjacency matrix 
%   C: labels for each vertex denoting the number of connected componet it
%     belongs to
%   ixs: cell array of lists of vertex indices for each actual connected 
%     component (after the removal of separator edges) 
%   eXp(1,:): expansion ratios { sum(Pcut,k)/sum(Pij), i,j are the buses 
%     in the kth island, k: (1st island) .. (last island) } of all islands
%     Hammer47!

%   eXp(2,:): expansion ratios in terms of RatioCut { sum(Pcut,k)/|Vk|, 
%      Vk are the buses in the kth island, k: (1st island)..(last island) } 
%      of all islands 
%   eXp(3,:): number of buses in each island
%   eXp(4,:): similarity-based modularity as [1]
%   pcut: power flow disruption of each island in units of g.adj (usually 
%      in MW)
% 
% Author: Ilya Tyuryukanov
% Date of first version: 27 September 2016 
% 
% [1] Feng Z., Xu X., Yuruk N., Schweiger T.A.J. (2007) A Novel 
% Similarity-Based Modularity Function for Graph Partitioning. 
% In: Song I.Y., Eder J., Nguyen T.M. (eds) Data Warehousing and Knowledge 
% Discovery. DaWaK 2007. Lecture Notes in Computer Science, vol 4654. 
% Springer, Berlin, Heidelberg

% Find cutset and total power flow disruption
cutind(cutind==2) = 0;  % neglect inner branches of islands, take only cutsets
cut = any(cutind, 1);
cut_idx = GraphUtils.edges2adj(adj, inc, cut);

% Separate graph into connected components
cut_idx = [cut_idx; cut_idx(:,2), cut_idx(:,1)];  % adjust to the symmetry of adj
cut_lind = sub2ind(size(adj), cut_idx(:,1), cut_idx(:,2));
adj_sep = adj;
adj_sep(cut_lind) = 0;
[S, C] = GraphUtils.alecconncomp(adj_sep);
vol_t = sum(sum(adj,2));

% Calculate expansions for each connected component
eXp = zeros(3, S);
ixs = cell(1, S);
pcut = NaN(1, S);
for i = 1:S
  idx = C == i;  % find vertices of kth connected component
  cut_i = idx*inc;  % find line cutset for this component
  cut_i(cut_i==2) = 0;
  cut_i = logical(cut_i);
  vol_i = sum(sum(adj_sep(idx,idx)));  % INTERNAL volume of ith conn. comp.
  vol_e = sum(vw(idx));                % total cluster node weight
  [~, bnd_i] = GraphUtils.edges2adj(adj, inc, cut_i);
  bnd_i = sum(bnd_i);  % cut of ith connected component
  pcut(i) = bnd_i;
  eXp(1, i) = abs(bnd_i/vol_e);  % expansion (NCut) of ith connected component
  eXp(2, i) = abs(bnd_i/nnz(idx));  % RCut of ith connected component
  eXp(3, i) = nnz(idx);  % number of vertices in ith connected component
  eXp(4, i) = vol_i/vol_t - ((bnd_i+vol_i)/vol_t)^2; % similarity-based modularity
  eXp(5, i) = bnd_i;  % cut of ith connected component
  eXp(6, i) = vol_i;  % internal volume of ith connected component
  ixs{i} = sort(find(idx));
end

end