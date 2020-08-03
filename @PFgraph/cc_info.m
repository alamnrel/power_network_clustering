function [adj_sep, C, ixs, eXp, pcut] = cc_info(g, cutind)
% 
% Syntax: [adj_sep, C, ixs, eXp, pcut] = cc_info(g, cutind)
% 
% Purpose: Separate graph into a set of connected components according to 
%   the index vectors cutind into the ordered set of graph edges (defined   
%   by the incidence matrix of the input graph g) and return some 
%   characteristics of each connected component.
%
% Input:
%   g: a PFgraph object (incidence and adjacency matrices + constraints)
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
%   eXp(2,:): expansion ratios in terms of RatioCut { sum(Pcut,k)/|Vk|, 
%      Vk are the buses in the kth island, k: (1st island)..(last island) } 
%      of all islands 
%   eXp(3,:): number of buses in each island
%   pcut: power flow disruption of each island in units of g.adj (usually 
%      in MW)
% 
% Author: Ilya Tyuryukanov
% Date of first version: 27 September 2016 

adj = g.adj;
inc = g.inc;
vw = g.vw;

if nargout < 3
  [adj_sep, C] = GraphUtils.cc_info(adj, inc, vw, cutind);
else
  [adj_sep, C, ixs, eXp, pcut] = GraphUtils.cc_info(adj, inc, vw, cutind);
end

end