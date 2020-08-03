function [cutind, busind, T] = final_cutset(g, T, merg)
%
% cutind: ith row is an index vector into the network branch list for showing 
%   the branches that should be cut to disconnect the ith identified island
% busind: ith row is an index vector into the network bus matrix for showing   
%   the buses belonging to the ith identified island        
% 
% Author: Ilya Tyuryukanov
% Date of first version: 01 March 2017

% Project the partition of a contracted graph to the initial graph
if nargin > 2
  for i = 1:1:numel(merg)
    T = T(merg{i});
  end
end

% Form islands
m = size(g.adj, 1);
busind = sparse(T, 1:m, true);
cutind = g.cutset(busind);
end