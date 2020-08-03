function [cutind, busind, valid, eXp, pcut, n_outl_gen, dif] = metis_info(g, map)

% 
% Purpose:
%   Having as input the output of Metis, hMetis and KaHIP creates the 
%   partitions of these algorithms.
%
% Input:
% g: a (sensible) PFgraph object.
% map: vertex that shows in which partition each vertex belongs to (the 
%   number of partitions starts from 0)
%
% Output:
% busind: ith row is an index vector into the network bus matrix for showing   
%   the buses belonging to the ith identified island 
% cutind: ith row is an index vector into the network branch list for showing 
%   the branches that should be cut to disconnect the ith identified island   
% valid: parameter that shows if the initial number of areas created are
% different than the number of connected components.
%   valid=0: bad results (partition is not valid).
%   valid=1: good results (partition is valid).
% eXp: expansions from the partitions created.
% n_outl_gen: number of coherent generators that do not belong to their
%   dominating connected component.
% dif: if a partition is non-contiguous sums up the number of vertices in
%   minor connected components.
% 

busind = sparse(map+1, 1:size(g.adj, 1), true);
cutind = g.cutset(busind);
np = size(cutind, 1);
if nargout>5
  [eXp, pcut, ~, ~, ~, ~, n_outl_gen] = g.ici_info( cutind );
else
  [~, ~, ~, eXp, pcut] = g.cc_info( cutind );
end

S = size(eXp, 2);   % number of connected components

if S==np
    valid=1;
    dif=0;
elseif S>np
    valid=0;  
    dif=sum( eXp(3,np+1:end) );
else
    valid = 0;
    dif = 0;
end
end