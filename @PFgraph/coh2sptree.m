function trees = coh2sptree(g)
% 
% This method aims to find subnetworks connecting coherent generators in
% each group with each other. For that a shortest path algorithm is used.
% The graph adjacency matrix thus should represent a distance matrix.
% 
% The test condition for success is that every returned subnetwork is
% connected (while containing all requested generator buses). If this
% condition is violated, the method terminates. No restarts are used. 
% The underlying problem is, in fact, NP-complete.
% 
% This method just represents one simple solution attempt (i.e. the method
% is simple by its nature and e.g. intermediate restarts most likely won't
% fix the fundamental issues).
%
% Author: Ilya Tyuryukanov
% Date of first version: 01 February 2016
% Last revision: 23 May 2016
% 
  
adj = g.adj;
adj_bin = double(logical(adj));
adj_dist = adj_bin;  % the alternative GraphUtils.simadj2distadj( adj, 'inv2' ) minimizes the "collapsed" power flow, but not the "collapsed" number of links -> less robust
siz = size(adj_dist);
coh_ml = g.coh;
epsilon = eps('double');

% Set existing pairwise must-links to be the most preferrable connectors
if ~isempty(g.ml)
  bran_ml = g.ml;
  bran_ml = [bran_ml; bran_ml(:,2), bran_ml(:,1)];
  idx = sub2ind(siz, bran_ml(:,2), bran_ml(:,1));
  adj_dist(idx) = 10*epsilon;
end

trees = {};   % coherency subnetworks will be stored as groups of nodes
coh_ml_id = unique(coh_ml(2, :));
adj_bin_tmp = adj_bin;
for i = coh_ml_id
  idx = coh_ml(2, :) == i;
  bus = coh_ml(1, idx);
  [~, hub] = max(sum(adj_bin(bus, :), 2));
  hub_bus = bus(hub);
  rem_bus = setdiff(bus, hub_bus);
  
  % Exclude the nodes of other coherent generator groups from the search 
  rem_idx = ~idx;
  rem_coh = coh_ml(1, rem_idx);
  adj_bin_tmp(rem_coh,:) = 0;
  adj_bin_tmp(:,rem_coh) = 0;
  if ~isempty(rem_bus)
    [dist, paths] = dijkstra(adj_bin_tmp, adj_dist, hub_bus, rem_bus, false);
    assert(~any(isinf(dist)), [ mfilename, ':AlgorithmNotConverged'],...
      ['[%s] The previously taken connectors have cut all possible links ',...
      'between the currently considered coherent generator buses %d and ',...
      '%d (the shown bus indices may be from a reduced graph)'],...
      mfilename, hub_bus, rem_bus(find(isinf(dist),1))); 
    
    % Remember the identified connectors
    paths_i = [paths{:}];
    used_vs = unique(paths_i(:));  % all used vertices (of the ith subnetwork)
  else
    used_vs = hub_bus;
  end  
  trees = [trees; used_vs];
  
  % Isolate used vertices in order to ban any intersection of the ith group
  %  with the following groups
  adj_bin(used_vs,:) = 0;
  adj_bin(:,used_vs) = 0;
  adj_bin_tmp = adj_bin;  % update adj_bin_tmp
end
