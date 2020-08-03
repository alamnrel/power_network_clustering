function g_cntrct = merge_nodes(g, ixs)
%
% This method compresses the original nodes information and returns a new 
% PFgraph object g_cntrct. 
% If it's undesirable (for any reason, e.g. in order to save time) to
% update the coherency or must-link branch information, the corresponding
% fields of the input graph object should be set to [].
% 
% Original authors: (C) 2012 Syama Sundar Rangapuram and Matthias Hein
% Date of first version: 28 January 2015
% Last revision: 01 February 2015

if isempty(ixs) || isempty(ixs{1})
  g_cntrct = g;
  return;
end

adj_in = g.adj;
vw_in = g.vw;
[adj_out, vw_out, map, ixsc] = GraphUtils.merge_nodes(adj_in, ixs, vw_in);
n_out = size(adj_out, 2);
n_ixsc = length(ixsc);

coh_in = g.coh;
if ~isempty(coh_in) && Utils.isint(coh_in)
  % Correct the indices of the unmerged gens for the reduced graph
  [~, ~, idx_free_gens] = intersect(ixsc, coh_in(1, :));  % identify generators which haven't been merged
  coh_out = coh_in(:, idx_free_gens);
  coh_out(1, :) = map(coh_out(1,:));
  
  % Add merged gens to coh_out for the reduced graph
  ngen_unmerged = numel(idx_free_gens);
  k_merged_gens = 1;
  for i = 1:1:length(ixs)
    ixs_curr = ixs{i};
    [~, ~, idx_merged_gens] = intersect(ixs_curr, coh_in(1,:));
    if ~isempty(idx_merged_gens)
      group_label = coh_in(2, idx_merged_gens);
      assert(all(group_label == group_label(1)),...
        [ mfilename, ':InconsistentResults'], ['[%s] The ',...
        'generator buses in the %d merged bus don''t belong to the same group'],...
         mfilename, i);
      coh_out(1, ngen_unmerged + k_merged_gens) = n_ixsc + i;  % the number of contracted bus as it was
      coh_out(2, ngen_unmerged + k_merged_gens) = group_label(1);
      k_merged_gens = k_merged_gens + 1;
    end
  end
  
  % Set the output generator vector
  gen_out = false(1, n_out);
  gen_out(coh_out(1,:)) = true;
else
  coh_out = [];
  gen_out = [];
end

% Find must-link branches remaining inside of the merged nodes (to delete)
% and rename the rest
ml_in = g.ml;
if ~isempty(ml_in) && Utils.isint(ml_in)
  ml_islands = zeros(length(ixs), size(adj_in, 1));  % represent each group to collapse with its node indicator vector
  for i = 1:1:length(ixs)
    ixs_curr = ixs{i};
    ml_islands(i, ixs_curr) = 1;
  end
  cut = cutset(g, ml_islands);
  cut_int = cut; 
  cut_int(cut_int == 1) = 0; cut_int = any(cut_int, 1);
  ml2delete = g.edges2adj( cut_int );
  ml2delete = [ml2delete; ml2delete(:,2), ml2delete(:,1)];  % account for any possible order of buses
  [~, idx2delete, ~] = intersect(ml_in, ml2delete, 'rows');
  ml_in(idx2delete, :) = [];
  ml_out(:,1) = map(ml_in(:,1));
  ml_out(:,2) = map(ml_in(:,2));
  idx2delete = (ml_out(:,1) == ml_out(:,2));  % these are mapped to a single node
  ml_out(idx2delete,:) = [];  
  ml_out = sort(ml_out, 2, 'ascend');
  ml_out = unique(ml_out, 'rows');
else
  ml_out = [];
end

g_cntrct = PFgraph('adj', adj_out, 'merge_map', map, 'vw', vw_out,...
  'coh', coh_out, 'ml', ml_out, 'gen', gen_out);
end