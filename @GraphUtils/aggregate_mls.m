function ixs = aggregate_mls( bran_ml, coh_ml )
% aggregate_mls This function aggregates must-links into pairwise must-link 
%   connected components. This is in particular useful to do before initiating 
%   the merging of pairwise must-links in order to aggregate them into largest
%   possible components before the merging and thus optimize the merging.
%   This also helps to avoid repeated node indices during the merging stage 
%   which may cause problems during actual vertex contractions (the same 
%   vertex can be merged only once).
%
%   bran_ml should be a two column matrix of whole numbers, with the pairs
%   of indices of buses to be merged in each column indicating the lines to
%   be collapsed (i.e. the same format as the PFgraph class 'ml' property).
%
%   coh_ml should be a two row matrix of whole numbers, in the same format
%   as the PFgraph class 'coh' property.
%
%   coh_ml input can be included ONLY in order to place the coherent generator 
%   subnetworks at the end of ixs.
%
% Author: Ilya Tyuryukanov
% Date of first version: 09 Sep 2015
% Last revision: 02 February 2016 (convertion to a PFgraph class static method)

switch nargin
  case 1
    n_bus = max(bran_ml(:));
  case 2
    n_bus = max([max(bran_ml(:)), max(coh_ml(1,:))]);
end

if nargin > 1 && ~isempty(coh_ml) && Utils.isint(coh_ml)
  bias = min(coh_ml(2,:)) - 1; % ensure that counting starts from one (just in case!)
  coh_ml(2,:) = coh_ml(2,:) - bias;
  coh_ml_bins = histcounts(coh_ml(2,:), 'BinMethod', 'integers');  % number of must-link buses in each grpML-group
  sngl_gen = find(coh_ml_bins == 1);
  if ~isempty(sngl_gen)
    for i = 1:numel(sngl_gen)
      idx = coh_ml(2,:) == sngl_gen(i);
      sngl_gen(i) = coh_ml(1, idx);  % save single generators
      coh_ml(:, coh_ml(1,:) == sngl_gen(i)) = [];  % remove single generators
    end
    coh_ml_bins(:, coh_ml_bins<2) = 0;  % update 'coh_ml_bins' to correspond to 'coh_ml'
  else
    sngl_gen = [];
  end
  
  gens_all = sort([sngl_gen(:); coh_ml(1,:)'], 'ascend');  % all generator bus indices
  gens = sparse(1, gens_all, true, 1, n_bus);
  gens_assigned = false(1, n_bus);  % gens (among all buses) assigned to any ML-node [ixs1{}]
  
  ml = NaN(sum(coh_ml_bins)-size(coh_ml_bins, 2), 2);  % NaN is for easy bug discovery
  enD = 0;
  for i = unique(coh_ml(2,:))
    cohIDX = coh_ml(1, coh_ml(2,:) == i)';  % matrix indices of coherent buses
    ml_part = [repmat(cohIDX(1), numel(cohIDX)-1, 1), cohIDX(2:end)];
    start = enD + 1;
    enD = start + coh_ml_bins(i) - 2;
    ml(start:enD, :) = ml_part;
  end
else
  ml = []; 
  sngl_gen = [];
  gens = false(1, n_bus);
  gens_assigned = false(1, n_bus);
  gens_all = [];
end
ml = [ml; bran_ml];

adj_ml = sparse(ml(:,1), ml(:,2), 1, n_bus, n_bus) + ...
  sparse(ml(:,2), ml(:,1), 1, n_bus, n_bus);

% Find connected components
[~, C] = GraphUtils.alecconncomp(adj_ml);  % alecconncomp is fast
large_idx = Utils.large_comp(C, 1);

% Assign must-link groups
num_ml = numel(large_idx) + numel(sngl_gen);  % max. number of must-link groups
ixs1 = cell(num_ml, 1);  % must-link groups containing generators
ixs2 = cell(num_ml, 1);  % must-link groups with no generators
n_ixs1 = 0;
n_ixs2 = 0;
for i = 1:1:numel(large_idx)
  large_nodes = C == large_idx(i);
  gens_ccmp_i = large_nodes & gens;
  gens_assigned = gens_assigned | gens_ccmp_i;
  is_GenCmpt = any(gens_ccmp_i);
  if is_GenCmpt
    n_ixs1 = n_ixs1 + 1;
    ixs1{n_ixs1} = find(large_nodes)';
  else
    n_ixs2 = n_ixs2 + 1;
    ixs2{n_ixs2} = find(large_nodes)';
  end
end
gens_assigned = find(gens_assigned);  % back to matrix indices
gens_separate = setdiff(gens_all, gens_assigned);
if ~isempty(gens_separate)
  for i = 1:1:numel(gens_separate)
    n_ixs1 = n_ixs1 + 1;
    ixs1{n_ixs1} = gens_separate(i);  % just a single node to "merge with itself"...
  end
end
ixs1 = ixs1(1:n_ixs1);
ixs2 = ixs2(1:n_ixs2);
ixs = [ixs2; ixs1];  % last n_ixs1 groups are the generator groups

end

