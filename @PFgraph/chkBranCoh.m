function [ok, bran_viols, coh_viols] = chkBranCoh(g)
% 
% Each Must-Link group in g.coh is Cannot-Link to all other Must-Link groups
% Therefore:
% ok = true if no pair of nodes from any separate Must-Link groups in g.coh 
%   is connected by an edge from g.ml (pairwise must-link constraints) or 
%   if any branch from g.ml does not exist in g.adj
% ok = false otherwise
% 
% Author: Ilya Tyuryukanov
% Date: 4 Oct 2015
% Last revision: 22 January 2015 (adjusting it as a PFgraph method)

% Initialize variables
adj = g.adj;
bran_ml = g.ml;
num_bran_ml = size(bran_ml, 1);
coh_ml = g.coh;
num_coh_ml = max(coh_ml(2,:));  % number of must-link-groups

% Initialize output
ok = true;
bran_viols = false(num_bran_ml, 1);

% Check if all pairwise must-links represent existing edges
for i = 1:1:num_bran_ml
  if (adj( bran_ml(i,1), bran_ml(i,2) ) == 0) && (adj( bran_ml(i,2), bran_ml(i,1) ) == 0)
    ok = false;
    bran_viols(i) = true;
  end
end

if ( isempty(coh_ml) || isempty(bran_ml) )
    coh_viols = [];
    return;
end

% Check if bran_ml contradict with inter-group cannot-links between coh_ml node groups
bias = min(coh_ml(2,:)) - 1; % to ensure that counting starts from one
coh_ml(2,:) = coh_ml(2,:) - bias;
coh_ml_bins = histcounts(coh_ml(2,:), 'BinMethod', 'integers'); % number of must-link buses in each coh_ml-group  
grp_pairs = combnk(1:num_coh_ml, 2);
tmp = coh_ml_bins(grp_pairs);
tmp = tmp(:,1).*tmp(:,2);  % number of CLs between each pair of coh_ml-groups
num_coh_cl = sum( tmp, 1);  % total number of CLs between coh_ml-groups
coh_cl = NaN(num_coh_cl, 2);  % preallocate cannot-links from coh_ml-groups 
if nargout > 2, coh_viols = NaN(num_coh_cl, 2); end  % preallocate the initial data for the 3rd output

k = 1;
for i = 1:1:size(grp_pairs,1)
  i1 = grp_pairs(i,1);
  i2 = grp_pairs(i,2);
  idx1 = coh_ml(2,:) == i1;
  idx2 = coh_ml(2,:) == i2;
  vec1 = coh_ml(1, idx1);
  vec2 = coh_ml(1, idx2);
  coh_cl( k : k+tmp(i)-1, : ) = allcomb(vec1, vec2);
  if nargout > 2
    idx1 = find(idx1);
    idx2 = find(idx2);
    coh_viols( k : k+tmp(i)-1, : ) = allcomb(idx1, idx2);
  end
  k = k + tmp(i);
end

% Take the possible reversed order of coh_cl's into account
coh_cl = [coh_cl; [coh_cl(:,2), coh_cl(:,1)]];
if nargout > 2, coh_viols = [coh_viols; [coh_viols(:,2), coh_viols(:,1)]]; end

% Write "violated" must-link branches
% bran_viols = intersect(bran_ml, coh_cl, 'rows'); % more optimal if 2nd output is optional
bran_viols = ismember(bran_ml, coh_cl, 'rows') | bran_viols;
ok = ~any(bran_viols);

% Write indices (into g.coh) of "violated" coherent generators (this is  
%   complementary to bran_viols if all bran_mls represent existing edges)
if nargout > 2
  idx_cl_viols = ismember(coh_cl, bran_ml, 'rows');
  coh_viols = coh_viols(idx_cl_viols,:);
end


