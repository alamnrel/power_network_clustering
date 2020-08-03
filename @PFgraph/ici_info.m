function [eXp, pcut, ixs, adj_sep, num_viols, viols, tot_outl_gen] = ici_info( g, cutind )
%
% Syntax: [eXp, pcut, ixs, adj_sep, num_viols, viols, tot_outl_gen] = ici_info( g, cutind )
%
% Purpose: Given the original power network graph (as a PFgraph object g)
%   and the islanding solution (cutind), find the quality metrics of this
%   partition (expansion ratios of conn. components, total power flow
%   disruption). The overall cutset is returned as well.
%
% Input:
%   g: a PFgraph object (incidence and adjacency matrices + constraints)
%   cutind: ith row is an index vector into the network branch list (or 
%     into the incidence matrix), for showing the branches that should be 
%     cut to disconnect the ith identified island.
%
% Output:
%   eXp(1,:): expansion ratios { sum(Pcut,k)/sum(Pij), i,j are the buses
%     in the kth island, k: (1st island) .. (last island) } of all islands
%   eXp(2,:): expansion ratios in terms of RatioCut { sum(Pcut,k)/|Vk|,
%     Vk are the buses in the kth island, k: (1st island)..(last island) }
%     of all islands
%   eXp(3,:): number of buses in each island
%   pcut: total power flow disruption in units of g.adj (usually in MW)
%   ixs: cell array of lists of vertex indices for each actual connected
%     component (after the removal of separator edges).
%   adj_sep: disconnected adjacency matrix
%   num_viols: number of constraints violations (NaN if constraints were 
%     not specified in the input graph object)
%   viols: cell array of strings describing each constraint violation
%   tot_outl_gen: number of coherent generators that do not belong to their
%     dominating connected component
% 
% A note: indicator vector into the network branch list (or into the incid.
%   matrix) showing the TOTAL islanding cutset (all lines to be tripped)
%   can be trivially obtained by:
%      cutind(cutind==2) = 0;
%      cut = any(cutind, 1);
%   As of September 2016 this 3rd output entry is replaced by the lists of
%   vertex indices for each actual connected component found.
%
% Author: Ilya Tyuryukanov
% Date of first version: 12 June 2015
% 03 February 2016: rewritten as a PFgraph method
% 22 September 2016: the 3rd output is now a cell array of node indices for
%   each cluster
% 27 September 2016: separating the calculation of the first four arguments
%   into a separate function, as this calculation is valuable on its own


% Return 1st four arguments
[adj_sep, C, ixs, eXp, pcut] = cc_info(g, cutind);
pcut = sum(pcut)/2;  % while summing individual pcuts of each island, each edge is counted twice

% Issue a warning in case something is suddenly wrong
S = numel(ixs);  % number of the actual connected components
if S ~= size(cutind, 1)
  warning([ mfilename, ':BadResults'],...
    ['[%s] The actual number of connected components is %d, but initial ',...
    'number of areas is %d'], mfilename, S, size(cutind,1));
end

% Analyze constraints (if requested)
if nargout > 4
  bran_ml = g.ml;
  coh_ml = g.coh;
  busmap = g.bus;
  num_coh_ml = max(coh_ml(2,:));  % number of must-link-groups
  num_viols = 0;
  viols = cell(0,0);
  
  % Check pairwise branch must-links
  for i = 1:1:size(bran_ml,1)
    if (adj_sep( bran_ml(i,1), bran_ml(i,2) ) == 0) && (adj_sep( bran_ml(i,2), bran_ml(i,1) ) == 0)
      num_viols = num_viols + 1;
      msg = sprintf(['Must-link branch between buses %d - %d (node indices ',...
        '%d - %d) was cut'], busmap(bran_ml(i,1)), busmap(bran_ml(i,2)),...
        bran_ml(i,1), bran_ml(i,2));
      viols(num_viols) = {msg};
    end
  end
  
  % Check group must-links (using the connected components found above)
  cc_of_grp = cell(num_coh_ml, 1);
  grp_of_cc = cell(S, 1);
  gen_assgn = zeros(S,num_coh_ml);
  for j = 1:1:num_coh_ml
    idx = coh_ml(2,:) == j;
    grp = coh_ml(1,idx);
    grp_cc = C(grp);  % indices of connected components for the ith must-link group
    grp_ct = tabulate(grp_cc);
    cc = grp_ct(:,1);
    ct = grp_ct(:,2);
    gen_assgn(cc,j) = ct;
    grp_cc = unique(grp_cc);
    cc_of_grp(j) = {grp_cc};
    for k = 1:1:numel(grp_cc)
      grp_of_cc{grp_cc(k)} = [grp_of_cc{grp_cc(k)}, j];
    end
  end
  tot_outl_gen = 0;    
  if S>=num_coh_ml
    num_col = num_coh_ml;   % if enough components, test group-wise
  else  
    gen_assgn = gen_assgn'; % if not enough components, test component-wise
    num_col = S;
  end
  for j = 1:1:num_col
    tot_outl_gen = tot_outl_gen + sum(gen_assgn(:,j)) - max(gen_assgn(:,j));
  end  
  
  % Each must-link nodes group should be within a single connected component
  n_cc_in_grp = cellfun(@numel, cc_of_grp, 'UniformOutput', true);
  viols_cc_in_grp = n_cc_in_grp > 1;
  disconn_grps = find(viols_cc_in_grp);
  n_cc_in_viol_grp = n_cc_in_grp(viols_cc_in_grp);
  for i = 1:1:nnz(viols_cc_in_grp)
    num_viols = num_viols + 1;
    msg = sprintf('Nr. %d must-link nodes group is spread among %d (out of %d) connected components', disconn_grps(i), n_cc_in_viol_grp(i), S);
    viols(num_viols) = {msg};
  end
  % Each connected component should contain only one must-link nodes group
  n_grp_in_cc = cellfun(@numel, grp_of_cc, 'UniformOutput', true);
  viols_grp_in_cc = n_grp_in_cc > 1;
  incoherent_cc = find(viols_grp_in_cc);
  n_grp_in_viol_cc = n_grp_in_cc(viols_grp_in_cc);
  for i = 1:1:nnz(viols_grp_in_cc)
    num_viols = num_viols + 1;
    msg = sprintf('Nr. %d connected component contains nodes from %d (out of %d) must-link groups', incoherent_cc(i), n_grp_in_viol_cc(i), num_coh_ml);
    viols(num_viols) = {msg};
  end
end