function [T_out, minor_cc] = contig(g, T, varargin)
%
% Syntax: [T_out, minor_cc] = contig(g, T, varargin)
%
% Purpose: Given the original power network graph (as a PFgraph object g)
%   and the clustering solution (cut_in), remove any possible minor
%   connected components to make the number of partitions equal to the
%   number of connected components. Then the new cutset is returned.
% 
% Input:
% g: a PFgraph object (incidence and adjacency matrices + constraints)
% T: input partition
% varargin: key-value pairs for this function (a struct or an isomorphic 
%   cell array). The keys and their possible values are listed below.
% 
% SplitLargeMinor: a bool flag splitting (if true) or not splitting of 
%  "too" large minor connected components into smaller ones (via Fiedler
%  bisection). 
% Objective ({'minphi'}/{'minimb'}): selects one of predefined objectives 
%  to guide the assignment of small components.
% 
% NumberOfConnCmp: the required number of connected components.
% 
% Output:
% T_out: output partition.
% 
% MAIN REFERENCES
% The main idea is from the METIS 5.1.0 contig.c function
% 
% Author: Ilya Tyuryukanov
% Date of first version: 28 November 2016

adj = g.adj;
m   = size(adj,1);

% Find all connected components and their characteristics
bus_in = sparse(T, 1:m, true);
cut_in = g.cutset(bus_in);
[~, c, ixs, eXp] = cc_info(g, cut_in);

% Return if we already have the required number of connected components
nc = numel(ixs);        % number of the actual connected components
np = size(cut_in, 1);   % number of partitions/islands


% varargin processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('SplitLargeMinor', false, @islogical);
p.addParameter('NumberOfConnCmp', np, @(x) Utils.isint(x) && isscalar(x));
p.addParameter('Objective', 'minphi', @ischar);
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
   'spelling. The allowed keys are: SplitLargeMinor.'], mfilename);
split = varinp.SplitLargeMinor;
np = varinp.NumberOfConnCmp;
obj = varinp.Objective;
if strcmp(obj,'minimb')
  vw = g.vw;
end

[~, idx] = sort(cellfun(@numel, ixs), 'descend');
ixs = ixs(idx);
minor_cc = ixs(np+1:nc);
if np >= nc
  T_out = T;
  return;
end

% Find the largest connected components   
[~, idx] = sortrows(eXp', [3,1], 'descend');
idx = idx';

% Subdivide large minor connected components into smaller ones
if split
  card_crit = max(2, ceil(eXp(3,idx(np))*0.1));
  c_large = ( eXp(3,:)>card_crit ) & ( eXp(3,:)<eXp(3,idx(np)) );
  c_large = find(c_large);
  n_large = numel(c_large);
  while n_large>0
    for i = 1:1:n_large
      c_curr = c_large(i);
      busidx = c == c_curr;
      g_tmp = g.extract_subgraph(busidx);
      busind = g_tmp.fiedler_bisect();
      init_idx = find(g_tmp.merge_map);
      c(init_idx(busind(2,:))) = nc+i;
    end
    bus_cc = sparse(c, 1:m, true);
    cut_cc = cutset(g, bus_cc);
    cut_cc(cut_cc~=1) = 0;
    [~, c, ~, eXp] = cc_info(g, cut_cc);  % update connected components
    nc = length(eXp(1,:));
    [~, idx] = sort(eXp(3,:), 'descend');
    c_large = ( eXp(3,:)>card_crit ) & ( eXp(3,:)<eXp(3,idx(np)) ) & ( eXp(3,:)>1 );
    c_large = find(c_large);
    n_large = numel(c_large);
  end
end

bus_cc = sparse(c, 1:m, true);
cut_cc = cutset(g, bus_cc);
cut_cc(cut_cc~=1) = 0;
c_minor = idx(np+1:end);
n_minor = numel(c_minor);

% Define the connection order (connect in order of max total weight to one
%  of the major connected components)
c_majr = idx(1:np);
pcut_l2m = zeros(n_minor, np);
for i = 1:1:n_minor
  c_curr = c_minor(i);
  pflow_cc = wgt_small_major(g, cut_cc, c, c_majr, c_curr);
  pcut_l2m(i, :) = pflow_cc;
end
maxcut_l2m = max(pcut_l2m, [], 2);   % the goal is max reduction of cut btw major components
[~, idx] = sort(maxcut_l2m, 'descend');

todo = c_minor(idx);  % connect well interconnected minor components first
pcut_l2m = pcut_l2m(idx,:);  % update the order of minor to major cuts as well
n_todo = n_minor;
itr_lft = n_minor+1;  % iteration limit of the while loop
% eXp_sep = eXp(1,:);
pcut_i = pcut_l2m(1,:);
while n_todo>0
  
  c_curr = todo(1);
  c_cand = find(pcut_i);  % pcut_i is cut to each major cc, the indices of which are stored in c_majr
  n_cand = numel(c_cand);
  if n_cand == 1
    cut_i = merge_cut(cut_cc, c_curr, c_majr(c_cand(1)));
  else
    fitnes1 = zeros(1,n_cand);
    fitnes2 = zeros(1,n_cand);
    for j = 1:1:n_cand
      cut_tmp = merge_cut(cut_cc, c_curr, c_majr(c_cand(j)));
      [~, ~, ~, eXp_tmp] = cc_info(g, cut_tmp);
      % eXp_dif = setdiff(eXp_tmp, eXp_sep);  % input order is important to extract the newly formed expansion
      [~, idx] = sort(eXp_tmp(3,:), 'descend');
      eXp_tmp = eXp_tmp(:, idx);
      eXp_tmp = eXp_tmp(1, 1:np);
      fitnes1(j) = max(eXp_tmp);
      fitnes2(j) = sum(eXp_tmp)/np;  % the updated NCut value
    end
    [~, c_join] = min(fitnes1);
    if nnz(fitnes1==fitnes1(c_join))>1
      round2 = find(fitnes1==fitnes1(c_join));
      [~, idx] = min(fitnes2(round2));
      c_join = round2(idx);
    end
    cut_i = merge_cut(cut_cc, c_curr, c_majr(c_cand(c_join)));
  end
  [~, c, ~, eXp] = cc_info(g, cut_i);
    
  bus_sep = sparse(c, 1:m, true);
  cut_cc = cutset(g, bus_sep);
  cut_cc(cut_cc~=1) = 0;
  [~, idx] = sort(eXp(3,:), 'descend');
  if numel(idx)>np
    c_majr = idx(1:np);
    c_todo = idx(np+1:end);
    n_todo = numel(c_todo);
    pcut_l2m = zeros(n_todo, np);
    for i = 1:1:n_todo
      c_curr = c_todo(i);
      pflow_cc = wgt_small_major(g, cut_cc, c, c_majr, c_curr);
      pcut_l2m(i, :) = pflow_cc;
    end
    maxcut_l2m = max(pcut_l2m, [], 2);  % the goal is max reduction of cut btw major components
    [~, idx] = sort(maxcut_l2m, 'descend');    
    todo = c_todo(idx);  % connect well interconnected minor components first
    pcut_l2m = pcut_l2m(idx,:);
    % eXp_sep = eXp(1,:);
    pcut_i = pcut_l2m(1,:);
  else
    n_todo = 0;
  end
  
  itr_lft=itr_lft-1;
  if itr_lft <= 0
    warning([ mfilename(), ':IterationLimitReached'],...
      ['[%s] The maximal number of iterations has been reached, and ',...
      'possibly not all minor connected components have been eliminated. ',...
      'Returning the partitionig reached so far.'], mfilename);
    break;
  end
end

T_out = c;
end

function cut_out = merge_cut(cut_cc, c1, c2)  
%
% Find the graph cut (e.g. in the original graph) after merging two
% components with indices (cut numbers) c1 and c2.
%
% cut_cc represents the cut index vectors of every connected component
% c1 repesents the cut number of the 1st connected component to merge
% c2 repesents the cut number of the 2nd connected component to merge
%
% cut_out represents the modified total cut index vector

cut1 = cut_cc(c1,:);
cut2 = cut_cc(c2,:);
common_cut = cut1 & cut2;

cut_out = cut_cc;
cut_out(c1,:) = xor(cut_cc(c1,:), common_cut);
cut_out(c2,:) = xor(cut_cc(c2,:), common_cut);
cut_out = any(cut_out,1);
end

function pflow_cc = wgt_small_major(g, cut_cc, c, c_majr, c_curr)
%
% Find connection weight between the minor connected component with index
% c_curr and the major connected components with indices 1 to n_large. The
% cut_cc matrix contains in its rows the cut index vectors of every
% connected component (in the order corresponding to the above mentioned
% indices). The original connected graph g is used to fetch the original
% weight data for requested edges.

cut_i = logical(cut_cc(c_curr,:));
[vx_pairs, pflow] = g.edges2adj(cut_i, true);
cc1 = c(vx_pairs(:,1));
cc2 = c(vx_pairs(:,2));
cc1(cc1==c_curr) = 0;
cc2(cc2==c_curr) = 0;
cc = cc1 + cc2;
cc = cc(:);
pflow_cc = accumarray(cc, pflow);
i_adj_cc = numel(pflow_cc);
i_majr = max(c_majr);
if i_adj_cc > i_majr
  pflow_cc = pflow_cc(1:i_majr);
else
  pflow_cc = [pflow_cc; zeros(i_majr-i_adj_cc,1)]; 
end
pflow_cc = pflow_cc(c_majr);
end


