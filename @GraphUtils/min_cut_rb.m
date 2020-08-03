function T = min_cut_rb(adj, vs, vw, debug)
%
% Given a reduced graph with its k nodes representing merged cluster cores
% (i.e., the "cluster core nodes")...
% Choose randomly a cluster core node, and compute k-1 mincuts between it
% and the remaining cluster core nodes. Select the smallest mincut or Ncut
% depending on the input. Repeat this process on the cluster core nodes 
% belonging to the resulting two parts again and again, until all cluster 
% core nodes are separated from each other. 
% This process resembles Gomory-Hu mincut trees, but it is computed  
% multiple times and the best Ncut may be chosen at each iteration as well.
%

old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party'),...
  fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party', 'matlab_bgl'));

switch nargin
  case 2
    debug = false;
    ncut = false;
  case 3
    if islogical(vw)
      debug = vw;
      ncut = false;
    else
      debug = false;
      ncut = true;
    end
  case 4
    ncut = true;
end

if (~isequal(adj, adj'))
  error('min_cut_rb:invalidParameter',...
    'the matrix must be symmetric.');
end

if (min(min(adj)) < 0)
  error('min_cut_rb:invalidParameter',...
    'the matrix cannot contain negative weights.');
end

m0 = size(adj,1);
k = numel(vs);

if k==2
  [adj_int] = round_adj(adj);
  [flow, ci] = max_flow(adj_int, vs(1), vs(2));
  T = zeros(m0, 1);
  T(ci==1) = 1;
  T(ci==-1) = 2;
  return;
end

% Create graph with constraints
adj0 = adj;
gen = false(1,m0);
gen(end-k+1:end) = true;
g = PFgraph('adj',adj,'gen',gen);
map = 1:m0;
g.coh = [map(end-k+1:end); 1:k];

% Create tree for the rb
g_curr = g;
bipart = tree( g_curr );  % graph bipartitioning tree
leaves = bipart.findleaves();
num_leaf = numel(leaves);
while num_leaf<k
  for j = 1:1:num_leaf
    g_curr = bipart.Node{leaves(j)};
    if isa(g_curr,'PFgraph')
      p_curr = leaves(j);  % current tree parent node index
      break;
    end
  end
  n_core = size(g_curr.coh, 2);
  adj = g_curr.adj;
  m = size(adj,1);
  vs = m-n_core+1:1:m;
  [adj_int, c_mult] = round_adj(adj); 
  flows = Inf(n_core-1, 1);
  part_inds = false(m, n_core-1);
  for j = 2:1:n_core
    [flow, ci] = max_flow(adj_int, vs(1), vs(j));
    part_ind = ci==sign(vs(1));
    part_inds(:,j-1) = part_ind;
    flows(j-1) = flow/c_mult;
    if ncut
      flows(j-1) = flows(j-1)*(1/vw(vs(1)) + 1/vw(vs(j)));
    end
    if debug
      assert( flow == sum(sum(floor(adj_int(part_ind,~part_ind)))) );
      cut_curr = sum(sum( adj(part_ind,:) )) - sum(sum( adj(part_ind, part_ind) ));
      assert( (round(cut_curr)==0 && round(flow/c_mult)==0) || abs(round(flow/c_mult) - round(cut_curr))/round(cut_curr)<1e-5);
    end
  end
  [minflow, idx_min] = min(flows);
  part_ind = part_inds(:,idx_min);
  % [mv, mw] = stoer(adj, vs);
 
  g1 = g_curr.extract_subgraph( part_ind); 
  g2 = g_curr.extract_subgraph(~part_ind);  
    
  % Add 2 child nodes to the currently bisected graph
  if size(g1.coh,2)>1
    [bipart] = bipart.addnode(p_curr, g1);
  else
    bipart = bipart.addnode(p_curr, part_ind);  
  end  
  if size(g2.coh,2)>1
    [bipart] = bipart.addnode(p_curr, g2);
  else
    part_ind = ~part_ind;
    bipart = bipart.addnode(p_curr, part_ind);     
  end
  
  leaves = bipart.findleaves();
  num_leaf = numel(leaves);
end

% Traverse the bipart tree starting from the leaf nodes upwards in order
% to recover the bus indicator vectors (busind) for the original graph
leaves = bipart.findleaves();
g_root = bipart.Node{1};
swap_map = g_root.merge_map;
busind = false(numel(leaves), size(g_root.adj,1));
for i = 1:1:numel(leaves)
  leaf = leaves(i);
  p_curr = bipart.Parent(leaf);
  currind = bipart.Node{leaf};
  while p_curr~=1
    idx = find(currind);
    g_curr = bipart.Node{p_curr};
    map = g_curr.merge_map;
    currind = ismembc(map, idx);
    p_curr = bipart.Parent(p_curr);
  end
  busind(i,:) = currind(swap_map);
end
T = GraphUtils.indvec2vxlbl(busind);
end


function [adj, c_mult] = round_adj(adj)
m = size(adj, 1);
[i, j, w0] = find(adj);
max_int = 2^28;
c_mult = max_int/max(w0);
w = w0*c_mult;
w(w<1.00000001) = 1;
adj = sparse(i, j, w, m, m);
end


% % function [adj, int_infinity] = round_adj(adj, vs)
% % m = size(adj, 1);
% % [i, j, w0] = find(adj);
% % int_infinity = 2^33;
% % mf = 1e6/max(w0);
% % max_int = 2^32;
% % while int_infinity>=max_int
% %   mf = mf/2;
% %   w = w0*mf;
% %   w(w<1.00000001) = 1;
% %   adj = sparse(i, j, w, m, m);
% %   int_infinity = sum(sum(adj))+2*sum(sum(adj(vs,:)))+1;  % this should be larger than any conceivable flow...
% % end
% % end