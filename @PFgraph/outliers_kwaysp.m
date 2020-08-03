function [g_pp, merges] = outliers_kwaysp(g, varargin)
% 
% Find and merge graph outliers (both nodes and clusters) satisfying the 
% specified input parameters.
% 
% Reference: I. Tyuryukanov, M. Popov, M. A. M. M. van der Meijden, and V. 
% Terzija. "Spectral MST-based Graph Outlier Detection with Application to 
% Clustering of Power Networks". In: Proc. 20th Power Systems Computation 
% Conference (PSCC 2018). Dublin, Ireland, 2018, pp. 1-8.
% 

%% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addRequired('nv', @(x) isscalar(x) && Utils.isint(x) && x>=1); 
p.addRequired('max_card',  @(x) isscalar(x) && Utils.isint(x) && x>=1);
p.addRequired('score_max', @(x) isscalar(x) && isnumeric(x) && Utils.isdbl(x) && x>0);
p.addParameter('debug', 0, @(x) isscalar(x) && isnumeric(x) && x>=0);
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check ',...
  'your spelling. The allowed keys are: debug.'], mfilename);
nv = varinp.nv;
max_card = varinp.max_card;
score_max = varinp.score_max;
debug = varinp.debug;

% Higher dbg_lvl corresponds to longer checks
dbg_lvl_1 = 1;
dbg_lvl_2 = 2;
dbg_lvl_3 = 3;

adj0 = g.adj;
inc = g.inc;
m = size(adj0,1);         % number of buses
n = size(inc, 2);         % number of branches
vd = full(sum(adj0,2));   % weighed degrees of each vertex for VOLUME calculation
epsilon = eps('double');

%% Construct SPEC-MST
L_sym = GraphUtils.createLaplacian(adj0, 'sym');
L_sym(1:m+1:m*m) = 0; L_sym = -L_sym;  % Anorm for > speed and accuracy
L_rwk = GraphUtils.createLaplacian(adj0, 'rwk');
L_rwk(1:m+1:m*m) = 0; L_rwk = -L_rwk;  % for single-node outlier scores
% [Y, lam] = GraphUtils.spClust(L_sym, 'la', nv);  % LARGEST ALGEBRAIC
% Y = bsxfun(@times, 1./sqrt(vd), real(Y(:,2:nv)));
% adjD = GraphUtils.buildDistanceGraph(Y, adj0, 'rwk');
[Y, lam] = GraphUtils.spClust(L_sym, 'la', nv);  % <- (instead of the above three lines)
Y = Utils.normalize_rows(Y);
adjD = GraphUtils.buildDistanceGraph(Y, adj0, 'sym');
root = 1; % set root to node 1 for a robust implementation
%[tree, pred] = graphminspantree(adjD, 'Method', 'Prim');   % (Prim) - MATLAB mst is very slow
mst_sp = kruskal_mst(adjD);
[~, ~, pred, children] = mst_bfs(mst_sp,root);
pred(root) = 0;   % to make compatible with bfs_mex
[f, t, wx] = find(triu(mst_sp));
[wx, idx] = sort(wx, 'ascend');   % 'ascend' is for the KRUSKAL stage below, as we decide  to examine ALL fundamental cuts anyway
f = f(idx); t = t(idx);
wx(wx<3*epsilon) = 3*epsilon;   % wx appears sorted!
mst_sp_bin = sparse([f;t], [t;f], true(2*m-2,1), m, m);
mst_sp = sparse([f;t], [t;f], [wx;wx], m, m);

% Retrieve the original MST-weights
edg_idx = sub2ind([m, m], f, t);
w0 = adj0(edg_idx);
w0 = full(w0); 

if debug>=dbg_lvl_2
  stack = zeros(m,1);
  assert(pred(root)==0);
  stack(1) = root;  % the root node enters first
  curr_node = 2;
  children_tst = cell(m,1);
  pred = pred(:);
  for i = 1:1:m
    childs = find(pred==stack(i));
    children_tst{stack(i)} = childs(:);
    stack(curr_node : curr_node+numel(childs)-1) = childs;
    curr_node = curr_node+numel(childs);
  end
  stack = flipud(stack);
  children_tst = cellfun(@sort, children_tst, 'unif', false);
  children = cellfun(@sort, children, 'unif', false);
  assert(isequal(children_tst,children));
  
  desc_count_tst = ones(m,1);
  i = 1;
  while pred(stack(i))~=0
    desc_count_tst( pred(stack(i)) ) = desc_count_tst( pred(stack(i)) ) + desc_count_tst( stack(i) );
    i = i+1;
  end
end

% Find the map of parents to total children
descendants = num2cell((1:m)');
descendants = getdescendants(root, descendants, children);
desc_count = cellfun(@numel, descendants);
if debug>=dbg_lvl_2
  assert(isequal(desc_count,desc_count_tst));
end

%% Examine the fundamental cuts
outl_fund = cell(m,1);
outl_sngl = cell(m,1);
score_sngl = Inf(m,1);
score_fund = Inf(m,1);
nodes_all = (1:m)';
k_ou = 1;
k_ou_sngl = 1;
num_tst = 1;   % purely for UNCONDITIONAL merging of "suspicious singletons"   // m-ceil(0.1*m)
k_start = min(num_tst, m-ceil(0.9*m));
for k = k_start:1:m-1
  
  if t(k) == pred(f(k))
    node_down = f(k);
    edge_curr = [f(k),t(k)];
  elseif f(k) == pred(t(k))
    node_down = t(k);
    edge_curr = [t(k),f(k)];
  else
    error('This branch is not justified!');
  end
  desc_down = desc_count(node_down);
  desc_rest = m - desc_down;
  [desc_curr, idx_bfs] = min([desc_down, desc_rest]);
  
  if debug>=dbg_lvl_3
    mst_tst = mst_sp_bin;
    mst_tst(f(k),t(k))=0; mst_tst(t(k),f(k))=0;   % #disconnect the tree
    [S, T] = GraphUtils.alecconncomp(mst_tst);
    assert(S == 2);
    bfs1 = find(bfs(mst_tst, f(k))~=-1);
    bfs2 = find(bfs(mst_tst, t(k))~=-1);
    
    num1 = nnz(T==1);
    num3 = numel(bfs1);
    num4 = numel(bfs2);
    
    if num1 == num3
      assert(isequal(find(T(:)==1), sort(bfs1)));
      assert(isequal(find(T(:)==2), sort(bfs2)));
    elseif num1 == num4
      assert(isequal(find(T(:)==1), sort(bfs2)));
      assert(isequal(find(T(:)==2), sort(bfs1)));
    else
      error(';;;');
    end
    
    bfs_start = edge_curr(idx_bfs);
    bfs_dist = bfs_mex(mst_tst, bfs_start, 0);
    card_merg_bfs = find(bfs_dist~=-1);
  end
  
  if desc_curr<=max_card
    if desc_curr==desc_down
      card_merg = descendants{node_down};
    else
      card_merg = my_setdiff(nodes_all, descendants{node_down});
    end
    if debug>=dbg_lvl_3
      assert( isequal(sort(card_merg_bfs), sort(card_merg)) );
      assert(numel(card_merg)==desc_curr);
    end
    if (desc_curr==1 && k>=num_tst)
      outl_sngl{k_ou_sngl} = card_merg;
      %score_sngl(k_ou_sngl) = mst_sp(f(k),t(k));
      score_sngl(k_ou_sngl) = max(L_rwk(:,card_merg));
      k_ou_sngl = k_ou_sngl + 1;
      continue;
    end    
    vol = sum(vd(card_merg));
    vol_tst = 0;
    
    phi_pre = w0(k)/vol;   % lower estimate of expansion
    if phi_pre<score_max
      links = sum(adj0(:,card_merg), 2);
      links = sum(links(card_merg));
      %links = sum(sum( adj0(card_merg,card_merg) ));
      phi = 1-links/vol;
      if phi<score_max
        score_fund(k_ou) = phi;
        outl_fund{k_ou} = card_merg(:);
        k_ou = k_ou + 1;
      end
    end
  else
    continue;
  end
end
outl_fund(k_ou:end) = [];
score_fund(k_ou:end) = [];
outl_sngl(k_ou_sngl:end) = [];
score_sngl(k_ou_sngl:end) = [];

%% Merge MST successfully
if debug>=dbg_lvl_1
  [wx, idx] = sort(wx, 'ascend');
  assert( isequal( (idx(:)'), 1:1:m-1 ) );
end

k_ou = 1;
outl_krus = cell(m,1);
score_krus = Inf(m,1);
[prev, siz, memb] = makeset(m);   % start with every vertex in its connected component
for k = 1:1:m-1
  i = f(k);
  j = t(k);
  
  %Inlined ROOT()
  while i ~= prev(i)
    i = prev(i);
  end
  r_i = i;
  while j ~= prev(j)
    j = prev(j);
  end
  r_j = j;
  assert(r_i~=r_j);  %#DEBUG-related!
  
  if (siz(r_i)+siz(r_j)<=max_card)
    
    % Inlined UNION()
    if siz(r_i)<siz(r_j)
      prev(r_i) = r_j;
      memb{r_j} = [memb{r_j};memb{r_i}];
      memb{r_i} = [];
      siz(r_j) = siz(r_j) + siz(r_i);
      r_out = r_j;
    else
      prev(r_j) = r_i;
      memb{r_i} = [memb{r_i};memb{r_j}];
      memb{r_j} = [];
      siz(r_i) = siz(r_i) + siz(r_j);
      r_out = r_i;
    end
    
    % Calculate scores and possibly save
    card_merg = memb{r_out};
    vol = sum(vd(card_merg));
    links = sum(adj0(:,card_merg), 2);
    links = sum(links(card_merg));
    % links = sum(sum( adj0(card_merg,card_merg) ));
    phi = 1-links/vol;
    if phi<score_max
      score_krus(k_ou) = phi;
      outl_krus{k_ou} = card_merg;
      k_ou = k_ou + 1;
    end
  end
end
outl_krus(k_ou:end) = [];
score_krus(k_ou:end) = [];


%% Combine all methods
outl = [outl_fund; outl_krus];
ou_score = [score_fund; score_krus];

[score_sngl, idx] = sort(score_sngl, 'ascend');
outl_sngl = outl_sngl(idx);
outl_sngl(ceil(0.01*m)+1:end) = [];
score_sngl(ceil(0.01*m)+1:end) = [];

if isempty(outl) && isempty(outl_sngl)
  g_pp = g;
  merges = {1:1:m};
  return;
end

%% Process outl_mult to remove clusters with intersecting node sets
%outl = cellfun(@sort, outl, 'unif', false);   % convenient, but inefficient if many "outl"
%[outl, ou_uniq] = uniquecell(outl);
%ou_score = ou_score(ou_uniq);
[ou_score, idx] = sort(ou_score, 'ascend');
outl = outl(idx);   % sort so that more severe come first!

num_ou = numel(outl);
if num_ou>0
  node_tabl = false(num_ou,m);
  for i = 1:1:num_ou
    node_tabl(i,outl{i}) = true;
  end
  idx_del = false(num_ou,1);   % clusters staged for deletion (due to obvious overlaps)
  idx_exl = false(num_ou,1);   % processed clusters (either to be deleted, or already processed)
  card_ou = sum(node_tabl,2); % number of entries in each "outlier" 
  curr = 1;
  while any(curr)  
    overlap  = node_tabl(curr,:);    
    num_curr = nnz(overlap);    
    num_same = card_ou == num_curr;
    incl_all = all( node_tabl(:,overlap), 2);
    incl_crs = any( node_tabl(:,overlap), 2);
    curr_del = (incl_crs & ~incl_all) | (incl_all & num_same);  
    row_excl = curr_del;
    node_tabl(row_excl,:) = false;  % remove intersecting clusters from further consideration    
    card_ou(row_excl,:) = 0;
    idx_exl = idx_exl | row_excl;
    curr_del(curr) = false;  % but don't stage the current ones for deletion!
    idx_del = idx_del | curr_del;
    curr = find(~idx_exl,1);
  end
  outl(idx_del) = [];
  ou_score(idx_del) = [];
end
%}

%% Re-unite
if debug>=dbg_lvl_1
  [ou_score, idx] = sort(ou_score, 'ascend');   % Re-sort (for fun..)
  assert(isequal(idx(:)', 1:1:numel(idx)));
end
outl = [outl;outl_sngl];

%% Sum up all valid clusters
outl = process_clust(outl, adj0, inc);   % always merge starting intersections!
if debug>=dbg_lvl_2
  outlY = process_clustY(outl);
  outl = process_clust(outl, adj0, inc);
  outlY = cellfun(@sort, outlY, 'unif', false);
  outl = cellfun(@sort, outl, 'unif', false);
  [outlY, ou_uniq] = uniquecell(outlY);
  [outl, ou_uniq] = uniquecell(outl);
  assert(isequal(outlY,outl))
end

% This is for infos only:
n_ou = numel(outl);
vx_red = sort(vertcat(outl{:}), 'ascend');   % all "outlier nodes"
if debug>=dbg_lvl_2
  phi = zeros(n_ou,1);
  for k = 1:1:n_ou
    vol(k) = sum(vd(outl{k}));
    links(k) = sum(sum( adj0(outl{k},outl{k}) ));
    phi(k) = 1-links(k)/vol(k);
  end
  [phi,idx] = sort(phi, 'ascend');
  outl = outl(idx);
end

%% Augment the outliers
card_ou = cellfun(@numel, outl);
[~, idx] = sort(card_ou, 'ascend');
outl = outl(idx);
ixs_ext = cell(numel(outl), 1);
for i = 1:1:n_ou
  idx_out = outl{i};
  ml_cand = adj0(idx_out,:);
  ml_cand_1 = ml_cand;
  ml_cand_2 = ml_cand;
  ml_cand_1(:,vx_red) = 0;
  ml_cand_2(:,idx_out) = 0;
  
  % First try to avoid merging clusters in "outl"
  [~,k,w] = find(ml_cand_1);
  if ~isempty(w)
    [~, idx] = max(w);
    idx1 = k(idx);
    vx_red = [vx_red;idx1];     % update the list of ou-included-vetrices
  else
    [~,k,w] = find(ml_cand_2);
    [~, idx] = max(w);
    idx1 = k(idx);
  end
  ixs_ext{i} = [outl{i}; idx1];  % add the strongest link as must-link
end

% Merge clusters if they share common vertices
ixs_ext = process_clust(ixs_ext, adj0, inc);

%% Wrap it up (with the new weights)
if nargout>0
  g_xx = g;  
  g_xx.adj = adj0;
  g_pp = g_xx.merge_nodes(ixs_ext);
  merges{1} = g_pp.merge_map;
  clear('g_xx');
end
end

function descendants = getdescendants(node, descendants, children)
if ~any(children{node})
  % descendants{node} = node;
  return;
else
  for k = children{node}'
    descendants = getdescendants(k, descendants, children);
    descendants{node} = [descendants{node}; descendants{k}];
  end
end
end

function [prev, siz, memb] = makeset(N)
prev = (1:N)';
siz  = ones(N,1);
memb = num2cell((1:N)');
end

function o = ROOT(i,prev)
while i ~= prev(i)
  i = prev(i);
end
o = i;
end

function [prev, siz, memb, r_out] = UNION(i, j, prev, siz, memb)
r_i = ROOT(i,prev);
r_j = ROOT(j,prev);
if r_i == r_j, return; end

if siz(r_i)<siz(r_j)
  prev(r_i) = r_j;
  memb{r_j} = [memb{r_j};memb{r_i}];
  memb{r_i} = [];
  siz(r_j) = siz(r_j) + siz(r_i);
  r_out = r_j;
else
  prev(r_j) = r_i;
  memb{r_i} = [memb{r_i};memb{r_j}];
  memb{r_j} = [];
  siz(r_i) = siz(r_i) + siz(r_j);
  r_out = r_i;
end
end

function ixs_out = process_clust(clust, adj, inc)
m = size(adj, 1);
n = size(inc, 2); % number of branches

links = false(1,n);
for i = 1:1:numel(clust)
  vx_curr = clust{i};
  links_curr = sum(inc(vx_curr,:), 1)==2;
  links = links | links_curr;
end

[idx_ou] = GraphUtils.edges2adj(adj, inc, links, false);
clear('adj', 'inc');
f = idx_ou(:,1); t = idx_ou(:,2); w = ones(size(idx_ou,1), 1);

adj_sep = sparse([f;t], [t;f], [w;w], m, m );
[S, T] = GraphUtils.alecconncomp(adj_sep);

vx_all = unique(vertcat(clust{:}));
vx_rst = my_setdiff(1:m, vx_all);
T(vx_rst) = 0;   % these component indicators are irrelevant
comp = nonzeros( unique(T) );
comp_num = numel(comp);
ixs_out = cell(comp_num,1);
for i = 1:1:comp_num
  ixs_curr = find(T==comp(i));
  ixs_out{i} = ixs_curr;
end
ixs_out = cellfun(@(x) x(:), ixs_out, 'unif', false);
end


function ixs_out = process_clustY(clust)  % -> slow equivalent of process_clust(), for debug purpose
n_vx = cellfun(@numel, clust);
[~, idx] = sort(n_vx, 'descend');
clust = clust(idx);
clust = cellfun(@(x)sort(x,'ascend'), clust, 'unif', false);  % for ismembc

n_del  = 0;
n_done = 0;
n_orig = numel(clust);
j = 1;
while (n_del + n_done) < n_orig
  memberships = cellfun(@(x) ismembc(x,clust{j}), clust, 'unif', false);
  overlap = my_setdiff(find(cellfun(@any, memberships)), j);
  last = numel(overlap);
  while last>0
    appnd = my_setdiff(clust{overlap(1)}, clust{j});
    appnd = appnd(:);
    clust_merge = [clust{j}; appnd];
    %if numel(clust_merge)<=max_card+2 % +2 is to account for singleton vertices that always should be merged
    clust{j}(1:numel(clust_merge)) = sort(clust_merge, 'ascend');
    n_del = n_del + 1;
    clust(overlap(1)) = [];
    memberships = cellfun(@(x) ismembc(x,clust{j}), clust, 'unif', false);
    overlap = my_setdiff(find(cellfun(@any, memberships)), j);
    last = numel(overlap);
    %end
  end
  j = j+1;
  n_done = n_done + 1;
end
ixs_out = clust;
end

function eXp = calc_eXp(adj, clust)  % #debug
links = sum(sum( adj(clust,clust) ));
vol = sum(sum( adj(:,clust) ));
eXp = 1 - links/vol;
end








