function trees = hclust_constrained_path(g)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% Purpose:                                                                %
% Perform hierarchical clustering, starting with the generators as        %
% clusters and trying to merge these clusters if possible by using        %
% shortest paths.                                                         %
%                                                                         %
% Input:                                                                  %
% g: graph object with adjacency distance matrix                          %
%                                                                         %
% Output:                                                                 %
% trees: groups of nodes in each group tree                               %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

debug = false;
javaaddpath(fileparts(mfilename('fullpath')));

old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath( fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party','matlab_bgl'),...
    fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party') );

adj = g.adj;
coh = g.coh;
m = size(adj,2);
epsilon = eps('double');
gen_all = coh(1,:);
grp_gen = coh(2,:);
[grp_gen,idx] = sort(grp_gen,'ascend');
gen_all = gen_all(idx);
grp = unique(grp_gen);
num_grp = numel(grp);
num_gen = numel(gen_all);

ixs = cell(num_grp,1);    % summary of groups' target nodes
pred = cell(num_grp,1);   % hierarchical component tree of each group
memb = cell(num_grp,1);   % members of each component of each subtree
siz = cell(num_grp,1);    % size of each component in each subtree (needed for UNION-FIND)
for g_cur = 1:1:num_grp
  ixs_horz = gen_all(grp_gen==grp(g_cur));
  ixs{g_cur} = ixs_horz(:);
  pred{g_cur} = 1:1:numel(ixs{g_cur});  % components labels of each subtree start from one
  memb{g_cur} = num2cell(ixs{g_cur});  % each component is on its own
  siz{g_cur} = ones(numel(ixs{g_cur}),1);
end
grp_siz = cellfun(@numel, ixs);

% Initialize distance and predecessor matrices
adj = adj';
djk_opt.istrans = 1;
djk_opt.nocheck = 1;

% Initialize nodes and edges of each group and "taboo adjacency matrix"
lin_grp = cell(num_grp,1);
grp_all = 1:1:num_grp;
for g_cur = grp_all
  lin_idx = getadjedg(adj, ixs{g_cur});
  grp_rst = my_setdiff(grp_all,g_cur);
  for g_rst = grp_rst
    lin_grp{g_rst}=[lin_grp{g_rst};lin_idx];
  end
end

% Build the initial set of candidate paths-edges for merging
% 1: <path "from">; 2: <path "to">; 3: <path distance>; 4: <generator group>;
% 5: <group component of the from edge>; 6: <group component of the to edge>;
cmp_pth = zeros(num_grp,m);   % number of paths of each group going through a point
pth_ini = zeros(m,1);
siz_ini = round(m*m);
comp = PriorityComparator();
qcur = java.util.PriorityQueue(siz_ini,comp);

for g_cur = 1:1:num_grp
  memb_cur = memb{g_cur};
  cmp_idx = ~(cellfun(@isempty,memb_cur));  
  if sum(cmp_idx)==1, continue; end   % all subgroups are already in one group      
  adj_tmp = adj;
  adj_tmp( lin_grp{g_cur} ) = 0;
  cmp_idx = find(cmp_idx);
  for i = 1:1:numel(cmp_idx)-1
    srces = cmp_idx(i);
    sinks = cmp_idx(i+1:end);    
    n_i = numel(sinks);
    from = [memb_cur{srces}(1), srces];       
    v_dst = cellfun(@(x) x(1), memb_cur(sinks), 'unif', true);
    dest = zeros(n_i, 2); dest(:,1) = v_dst(:); dest(:,2) = sinks(:);    
    [new_new, pth_new, cmp_pth] = getnewpth(adj_tmp, g_cur, from, dest, cmp_pth, pth_ini, djk_opt);
    n_i = numel(pth_new);  
    for j = 1:1:n_i
      data = [new_new(j,:)'; pth_new{j}];
      add(qcur,data);
    end
 end   
end
%{
if any(sum(logical(cmp_pth), 1)>1)
  cprintf('magenta', 'Some of the initial paths contradict each other.\n');
end
%}

nlate = 0;
qlate = java.util.PriorityQueue(siz_ini,comp);
while num_gen>num_grp   % num_gen-num_grp merges are needed to merge all generators  
  if size(qcur)==0
    nlate = [nlate; size(qlate)];
    if nlate(end-1)==nlate(end) || nlate(end)==0   % if all prev. paths are done  
      pts_tbu = find(sum(logical(cmp_pth),1)>1);   % compute new tabu-points
      if numel(pts_tbu)==0
        break;
      else
        lin_idx = getadjedg(adj, pts_tbu);
      end
      for g_cur = 1:1:num_grp
        memb_cur = memb{g_cur};
        siz_cur = siz{g_cur};
        cmp_idx = ~(cellfun(@isempty, memb_cur));
        if sum(cmp_idx)==1, continue; end          % all subgroups are already in one group
        adj_tmp = adj;
        adj_tmp( lin_grp{g_cur} ) = 0;
        adj_tmp( lin_idx ) = 0;
        cmp_idx = find(cmp_idx);
        [~,srces] = max(siz_cur);
        sinks = cmp_idx(cmp_idx~=srces);
        n_i = numel(sinks);
        from = [memb_cur{srces}(1),srces];
        v_dst = cellfun(@(x) x(1), memb_cur(sinks), 'unif', true);
        dest = zeros(n_i, 2); dest(:,1) = v_dst(:); dest(:,2) = sinks(:);
        [new_new, pth_new, cmp_pth] = getnewpth(adj_tmp, g_cur, from, dest, cmp_pth, pth_ini, djk_opt);
        n_i = numel(pth_new);
        for j = 1:1:n_i
          data = [new_new(j,:)'; pth_new{j}];
          add(qcur,data);
        end                
      end
      if size(qcur)==0   % all paths are unreachable, nothing to add
        break;      
      end
    else
      qcur = qlate;
      qlate = java.util.PriorityQueue(siz_ini,comp);
    end
  end
  
  cur = poll(qcur);
  row_cur = cur(1:6);  % current components
  g_cur = row_cur(4); cf_cur = row_cur(5); ct_cur = row_cur(6);
  pts_pth = cur(7:end);  % path between current components      
  
  if debug
    T = zeros(1,m);
    for k = 1:1:numel(memb), tmp = vertcat(memb{k}{:}); T(tmp) = k; end  %part_viz(adj, T);    
  end
  
  % Extract tree components...
  tree_cur = pred{g_cur};
  while cf_cur ~= tree_cur(cf_cur)  % Inlined ROOT()
    cf_cur = tree_cur(cf_cur);
  end
  rf_cur = cf_cur;
  while ct_cur ~= tree_cur(ct_cur)
    ct_cur = tree_cur(ct_cur);
  end
  rt_cur = ct_cur;
  if rf_cur==rt_cur  % the path is superfluous
    cmp_pth(g_cur,pts_pth) = cmp_pth(g_cur,pts_pth)-1;
    continue;
  end
  
  % Decide if the components at hand are to be merged now
  if numel(pts_pth)==0
    goon = true;
  else
    idx_grp = any(cmp_pth(:,pts_pth),2);  % groups that use poits of this path
    if all(find(idx_grp)==g_cur)
      goon = true;
    else
      goon = false;
    end
  end
  if goon
    
    % Initialize the variables of g_cur'th subtree
    memb_cur = memb{g_cur};
    siz_cur = siz{g_cur};
    
    % Inlined UNION() [if the vertices are not already connected]
    if siz_cur(rf_cur) < siz_cur(rt_cur)
      tree_cur(rf_cur) = rt_cur;
      memb_cur{rt_cur} = [memb_cur{rt_cur};pts_pth;memb_cur{rf_cur}];
      memb_cur{rf_cur} = [];
      siz_cur(rt_cur) = siz_cur(rt_cur) + siz_cur(rf_cur);
      r_out = rt_cur;
    else
      tree_cur(rt_cur) = rf_cur;
      memb_cur{rf_cur} = [memb_cur{rf_cur};pts_pth;memb_cur{rt_cur}];
      memb_cur{rt_cur} = [];
      siz_cur(rf_cur) = siz_cur(rf_cur) + siz_cur(rt_cur);
      r_out = rf_cur;
    end
    
    % Add candidate paths-edges from the newly added points
    num_cmp = grp_siz(g_cur);
    cmp_idx = 1:1:num_cmp;   % all initial component indices of the current group    
    if numel(pts_pth)>0
      cmp_rst = find(cmp_idx~=r_out);
      new_pth = vertcat(memb_cur{cmp_rst});
      dest = zeros(numel(new_pth),2);   % dst vertices and their components
      dest(:,1) = new_pth;
      new_siz = cellfun(@numel,memb_cur(cmp_rst));
      eom = 1;
      for k = 1:1:numel(cmp_rst)
        dest(eom:eom+new_siz(k)-1,2) = cmp_rst(k);
        eom = eom+new_siz(k);
      end      
      adj_tmp = adj;
      adj_tmp( lin_grp{g_cur} ) = 0;
      for i = 1:1:numel(pts_pth)
        from = [pts_pth(i), r_out];  % from vertex and its component
        [new_new, pth_new, cmp_pth] = getnewpth(adj_tmp, g_cur, from, dest, cmp_pth, pth_ini, djk_opt);
        n_i = numel(pth_new);
        for j = 1:1:n_i
          data = [new_new(j,:)'; pth_new{j}];
          add(qcur,data);
        end  
      end
    end
            
    % Update group globals after updating group paths
    num_gen = num_gen-1;
    if numel(pts_pth)>0
      lin_idx = getadjedg(adj, pts_pth);
      grp_rst = my_setdiff(grp_all,g_cur);
      for g_rst = grp_rst      
        lin_grp{g_rst}=[lin_grp{g_rst};lin_idx];
      end
    end
    
    % Put all back
    pred{g_cur} = tree_cur;
    memb{g_cur} = memb_cur;
    siz{g_cur} = siz_cur;    
    
  else
    qlate.add(cur);
  end
end

% PRIORITY-MERGE THE PRESENT DISCONNECTED COMPONENTS
prpr = [];   % "priority pairs"
for g_cur = 1:1:num_grp
  memb_cur = memb{g_cur};
  siz_cur = siz{g_cur};
  cmp_idx = find(~(cellfun(@isempty, memb_cur)));
  if numel(cmp_idx)==1, continue; end   % all subgroups are already in one group
  [~,srces] = max(siz_cur);
  sinks = cmp_idx(cmp_idx~=srces);
  cmp_cur = [repmat(g_cur,numel(sinks),1),repmat(srces,numel(sinks),1),sinks,siz_cur(sinks)];
  prpr = [prpr; cmp_cur];
end
if ~isempty(prpr)
  prpr = sortrows(prpr,4,'descend');
end

for i = 1:1:size(prpr,1)
  g_cur = prpr(i,1);
  memb_cur = memb{g_cur};
  adj_tmp = adj;
  adj_tmp( lin_grp{g_cur} ) = 0;
  from = [memb_cur{prpr(i,2)}(1),prpr(i,2)];
  dest = [memb_cur{prpr(i,3)}(1),prpr(i,3)];
  [new_new, pth_new, cmp_pth] = getnewpth(adj_tmp, g_cur, from, dest, cmp_pth, pth_ini, djk_opt);
  n_i = numel(pth_new);
  if n_i==0
    continue;
  end
  
  tree_cur = pred{g_cur};
  siz_cur = siz{g_cur};
  pts_pth = pth_new{:};
  tree_cur(prpr(i,3)) = prpr(i,2);
  memb_cur{prpr(i,2)} = [memb_cur{prpr(i,2)};pts_pth;memb_cur{prpr(i,3)}];
  memb_cur{prpr(i,3)} = [];
  siz_cur(prpr(i,2)) = siz_cur(prpr(i,2)) + siz_cur(prpr(i,3));
  num_cur = numel(pts_pth);  assert(num_cur>0);
  lin_idx = getadjedg(adj, pts_pth);
  grp_rst = my_setdiff(grp_all,g_cur);
  for g_rst = grp_rst    
    lin_grp{g_rst}=[lin_grp{g_rst};lin_idx];
  end
  
  % Put all back
  pred{g_cur} = tree_cur;
  memb{g_cur} = memb_cur;
  siz{g_cur} = siz_cur;
end

% Merge groups and produce final clustering
for k = 1:1:num_grp
  idx_tree = ~(cellfun(@isempty, memb{k}));
  if sum(idx_tree)==1
    memb{k} = memb{k}{idx_tree};
  elseif sum(idx_tree)>1
    cmp_card = cellfun(@numel, memb{k});
    [~, idx_max] = max(cmp_card);
    memb{k} = memb{k}{idx_max};
  else
    error('How could it arrive here?')
  end
end
memb = cellfun(@unique, memb, 'unif', false);   %?!
trees = memb;
end  % end of main


function [new, pth, cmp_pth] = getnewpth(adj_tmp, g_cur, fro, dst, cmp_pth, pth_ini, djk_opt)
n_i = size(dst,1);  % number of "dst" vertices
new = zeros(n_i,6);
pth = cell(n_i,1);
vf = fro(1,1);
cf = fro(1,2);
vt = dst(:,1);
ct = dst(:,2);
[d_i, p_i] = dijkstra_sp(adj_tmp, vf, djk_opt);

eom = 1;
for i = 1:1:n_i
  if isinf(d_i(vt(i)))
    continue;
  end
  new(eom,1) = vf; new(eom,2) = vt(i); new(eom,3) = d_i(vt(i));
  new(eom,4) = g_cur; new(eom,5) = cf; new(eom,6) = ct(i);
  eop = 0; path = pth_ini; strt = vt(i);  % Inlined EXTRACT_PATH()
  while strt ~= 0
    eop = eop+1; path(eop) = strt; strt = p_i(strt);
  end
  pts = path(2:eop-1); pth{eom} = pts;
  cmp_pth(g_cur,pts) = cmp_pth(g_cur,pts)+1;  % update point-component statistics
  eom = eom+1;
end
new(eom:end,:) = [];
pth(eom:end,:) = [];
end


function lin_idx = getadjedg(adj, nod)
num_cur = numel(nod);
edg_nbr = [];
for k = 1:1:num_cur
  nod_nbr = find(any(adj(:,nod(k)),2));
  num_nbr = numel(nod_nbr);
  edg_nbr = [edg_nbr; [repmat(nod(k),num_nbr,1), nod_nbr]];
end
edg_nbr = [edg_nbr; [edg_nbr(:,2),edg_nbr(:,1)]];
lin_idx = sub2ind(size(adj), edg_nbr(:,1), edg_nbr(:,2));
%val = adj(lin_idx); val = full(val);  % sparsity can be confusing here
end

function path = extract_path(i, pred)
path = zeros(size(pred));
eov = 0;
while i ~= 0
  eov = eov+1;
  path(eov) = i;
  i = pred(i);
end
path = path(2:eov-1);
end
