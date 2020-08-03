function part_viz(adj, vx_lbl, name, dir)

% Amend MATLAB path & cleanup
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath('third_party');

if isa(adj, 'PFgraph')
  g = adj;
else
  g = PFgraph('adj', adj);
end
n = length(vx_lbl);

if Utils.isint(vx_lbl)   % if cluster labels
  grp = unique(vx_lbl);
  nc = numel(grp);
  col = distinguishable_colors(nc,[1 1 1; 0 0 0]);  
  blu = find(col(:,1)==0 & col(:,2)==0 & col(:,3)==1);
  col(blu,:) = [0, 1, 1];  % set to orange 1, 0.647, 0]
  vx_col = zeros(n,3);
  for i = 1:1:nc
    idx = vx_lbl==grp(i);
    grp_siz = nnz(idx);
    vx_col(idx,:) = repmat(col(i,:),grp_siz,1);
  end  
else   % if "nodal potentials"
  map = cool(100);   % colormap
  map = flipud(map);   % blue colors on bottom
  k = (100-1)/(min(vx_lbl)-max(vx_lbl));
  b = 1-k*max(vx_lbl);
  col_idx = round(k*vx_lbl + b);
  vx_col = map(col_idx,:);
end

if nargin<3, name = 'test_graph'; end
if nargin<4, dir = pwd; end
file = fullfile(dir, name);
g_py = g.pfgraph2igraph(vx_col);
py.igraphmod.igraph2graphviz(g_py, file);
GraphUtils.graphviz(dir);
end
