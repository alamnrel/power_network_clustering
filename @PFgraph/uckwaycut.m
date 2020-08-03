function T = uckwaycut(g, core_nodes, method, card_min)

% 
% For recursive bisection methods, cluster cores in core_nodes must be
% provided in the order of certainty that the corresponding "core" 
% represents a complete cluster.
% 

old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party'),...
  fullfile(fileparts(fileparts(mfilename('fullpath'))),'third_party', 'matlab_bgl'));

if nargin < 4
  card_min = 1;
end

%g.bus
if nargin < 3
  method = 2;
end

num_clst = numel(core_nodes);
g_red = g.merge_nodes(core_nodes);
merge_map = g_red.merge_map;
m = size(g_red.adj, 1);
g_red.coh = [1:m; ones(1,m)];
g_red.coh(2,end-num_clst+1:end) = 2:num_clst+1;

vw = tabulate(merge_map);
vw = vw(:,2);
core_card = cellfun(@numel, core_nodes);
assert( all(core_card==vw(end-num_clst+1:end)) );

switch method
  case 1   % multiway recursive-bisection Ncut
    vs = m-num_clst+1:m;
    T = GraphUtils.min_cut_rb(g_red.adj, vs, g_red.vw, true);    
    T = T(merge_map);    
  
  case 2   % multiway cut methods 
    %{
    % hMetis2
    GraphUtils.pfgraph2hmetis(g_red.adj, g_red.inc, ones(m,1), g_red.coh);
    system(sprintf('hmetis2.0pre1 test_hmetis.hgr %d -rtype=slow -ufactor=%d -nruns=30 -nvcycles=30 -kwayrefine -fixed=fixfile.txt >/dev/null 2>/dev/null', num_clst, 1000000));
    % system(sprintf('hmetis2.0pre1 test_hmetis.hgr %d -ptype=rb -ctype=edge1 -rtype=slow -ufactor=%d -nruns=30 -nvcycles=20 -kwayrefine -fixed=fixfile.txt', num_clst, 100000)); %-kwayrefine
    T = transpose(load(sprintf('test_hmetis.hgr.part.%d',num_clst))) + 1;
    %}
    %{
    % Multiway cut: (2-2/k)OPT approximation by the isolating cut heuristic
    vs = m-num_clst+1:m;
    T = GraphUtils.mwc_2OPT(g_red.adj, vs);  %, true
    %}
    
    % Multiway cut via recursive mincut
    vs = m-num_clst+1:m;
    T = GraphUtils.min_cut_rb(g_red.adj, vs, true);    
    %T = GraphUtils.max_flow_rb({g_red, g.adj}, vs, 10, true);
    %}
    
    T = T(merge_map);
  case 3   % simple k-way separation of k terminals based on spectral shortest paths
    g_red = g_red.createLaplacian('sym');
    assert(issymmetric(g_red.Lsym) && ~any(g_red.Lsym(1:m+1:m*m)));
    [V_red, v] = GraphUtils.spClust(g_red.Lsym, 'sym', num_clst);
	V_red = Utils.normalize_rows(V_red);
    emb_graph = GraphUtils.buildDistanceGraph(V_red, g_red.adj, 'sym');
    dm = GraphUtils.nShortestPaths(emb_graph, g_red.coh(1,end-num_clst+1:end), 'orig');
    [~, T] = min(dm, [], 2);
    T = T(merge_map);
  otherwise
    error('');
end
end