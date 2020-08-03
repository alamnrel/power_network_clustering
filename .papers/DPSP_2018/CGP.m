function varargout = CGP(pst, ncl_ref, mth_ref)
%
% Constrained graph partitioning.
%
% pst - Power network data in the PST format (Power System Toolbox by JH 
% Chow and G Rogers)
% ncl_ref - numbers of groups to split (can be a scalar or vector)
% mth_ref - partitioning methods to use (see the paper)
%     mth_ref(1) - TCRBGC  
%     mth_ref(2) - hMetis
%     mth_ref(3) - our tree heuristic
%
% Reference: I. Tyuryukanov, A.C. Karagiannis, M. Popov, M.A.M.M. van der 
% Meijden, and V. Terzija. "Generator grouping cutset determination based 
% on tree construction and constrained spectral clustering". In: The 
% Journal of Engineering 2018.15 (2018), pp. 1309-1314.
%

if nargout>1
  num_nc = numel(ncl_ref);
  clk_nc = Inf(2,num_nc);
end

% Iterate through the cases
case_stat = zeros(numel(ncl_ref), 3);
for k = 1:1:numel(ncl_ref)   % number of clusters to be created
  nc = ncl_ref(k);  
  [pst, Q, labels] = cohConstraints(pst, nc);
  g = pst2graph(pst, 'mode', 'apfg');   % aim for min power flow disruption
  warning('off'); 
  m  = size(g.adj, 1);
  ng = size(g.coh, 2);
  g.ml = [];  % deliberately disable as non-essential
  % m = size(g.adj, 1); T = zeros(1,m); T(g.coh(1,:)) = g.coh(2,:)+1; part_viz(g.adj, T);  
  
  % Constrained Spectral Clustering
  if any(mth_ref(3))
    f1 = @() GraphUtils.simadj2distadj( g.adj, 'inv1' );
    adjD = GraphUtils.simadj2distadj( g.adj, 'inv1' );
    gD = PFgraph('adj', adjD, 'inc', g.inc, 'merge_map', g.merge_map, 'coh',...
      g.coh, 'ml', g.ml, 'bus', g.bus, 'gen', g.gen);  % double(logical(g.adj)) 
    if nargout==1
      trees = hclust_constrained_path(gD);   % double(logical(g.adj))  
      T_cs = g.coGrReduClust(trees);   
    else
      f2 = @() hclust_constrained_path(gD);
      trees = hclust_constrained_path(gD);   % double(logical(g.adj))  
      T_cs = g.coGrReduClust(trees);
      clk_nc(1,k) = timeit(f1,1);
      clk_nc(2,k) = timeit(f2,1);
    end
    [cut_cs, ~, T_cs] = final_cutset(g, T_cs);
    [eXp_cs, pcut_cs, ~, adj_sep, num_viols_cs, viols_cs, tot_outl_gen_cs] = g.ici_info( cut_cs );
    %part_viz(adj_sep, T_cs);
  else
    tot_outl_gen_cs = -1;
  end  
  
  % TCRBGC
  if any(mth_ref(1))    
    n_try = 15;
    [DUMMY, T_bc] = evalc('ourTransductiveBalancedKCut(g, size(labels,2), labels, 1, ones(size(g.adj,1),1), 0, 12, 0);');   % ourT, n_try, false
    [cut_bc, ~, T_bc] = final_cutset(g, T_bc);
    [eXp_bc, pcut_bc, ~, ~, num_viols_bc, viols_bc, tot_outl_gen_bc] = g.ici_info( cut_bc );    
  else
    tot_outl_gen_bc = -1;
  end  
  
  % hMetis (from Linux)
  if any(mth_ref(2))
    [T_hm, UBfactor, contig] = ucHMETIS(g, 0, 1);
    [cut_hm, ~, T_hm] = final_cutset(g, T_hm);
    [eXp_hm, pcut_hm, ~, ~, num_viols_hm, viols_hm, tot_outl_gen_hm] = g.ici_info( cut_hm );
  else
    tot_outl_gen_hm = -1;
  end   
    
  if tot_outl_gen_cs<=0
    fprintf('islands = %d, fail_cs = %d | fail_hm = %d | fail_bc = %d\n', nc, tot_outl_gen_cs, tot_outl_gen_hm, tot_outl_gen_bc);
  else
    cprintf('red', 'islands = %d, fail_cs = %d | fail_hm = %d | fail_bc = %d\n', nc, tot_outl_gen_cs, tot_outl_gen_hm, tot_outl_gen_bc);
  end
  case_stat(k,:) = [tot_outl_gen_cs, tot_outl_gen_hm, tot_outl_gen_bc];
end

varargout{1} = case_stat;
if nargout>1
  varargout{2} = clk_nc;
end
end

