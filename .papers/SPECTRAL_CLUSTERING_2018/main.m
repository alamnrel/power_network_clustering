function main(cases)
% 
% This file performs some case studies from the reference paper for the 
% power network test cases from the MATPOWER toolbox.
% 
% The names of MATPOWER test cases combined in a cell array of strings are 
% the only required input. Some example inputs are given below: 
% 
% cases_small = {'case39', 'case118', 'case300'}; 
% cases_large = {'case1354pegase', 'case2869pegase', 'case2383wp', ...
%    'case6495rte', 'case3120sp',  'case9241pegase'};
% cases_all = {'case30', 'case39', 'case57', 'case118', 'case300', ...
%   'case2383wp', 'case2736sp', 'case2737sop', 'case2746wop', ...
%   'case2746wp', 'case3012wp', 'case3120sp', 'case89pegase', ...
%   'case1354pegase', 'case2869pegase', 'case9241pegase', 'case6495rte'};  % for reference 
% 
% Reference: I. Tyuryukanov, M. Popov, M.A.M.M. van der Meijden, and V. 
% Terzija. "Discovering Clusters in Power Networks From Orthogonal  
% Structure of Spectral Embedding". In: IEEE Trans. Power Syst. 33.6 
% (2018), pp. 6441-6451.
% 

close all;
if nargin<1
  error(['Please provide the names of the MATPOWER testcases in a cell',...
    ' array of strings.']);
end

% Disable selected warnings
old_warn = warning;   % save warning state
cleanup0 = onCleanup(@() warning(old_warn));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:eigs:SigmaNearExactEig');
warning('off', 'BaseIn:InvalidFolder');
warning('off', 'max_flow_mex:cutValueNotFlowValue');

% Amend MATLAB path
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath('netw_data',fileparts(fileparts(fileparts(mfilename('fullpath')))),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'zp_clustering'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'matlab_bgl'),... 
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'graclus1.2', 'matlab'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'balancedKCuts_V1_1'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'balancedKCuts_V1_1','Ncut_9'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'balancedKCuts_V1_1','GraphDemos'));

% INITIALIZE SEVERAL CONSECUTIVE TESTCASES 
rng('default'); 
study(1, numel(cases)) = MatpowerIn();
for i = 1:1:numel(cases)
  study(1, i) = MatpowerIn('caseid', cases{i}, 'n_pf', 1, 'f0', 60,...
     'dev_pq', 0, 'mean_pq', 1);
end

% Iterate through cases
pp_siz = 0.2;              % cluster size constraint (e.g., no less than 20% of average size)
k_start = 2;               % lowest number of clusters
n_tst = 38;                % highest number of clusters
x_width = 317;             % figure size
y_width = 200;             % figure size
n_stp = 2;                 % x-axis step in figures
spcopt.card_min = 1;      
spcopt.improv = true;
spcopt.thrs = sqrt(3)/2;
for obj = study
  mpw = mp2bgl( obj );
  mpw = ext2int( mpw );    % to ensure consistent indexing
  pst = mp2pst( mpw );    
  if ~isempty(pst.bus)
    g = pst2graph(pst, 'mode', 'dcpf');
  else
    obj.curr = obj.curr + 1;
    failed = failed + 1;
    continue;
  end
  assert(nnz(g.adj<0)==0);  % reject cases with capacitive admittances
  g = abs_adj(g);           % usual spectral clustering doesn't suppose negative weights

  m = size(g.adj, 2);
  g.ml = []; g.coh = [];    
  merg = {(1:1:m)'}; 	
  tic;    
  Asym = GraphUtils.createLaplacian(g.adj, 'sym');
  Asym(1:m+1:m*m) = 0;  Asym = -Asym;  % Anorm for > speed and accuracy
  assert(issymmetric(Asym) && ~any(Asym(1:m+1:m*m)))
  [V_full, v] = GraphUtils.spClust(Asym, 'la', n_tst);
  time_spc = toc;

  % THE CASE STUDY TO SHOW THE SAME SHAPE OF YS AND ZM OPTIMISATION OBJECTIVES
  %{
  % The final method ZM:
  tic; [Vr5, ZM5, nc5] = GraphUtils.nc_ev_sym(V_full, [1 1], [1 0 0 0], [], 2, true); t_exe12 = toc;
  fprintf(['Run time for 2..%d clusters with %s ', 'method was %f seconds.\n'], n_tst, ['(prev.rot. & rand.orth. & mult.YuShi -> ', dm2str([0 0 0 1]),')'], t_exe12);
  % The final method YS:
  tic; [Vr, YS5, nc5] = GraphUtils.nc_ev_sym(V_full, [1 1], [0 0 0 1], [], 2, true); t_exe11 = toc;
  fprintf(['Run time for 2..%d clusters with %s ', 'method was %f seconds.\n'], n_tst, ['(prev.rot. & rand.orth. & mult.YuShi -> ', dm2str([0 0 0 1]),')'], t_exe11);    
  Qt = NaN(1, numel(Vr));   
  for i = 1:1:numel(Vr)
    Qt(i) = ev_qt(Vr{i});
  end    
  figure;
  subplot(2,2,1); plot(2:1:n_tst, ZM5,  'ko-', 'LineWidth', 1.0, 'MarkerSize', 4);
  xlabel('Nr. of eigenvectors, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$J^*$, [-]', 'Interpreter', 'Latex', 'FontSize', 12);  % (\bf{Y}\bf{R}^*)
  xlim([2 n_tst]);    
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gca, 'xgrid', 'on');    
  subplot(2,2,2); plot(2:1:n_tst, YS5,  'ko-', 'LineWidth', 1.0, 'MarkerSize', 3);
  xlabel('Nr. of eigenvectors, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$Q^*$, [-]', 'Interpreter', 'Latex', 'FontSize', 12); 
  xlim([2 n_tst]);    
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gca, 'xgrid', 'on');     
  subplot(2,2,3); plot(2:1:n_tst, v(2:end),  'ko-', 'LineWidth', 1.0, 'MarkerSize', 2);
  xlim([2,n_tst]);
  ylim([min(v)-0.0005 1]);
  xlabel('Eigenvalue Nr., [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('Eigenvalue, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gca, 'xgrid', 'on');    
  v = 1 - v;
  subplot(2,2,4); plot(2:1:n_tst-1, abs(diff(v(2:end))./v(2:end-1)),  'ko-', 'LineWidth', 1.0, 'MarkerSize', 2);
  xlabel('Eigenvalue Nr., [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('Rel. Eigengap, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gca, 'xgrid', 'on'); 
  set(gcf, 'pos', [100 100 2.2*x_width 2.2*y_width]);
  close all;
  %}

  % THE CASE STUDY TO PARTITION IEEE 39 DIFFERENTLY 
  % ({'case39'} should be the input argument)
  %{
  tic; [Vr, YS1, nc1] = GraphUtils.nc_ev_sym(V_full, [1 1], [0 0 0 1], [], 2, true); t_exe12 = toc;
  fprintf(['Run time for 2..%d clusters with %s ', 'method was %f seconds.\n'], n_tst, ['(prev.rot. & rand.orth. & mult.YuShi -> ', dm2str([0 0 0 1]),')'], t_exe12);
  Qt = NaN(1, numel(Vr));
  for i = 1:1:numel(Vr), Qt(i) = ev_qt(Vr{i}); end 
  Vr = Vr{2};  
  core_nodes = spec_cores2(Vr, g.vw, g.adj, spcopt);
  [T2] = g.uckwaycut(core_nodes, 2, spcopt.card_min);
  [T3] = GraphUtils.sp_postpr(Vr, 'kmeans', 'sym+', size(Vr,2));
  [T2] = contig(g, T2, 'SplitLargeMinor', true);
  [T3] = contig(g, T3, 'SplitLargeMinor', true);
  [cut2, bus2, T2] = final_cutset(g, T2, merg);
  [cut3, bus3, T3] = final_cutset(g, T3, merg);
  [eXp2, pcut2] = g.ici_info( cut2 );
  [eXp3, pcut3] = g.ici_info( cut3 );
  %}

  % THE CASE STUDY TO COMPARE PARTITIONING QUALITY WITH KMEANS AND HSC
  phi_max = NaN(4,n_tst-1);
  run_exe = NaN(5,n_tst-1);
  cut_tot = NaN(4,n_tst-1);
  num_min = NaN(4,n_tst-1);
  nrm_cut = NaN(4,n_tst-1);
  eXp_clu = cell(4,n_tst-1);
  cut_clu = cell(4,n_tst-1);
  lbl_clu = cell(4,n_tst-1);      

  tic; [Vr, YS1, t_run] = GraphUtils.nc_ev_sym(V_full, [1 1], [0 0 0 1], [], 2, false); t_exe12 = toc;      
  fprintf(['Run time for 2..%d clusters with %s ', 'method was %f seconds.\n'], n_tst, ['(prev.rot. & rand.orth. & mult.YuShi -> ', dm2str([0 0 0 1]),')'], t_exe12);       
  CQ = ev_qt(Vr);   % alignment quality of normalized spectral embedding
  card_core = cell(1,n_tst-1);
  for k = 1:1:n_tst-1
    V = Vr{k};
    Vd = V > 0.91; 
    card_core{k} = sum(Vd, 1);   % lower estimates of each cluster size
  end      
  col_curr = k_start-1;
  col_max = n_tst-1;
  while col_curr<=col_max
    nc_curr = col_curr+1;
    spcopt.card_min = max(2, ceil(m/nc_curr*pp_siz));
    for p = col_curr:1:col_max
      sel = card_core{p};          
      if nnz(sel>=spcopt.card_min) >= nc_curr
        col_sel = p;
        sel = card_core{p};
        [val, sel] = sort(sel, 'descend');
        sel = sel(1:nc_curr);
        break;
      else
        col_sel = [];
      end
    end
    if isempty(col_sel)
      n_tst = col_curr;  % as n_tst-1 is what is actually used
      break;
    end
    V = Vr{col_sel};     
    tic;      
    core_nodes = spec_cores2(V(:,sel), g.vw, g.adj, spcopt);        
    [T2] = g.uckwaycut(core_nodes, 2, spcopt.card_min);              
    time_our = toc;            
    tic;
    [dummy1, dummy2] = GraphUtils.spClust(Asym, 'la', n_tst);
    time_spec = toc;                
    [T2] = contig(g, T2, 'SplitLargeMinor', true);
    [cut2, bus2, T2] = final_cutset(g, T2, merg);        
    [eXp2, pcut2, ixs2] = g.ici_info( cut2 );
    n_avg = m/nc_curr;
    run_exe(2,col_curr) = time_our;
    run_exe(5,col_curr) = time_spec;
    phi_max(2,col_curr) = max(eXp2(1,:));
    nrm_cut(2,col_curr) = sum(eXp2(1,:))/nc_curr;
    cut_tot(2,col_curr) = pcut2;
    num_min(2,col_curr) = min(eXp2(3,:))/n_avg*100;               
    [eXp, ixs, pcut] = GraphUtils.lbl_info(g.adj, g.vw, T2);
    eXp_clu{2,col_curr} = eXp(1,:);
    cut_clu{2,col_curr} = pcut;
    lbl_clu{2,col_curr} = T2(:);
    col_curr = col_curr+1;
  end            

  % Now repeat for the other 3 methods (kmeans, HSC, graclus)
  g.coh = ones(2);   % some random placeholder
  for k = k_start-1:1:n_tst-1
    V = Vr{k};           
    n_clust = k+1;
    g.coh(2,:) = n_clust;
    tic;        
    [T1] = GraphUtils.sp_postpr(V, 'average', 'sym-', n_clust, g.adj);
    time_AHC = toc;
    tic;
    [T3] = GraphUtils.sp_postpr(V, 'kmeans', 'sym+', n_clust);
    time_skm = toc;
    if isunix
      tic; [T4, grac_obj] = graclus(g.adj,n_clust,0,100,0);
      time_graclus = toc;
    else
      T4 = T3; time_graclus = 0;          
      %{
      tic;  % another method as an alternative under Windows...
      [DUMMY,T4] = evalc('balancedKCut(g.adj,n_clust,1,sum(g.adj,2),0,1);');
      time_graclus = toc;     
      %}         
    end                  
    [T1] = contig(g, T1, 'SplitLargeMinor', true);
    [T3] = contig(g, T3, 'SplitLargeMinor', true);
    [T4] = contig(g, T4, 'SplitLargeMinor', true);            
    [cut1, bus1, T1] = final_cutset(g, T1, merg);
    [cut3, bus3, T3] = final_cutset(g, T3, merg);
    [cut4, bus4, T4] = final_cutset(g, T4, merg);        
    [eXp1, pcut1, ixs1] = g.ici_info( cut1 );
    [eXp3, pcut3] = g.ici_info( cut3 );
    [eXp4, pcut4] = g.ici_info( cut4 );

    n_avg = m/(k+1);        
    run_exe(1,k) = time_AHC;
    run_exe(3,k) = time_skm; 
    run_exe(4,k) = time_graclus; 
    phi_max(1,k) = max(eXp1(1,:));       
    phi_max(3,k) = max(eXp3(1,:));
    phi_max(4,k) = max(eXp4(1,:));
    nrm_cut(1,k) = sum(eXp1(1,:))/(k+1);        
    nrm_cut(3,k) = sum(eXp3(1,:))/(k+1);
    nrm_cut(4,k) = sum(eXp4(1,:))/(k+1);
    cut_tot(1,k) = pcut1;        
    cut_tot(3,k) = pcut3;
    cut_tot(4,k) = pcut4;
    num_min(1,k) = min(eXp1(3,:))/n_avg*100;        
    num_min(3,k) = min(eXp3(3,:))/n_avg*100;
    num_min(4,k) = min(eXp4(3,:))/n_avg*100;                 
    [eXp, ixs, pcut] = GraphUtils.lbl_info(g.adj, g.vw, T3);
    eXp_clu{3,k} = eXp(1,:);
    cut_clu{3,k} = pcut;
    lbl_clu{3,k} = T3(:)';                        
    [eXp, ixs, pcut] = GraphUtils.lbl_info(g.adj, g.vw, T4); 
    eXp_clu{1,k} = eXp(1,:);
    cut_clu{1,k} = pcut;
    lbl_clu{1,k} = T4(:)';         
  end    
  t_run = t_run(:,1:n_tst-1);
  run_exe = run_exe(:,1:n_tst-1);
  run_exe(3,1) = run_exe(3,2);
  run_exe(1:3,:) = bsxfun(@plus, run_exe(1:3,:), run_exe(5,:));
  phi_max = phi_max(:,1:n_tst-1);
  nrm_cut = nrm_cut(:,1:n_tst-1);
  cut_tot = cut_tot(:,1:n_tst-1);
  num_min = num_min(:,1:n_tst-1);
  run_exe(2,:) = run_exe(2,:) + t_run;    

  hf = figure; hold on;
  plot(2:1:n_tst, num_min(3,:), 'ro--', 'MarkerSize', 4, 'LineWidth', 1);
  plot(2:1:n_tst, num_min(4,:), '*:', 'MarkerSize', 4, 'LineWidth', 1, 'Color', rgb('ForestGreen'));      
  plot(2:1:n_tst, num_min(1,:), 'bs-.', 'MarkerSize', 4, 'LineWidth', 1);
  plot(2:1:n_tst, num_min(2,:), 'k.-', 'MarkerSize', 10, 'LineWidth', 1);      
  xlabel('Number of clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  % ylabel('$\frac{\mathrm{min}_{1\leq i \leq k}(|V_i|)}{|V|/k}$, [-]', 'Interpreter', 'Latex', 'FontSize', 16);
  ylabel('$\varepsilon$, [\%]', 'Interpreter', 'Latex', 'FontSize', 13)      
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gcf, 'pos', [100 100 x_width y_width]);
  ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'off'); ax.GridAlpha = 0.25;
  columnlegend(2,{'k-means','Graclus','HSC', 'our'}, 'Box', 'off');
  close;      

  hf = figure; hold on;
  plot(2:1:n_tst, phi_max(3,:)./phi_max(2,:), 'ro--', 'MarkerSize', 4, 'LineWidth', 1);
  plot(2:1:n_tst, phi_max(4,:)./phi_max(2,:), '*:', 'MarkerSize', 4, 'LineWidth', 1, 'Color', rgb('ForestGreen'));      
  plot(2:1:n_tst, phi_max(1,:)./phi_max(2,:), 'bs-.', 'MarkerSize', 4, 'LineWidth', 1); % 'k.-',
  xlabel('Number of clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$\phi_{max} / \phi_{max}^{\dagger}$, [-]',  'Interpreter', 'Latex', 'FontSize', 13);
  set(gca, 'yscale', 'log');
  plot([2 n_tst], [1 1], 'k-', 'LineWidth', 1);
  set(gcf, 'pos', [100 100 x_width y_width]);
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  legend({'k-means','Graclus','HSC'}, 'Box', 'off', 'Location', 'Best');
  ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'on'); ax.GridAlpha = 0.5; ax.MinorGridAlpha = 0.5;            
  close;

  hf = figure; hold on;
  plot(2:1:n_tst, nrm_cut(3,:)./nrm_cut(2,:), 'ro--', 'MarkerSize', 4, 'LineWidth', 1);      
  plot(2:1:n_tst, nrm_cut(4,:)./nrm_cut(2,:), '*:', 'MarkerSize', 4, 'LineWidth', 1, 'Color', rgb('ForestGreen'));
  plot(2:1:n_tst, nrm_cut(1,:)./nrm_cut(2,:), 'bs-.', 'MarkerSize', 4, 'LineWidth', 1); % 'k.-'
  xlabel('Number of clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$Ncut / Ncut^{\dagger}$, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
  set(gca, 'yscale', 'linear');
  plot([2 n_tst], [1 1], 'k-', 'LineWidth', 1);
  set(gcf, 'pos', [100 100 x_width y_width]);
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  legend({'k-means','Graclus','HSC'}, 'Box', 'off', 'Location', 'Best');
  ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'off'); ax.GridAlpha = 0.25;
  close;

  hf = figure; hold on;
  plot(2:1:n_tst, cut_tot(3,:)./cut_tot(2,:), 'ro--', 'MarkerSize', 4, 'LineWidth', 1);  
  plot(2:1:n_tst, cut_tot(4,:)./cut_tot(2,:), '*:', 'MarkerSize', 4, 'LineWidth', 1, 'Color', rgb('ForestGreen'));      
  plot(2:1:n_tst, cut_tot(1,:)./cut_tot(2,:), 'bs-.', 'MarkerSize', 4, 'LineWidth', 1);
  xlabel('Number of clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$cut / cut^{\dagger}$, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
  set(gca, 'yscale', 'linear');      
  plot([2 n_tst], [1 1], 'k-', 'LineWidth', 1);
  set(gcf, 'pos', [100 100 x_width y_width]);
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  legend({'k-means','Graclus','HSC'}, 'Box', 'off', 'Location', 'Best');
  ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'off'); ax.GridAlpha = 0.25;      
  close;
  %}

  hf = figure; hold on;
  plot(2:1:n_tst, run_exe(3,:), 'ro--', 'MarkerSize', 4, 'LineWidth', 1);  
  plot(2:1:n_tst, run_exe(4,:), '*:', 'MarkerSize', 4, 'LineWidth', 1, 'Color', rgb('ForestGreen'));
  plot(2:1:n_tst, run_exe(1,:), 'bs-.', 'MarkerSize', 4, 'LineWidth', 1);
  plot(2:1:n_tst, run_exe(2,:), 'k.-', 'MarkerSize', 10, 'LineWidth', 1);
  columnlegend(2,{'k-means','Graclus','HSC', 'our'}, 'Box', 'off');
  xlabel('Number of clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
  ylabel('$t_{run}$, [s]', 'Interpreter', 'Latex', 'FontSize', 14);
  xlim([2 n_tst]);
  set(gca, 'XTick', 2:n_stp:n_tst);
  set(gcf, 'pos', [100 100 x_width y_width]);
  set(gca, 'yscale', 'log');      
  ax = gca; ax.GridColor = [0, 0, 0]; set(ax, 'ygrid', 'on'); set(gca, 'yminorgrid', 'on'); ax.GridAlpha = 0.5; ax.MinorGridAlpha = 0.5;                  
  close;
  %}    
end
end

function mth_str = dm2str(dm)
if dm(1) == 1
  mth_str = 'NAG_batch';
elseif dm(2) == 1
  mth_str = 'adam_MBGD';
elseif dm(3) == 1
  mth_str = 'ZMP_SGD';
elseif dm(4) == 1
  mth_str = 'Yu-Shi';
else
  mth_str = 'none';
end
end