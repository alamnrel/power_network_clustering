function runtimes()
%
% This file measures the runtime of MST-based graph outlier detection on 
% several test networks of various sizes.
% 
% Reference: I. Tyuryukanov, M. Popov, M. A. M. M. van der Meijden, and V. 
% Terzija. "Spectral MST-based Graph Outlier Detection with Application to 
% Clustering of Power Networks". In: Proc. 20th Power Systems Computation 
% Conference (PSCC 2018). Dublin, Ireland, 2018, pp. 1-8.
% 

% Disable selected warnings
old_warn = warning;  % save warning state
cleanup0 = onCleanup(@() warning(old_warn));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:eigs:SigmaNearExactEig');
warning('off', 'BaseIn:InvalidFolder');
warning('off', 'pst2graph:IllCondInput');

% Amend MATLAB path
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath(fullfile('..','..'),...
  fullfile('..','..','third_party'),...
  fullfile('..','..','third_party', 'matlab_bgl'),...
  fullfile('..','..','third_party', 'gaimc'));

% INITIALIZE SEVERAL CONSECUTIVE TESTCASES
n_test = 1;  
n_clust = [2*ones(1,n_test),3*ones(1,n_test),4*ones(1,n_test),...
  5*ones(1,n_test),6*ones(1,n_test),7*ones(1,n_test),8*ones(1,n_test)];
n_ovrl = numel(n_clust);

cases = {'case300', 'case1354pegase',  'case2383wp',...
  'case2869pegase', 'case3120sp', 'case6495rte', 'case9241pegase',...
  'case13659pegase'};   %  
n_coh = arrayfun(@(x) ones(1,x), n_ovrl*ones(1,10), 'UniformOutput', false);
n_pf = cellfun(@numel, n_coh, 'UniformOutput', true);
f0 = 60*ones(1, numel(n_coh));
study(1, numel(cases)) = MatpowerIn();
for i = 1:1:numel(cases)
  study(1, i) = MatpowerIn('caseid', cases{i}, 'n_pf', n_pf(i), 'n_coh',...
    n_coh{i}, 'f0', f0(i), 'dev_pq', 0*ones(1,n_pf(i)), 'mean_pq',...
    ones(1,n_pf(i)), 'seed', 2);
end

% Iterate through cases
timing = NaN(numel(study),numel(n_clust)+1);
score_max = 0.03;

for k = 1:1:numel(study)
  obj = study(k);
  failed = 0;   % number of failed powerflows  
  while obj.curr <= obj.n_pf    
    idx_curr = obj.curr - failed;
    mpw = mp2bgl(obj);
    pst = mp2pst(mpw);   
    if ~isempty(pst.bus)
      g = pst2graph(pst); 
    else
      obj.curr = obj.curr + 1;
      failed = failed + 1;
      continue;
    end
    m = size(g.adj, 2);
    
    num_clust = n_clust(idx_curr);  
    g.coh(2,:) = num_clust;
    max_card = max(2, round(0.1*m/median(n_clust)));
    
    f = @() g.outliers_kwaysp(num_clust, max_card, score_max, 'debug', 0);
    t_curr = timeit(f,0);  % 0 outputs is critical
    
    timing(k,idx_curr) = t_curr;
    obj.curr = obj.curr + 1;
  end
  timing(k,idx_curr+1) = m;
end
save('timing_all_nodes.mat', 'timing');
end
