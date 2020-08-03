function main(mode)
%
% This is a demo study to demonstate that graph outlier removal (both of  
% weakly connected nodes and too small clusters) allows to avoid a very 
% unbalanced graph partitioning if it is naturally present. 
%  
% However, it must be noted that small clusters resulting from graph
% partitioning cannot be entirely avoided by removing graph outliers. Graph
% partitioning results also depend on the input graph/network and the used
% partitioning method. Thus small clusters may occur even if graph outliers
% have been preprocessed. However, such small clusters usually have worse 
% cluster separation metrics than those required from outliers (i.e., they
% are usually not outliers but merely "small clusters"). 
%
% To automate this demo, the outlier size and separation metrics max_card 
% and score_max are computed for each required number of clusters using the 
% formulas from the reference paper. However, due to the reasons explained
% above (above all, the demo nature of this study) these parameters can be
% overridden for certain test networks (this is also mentioned in the 
% reference paper). 
%
% Input parameters:
%   MODE = 'rand' simulates test networks with edge weights changed to
%          random positive real values.
%   MODE = 'classic' simulates test networks with edge weights from random
%          power flows (same as in the reference paper)   
% 
% Reference: I. Tyuryukanov, M. Popov, M. A. M. M. van der Meijden, and V. 
% Terzija. "Spectral MST-based Graph Outlier Detection with Application to 
% Clustering of Power Networks". In: Proc. 20th Power Systems Computation 
% Conference (PSCC 2018). Dublin, Ireland, 2018, pp. 1-8.
%

if nargin<1
  mode = 'classic';
end

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
rng('default');

% INITIALIZE SEVERAL CONSECUTIVE TESTCASES
n_ovrl = 10050;
cases = {'case2383wp'};   % 'case89pegase', 'case39', 'case118', 'case2869pegase', 'case3120sp', 'case6495rte', 'case9241pegase'
n_coh = arrayfun(@(x) ones(1,x), n_ovrl*ones(1,10), 'UniformOutput', false);
n_pf = cellfun(@numel, n_coh, 'UniformOutput', true);
f0 = 60*ones(1, numel(n_coh));
study(1, numel(cases)) = MatpowerIn();
for i = 1:1:numel(cases)
  study(1, i) = MatpowerIn('caseid', cases{i}, 'n_pf', n_pf(i), 'n_coh',...
    n_coh{i}, 'f0', f0(i), 'dev_pq', 0.5*ones(1,n_pf(i)), 'mean_pq',...
    ones(1,n_pf(i)), 'seed', 2);
end

% Limit the number of tested cases at 150 cases per each nr. of clusters
n_test = 150; 
n_clust = [2*ones(1,n_test),3*ones(1,n_test),4*ones(1,n_test),...
  5*ones(1,n_test),6*ones(1,n_test),7*ones(1,n_test),8*ones(1,n_test),...
  9*ones(1,n_test),10*ones(1,n_test),11*ones(1,n_test),12*ones(1,n_test)];

% ITERATE THROUGH CASES
for obj = study
  failed = 0;   % number of failed powerflows  
  % bad_clu_ou: number of small clusters for each case (small regardless of 
  % coupling strength). This criterion is weakened slightly to 0.9*max_card
  % because the absence of any small clusters cannot be strictly enforced. 
  bad_clu_ou = NaN(6,numel(n_clust));                                             
  av_siz_ou = NaN(6,numel(n_clust));   % average small cluster size for each case
  min_siz_ou = NaN(6,numel(n_clust));  % minimal small cluster size for each case 
  pp_clu_ou = NaN(6,numel(n_clust));   % number of real outliers for each case (both small and weakly connected)   
  fail_ou = NaN(1,numel(n_clust));
  objcurr = NaN(1,numel(n_clust));

  while obj.curr <= obj.n_pf
    %k = k+1;
    %obj.curr = sel(k);    
    %failed = fal(k);    
    idx_curr = obj.curr - failed;
    mpw = mp2bgl(obj);
    pst = mp2pst(mpw);  
    if ~isempty(pst.bus)      
      g = pst2graph(pst, 'mode', 'apfg');
      g = abs_adj(g);    % ignore possible negative weights
    else
      obj.curr = obj.curr + 1;
      failed = failed + 1;
      continue;
    end
    
    % If mode=='rand', randomize edge weights (good for testing)
    if strcmp(mode, 'rand')
      rng(obj.curr, 'twister');
      [i, j, ~] = find(g.adj);
      n = numel(i); n1 = round(n/5); n2 = n - n1;
      w1 = abs(10*randn(n2, 1) + 10);
      w2 = abs(0.1*randn(n1, 1));   % trigger some outliers
      w = [w1; w2]; w =  w(randperm(n));
      adj = sparse(i, j, w);
      adj = adj + adj';
      g = PFgraph('adj', adj, 'gen', g.gen);
      g.coh = [1:20; -1*ones(1,20)];
    end             
    num_clust = n_clust(idx_curr);  
    g.coh(2,:) = num_clust;  % set the number of clusters for partitioning    
    m = size(g.adj, 2);
    
    % Find and merge graph outliers below score_max and max_card
    score_max = min(log2(num_clust)*0.02, 0.05);   
    % score_max = 0.05;   % forcefully override for some networks    
    max_card = max(2, round(0.1*m/num_clust));   
    [g_pp, merg] = g.outliers_kwaysp(num_clust, max_card, score_max, 'debug', 0);
    
    % Partition the preprocessed graph with two methods
    [T2] = g_pp.ucSpEmbGrPart('alg', 'sym+', 'postpr', 'kmeans');
    [T2] = contig(g_pp, T2, 'SplitLargeMinor', false);
    [cut2, bus2, T2] = final_cutset(g, T2, merg);
    [T3] = g_pp.ucSpEmbGrPart('alg', 'sym-', 'postpr', 'average');
    [T3] = contig(g_pp, T3, 'SplitLargeMinor', false);    
    [cut3, bus3, T3] = final_cutset(g, T3, merg);
    
    % Partition the original graph with the same two methods
    [T1] = g.ucSpEmbGrPart('alg', 'sym+', 'postpr', 'kmeans');
    [T1] = contig(g, T1, 'SplitLargeMinor', false);
    [cut1, bus1, T1] = final_cutset(g, T1);    
    [T4] = g.ucSpEmbGrPart('alg', 'sym-', 'postpr', 'average'); 
    [T4] = contig(g, T4, 'SplitLargeMinor', false);
    [cut4, bus4, T4] = final_cutset(g, T4);       
    
    % Get the partitioning metrics
    [eXp1, pcut1, ixs1] = g.ici_info( cut1 );
    [eXp2, pcut2, ixs2] = g.ici_info( cut2 );
    [eXp3, pcut3, ixs3] = g.ici_info( cut3 );
    [eXp4, pcut4, ixs4] = g.ici_info( cut4 );
    eXp1 = eXp1(1:3,:); eXp2 = eXp2(1:3,:);
    eXp3 = eXp3(1:3,:); eXp4 = eXp4(1:3,:);    
    
    % Identify clusters smaller than max_card    
    eXp1(3,:) = eXp1(3,:)/max_card; 
    eXp2(3,:) = eXp2(3,:)/max_card; 
    eXp3(3,:) = eXp3(3,:)/max_card; 
    eXp4(3,:) = eXp4(3,:)/max_card;      
    
    % Record statistics for outliers and small clusters 
    pp_clu_ou(1,idx_curr) = nnz(eXp1(1,eXp1(3,:)<=1)<=score_max);
    pp_clu_ou(2,idx_curr) = nnz(eXp2(1,eXp2(3,:)<=1)<=score_max);
    pp_clu_ou(3,idx_curr) = nnz(eXp3(1,eXp3(3,:)<=1)<=score_max);
    pp_clu_ou(4,idx_curr) = nnz(eXp4(1,eXp4(3,:)<=1)<=score_max);    
    av_siz_ou(1,idx_curr) = sum(eXp1(3,eXp1(3,:)<=1))*100/max(1, nnz(eXp1(3,:)<=1) );    
    av_siz_ou(2,idx_curr) = sum(eXp2(3,eXp2(3,:)<=1))*100/max(1, nnz(eXp2(3,:)<=1) );
    av_siz_ou(3,idx_curr) = sum(eXp3(3,eXp3(3,:)<=1))*100/max(1, nnz(eXp3(3,:)<=1) );
    av_siz_ou(4,idx_curr) = sum(eXp4(3,eXp4(3,:)<=1))*100/max(1, nnz(eXp4(3,:)<=1) );         
    min_siz_ou(1,idx_curr) = min(eXp1(3,:));
    min_siz_ou(2,idx_curr) = min(eXp2(3,:));
    min_siz_ou(3,idx_curr) = min(eXp3(3,:));
    min_siz_ou(4,idx_curr) = min(eXp4(3,:));
    bad_clu_ou(1,idx_curr) = nnz(eXp1(3,:)<0.9);
    bad_clu_ou(2,idx_curr) = nnz(eXp2(3,:)<0.9);
    bad_clu_ou(3,idx_curr) = nnz(eXp3(3,:)<0.9);
    bad_clu_ou(4,idx_curr) = nnz(eXp4(3,:)<0.9);     
    fail_ou(1,idx_curr) = failed;
    objcurr(1,idx_curr) = obj.curr;
        
    fprintf('|%d power flow cases (nc=%.0f):  %.2f, %.2f, %.2f, %.2f|\n',...  
      obj.curr, num_clust, [min(eXp3(3,:)), min(eXp2(3,:)), ... 
      min(eXp4(3,:)), min(eXp1(3,:))]);     
    obj.curr = obj.curr + 1;
    
    % Save the results for the current network (e.g. for plotting) & break
    if (idx_curr >= numel(n_clust))
      save(['ou_',obj.caseid,'_',mode,'.mat'], 'pp_clu_ou', 'av_siz_ou',...
        'min_siz_ou', 'bad_clu_ou', 'fail_ou', 'objcurr');              
      break;
    end
  end
end

end
