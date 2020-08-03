function [core_nodes] = quick_cores(V, lbl, adj, card_min, debug)
% 
% Works only on row-normalised and aligned eigenvector matrix V.
% lbl is for debug mostly (to see which buses have what degree of membership).
% 

old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(mfilename('fullpath')),'third_party'),...
  fullfile(fileparts(mfilename('fullpath')),'third_party', 'matlab_bgl'));

if nargin < 5
  debug = false;
end

nc = size(V,2);
core_nodes = cell(nc,1);
cut_offs = Inf(nc,1);
epsilon = eps('double');

for j = 1:1:nc
  [core, cut_off] = clu_est(V(:,j), lbl, adj, card_min, sqrt(2)/2+10*epsilon, debug);
  core_nodes{j} = core;
  cut_offs(j) = cut_off;
end
idx_del = isempty(core);

% here and now (!) we assume that small cores have been filtered in advance
assert(~any(idx_del));  
core_nodes(idx_del) = [];
V(:,idx_del) = [];
cut_offs(idx_del) = [];
end


function [core, cut_off] = clu_est(vec, lbl, adj, card_min, thrs_low, debug)
  m = numel(vec);
  if nargin < 6
    debug=false;
  end
  
  [vec_lrge, idx_lrge] = sort(vec, 'descend');    
  if vec_lrge(1)-0.1<thrs_low
    thrs_low = max(0.75, vec_lrge(1)-0.1);
  end
  idx_fnsh = find(vec_lrge<thrs_low, 1) - 1;
  idx_stop = find(vec_lrge>sqrt(2)/2+1e-12, 1, 'last');
  if idx_stop>=card_min
    idx_curr = idx_lrge(1:card_min);    
  else
    core = [];
    cut_off = 0;
    return;
  end    
  num_save = numel(idx_curr);
  idx_sweep = idx_lrge(card_min+1:idx_fnsh);
  num_save = num_save + numel(idx_sweep); 

  % Here it could be simplified to leaving one connected component only!  
  core = idx_lrge(1:1:num_save);
  % Graph vertices are extracted as in lbl(core): adj_core(2,3)<->{lbl(core(2)),lbl(core(3))}<->{core(2),core(3)}
  adj_core = adj(lbl(core),lbl(core));
  [n_cc, idx_cc] = GraphUtils.alecconncomp(adj_core);
  card_cc = tabulate(idx_cc);
  [card_max, icc_max] = max(card_cc(:,2));
  cc_max = card_cc(icc_max,1);  
  core = core(idx_cc==cc_max);
  adj_core = adj(lbl(core),lbl(core));
  [n_cc, idx_cc] = GraphUtils.alecconncomp(adj_core);  
  assert(n_cc==1);
  assert(size(adj_core,1)==card_max);
  
  if debug
    x = 1:1:num_save;
    plot(x, vec_lrge(1:num_save), 'k.-');
    xl = xlim;
    hold on;
    plot([xl(1),xl(2)], [vec(core(end)), vec(core(end))], 'k-.');
    plot(1:1:numel(core),sort(vec(core),'descend'),'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
    for k = 1:num_save
        if iscell(lbl(idx_lrge(k)))
          row_lbl = lbl{idx_lrge(k)};
        else
          row_lbl = lbl(idx_lrge(k));
        end
        text(k, vec_lrge(k), num2str(row_lbl), 'FontSize', 9, 'Color', 'k', 'FontWeight',...
            'bold', 'VerticalAlignment','top', 'HorizontalAlignment', 'center');
    end    
    hold off;
    ylim([max(0.707,vec_lrge(num_save+1)) 1]);
    % set(gca, 'Position', [0.1300    0.1762    0.6    0.2750])
    xlabel('Node number, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel('Eigenvector, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  end
  cut_off = min(vec_lrge(ismember(idx_lrge,core)));  
end

