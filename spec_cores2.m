function [core_nodes, eXp] = spec_cores2(Z, vw, adj, opts, lbl)
% 
% Works on axis-aligned eigenvector matrix Z obtained from the inv(M)*K 
% electromechanical model. It may work or may not work on similar models 
% of the type inv(D)*L ("random-walk Laplacian" models).
% lbl is for debug mostly (to see which buses have what degree of membership).
% 
vw = full(vw);
old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(mfilename('fullpath')),'third_party'),...
  fullfile(fileparts(mfilename('fullpath')),'third_party', 'matlab_bgl'));

if nargin > 4
  debug = true;
else
  debug = false;
end

if isstruct(opts)
  if ~isfield(opts,'card_min')
    card_min = 1;
  else
    card_min = opts.card_min;
  end
  if ~isfield(opts,'thrs')
    thrs = sqrt(2)/2;
  else
    thrs = opts.thrs;
  end 
  if ~isfield(opts,'improv')
    improv = true;
  else
    improv = opts.improv;
  end
elseif isscalar(opts)&&isnumeric(opts)&&opts>0
  card_min = opts;
  thrs = sqrt(3)/2;  % for back-compatibility
  improv = true;  % for back-compatibility
else
  error('[%s] The 4th input must be an options structure or a positive integer (minimum cluster size)', mfilename);
end
nc = size(Z,2);
 m = size(Z,1);
core_nodes = cell(nc,1);
eXp = ones(nc,1);
cut_offs = Inf(nc,1);
X = Utils.normalize_rows(Z);
mx_Z = max(Z,[],1);
[mx_Z,idx] = sort(mx_Z,'descend');
Z = Z(:,idx);
X = X(:,idx);

XX = X;
ZZ = Z;
for j = 1:1:nc
  if debug
    [core, phi, cut_off] = core_est(ZZ(:,j), XX(:,j), vw, adj, card_min, thrs, lbl);
  else
    [core, phi, cut_off] = core_est(ZZ(:,j), XX(:,j), vw, adj, card_min, thrs);
  end
  ZZ(core,:) = 0;
  XX(core,:) = 0;
  core_nodes{j} = core;
  eXp(j) = phi;
  cut_offs(j) = cut_off;
end
idx_del = eXp == 0;

% here and now (!) we assume that small cores have been filtered in advance
assert(~any(idx_del));  
core_nodes(idx_del) = [];
Z(:,idx_del) = [];
X(:,idx_del) = [];
cut_offs(idx_del) = [];
eXp(idx_del) = [];
nc = size(Z,2);

% Worst first
[eXp, idx] = sort(eXp, 'descend');
core_nodes = core_nodes(idx);   % "best order" (in terms of eXp)
cut_offs = cut_offs(idx); 
Z = Z(:,idx);
X = X(:,idx);

% Get integer adj (for max-flow-min-cut)
core_all = vertcat(core_nodes{:});
[i, j, w0] = find(adj);
int_infinity = 2^33;
mf = 1e6/max(w0);
while int_infinity>=2^32
  mf = mf/2;
  w = w0*mf;
  w(w<1.00000001) = 1;
  adj_int = sparse(i, j, w, m, m);
  int_infinity = sum(sum(adj_int))+2*sum(sum(adj_int(core_all,:)))+1;  % this should be larger than any conceivable flow...
end

% Repeat the cores, this time try to improve the initial ones...
if nc>2 && improv  % if nc=2, the subsequent uckwaycut() will be excellent 
  XX = X;
  ZZ = Z;
  for j = 1:1:nc
    core = core_nodes{j};
    others = my_setdiff(core_all, core);
    [core2, flow] = improve(adj_int, core, others, int_infinity, debug);
    vol_insd = sum(sum( adj(core2, core2) ));
    vol_full = sum(sum( adj(core2,:) ));
    wgt_full = sum(vw(core2));
    cut_curr = vol_full - vol_insd;
    if debug
      assert( abs(round(flow/mf) - round(cut_curr))/cut_curr<1e-3);
    end
    eXp2 = cut_curr/wgt_full;
    core2 = find(core2);
    
    if thrs>0.708
      if debug
        [core1, eXp1, cut_off] = core_est(ZZ(:,j), XX(:,j), vw, adj, card_min, 0.708, lbl);
      else
        [core1, eXp1, cut_off] = core_est(ZZ(:,j), XX(:,j), vw, adj, card_min, 0.708);
      end      
    else
      core1 = core_nodes{j};
      eXp1 = eXp(j);
    end
    
    if eXp1<eXp2
      core_nodes{j} = core1;
      eXp(j) = eXp1;
    else
      core_nodes{j} = core2;
      eXp(j) = eXp2;
    end
    
    ZZ(core_nodes{j},:) = 0;  % from now, not to be assigned to any other core
    XX(core_nodes{j},:) = 0;  % from now, not to be assigned to any other core
    core_all = vertcat(core_nodes{:});
  end
end
end


function [core, eXp, cut_off] = core_est(vecZ, vecX, vw, adj, card_min, thrs_low, lbl)
  persistent col;
  if isempty(col)
    col=1;
  else
    col=col+1;
  end
  m = numel(vecZ);
  if nargin>6
    debug=true;
  else
    debug=false;
  end
  clr = {'r','g','b'};
  
  [vecZ_lrge, idx_lrge] = sort(vecZ, 'descend');  
  vec_lrge = vecX(idx_lrge);
  % vecZ_lrge(1) is always included, even if vec_lrge(1)<0.707; 
  % This may in theory break the algorithm, but very unlikely;  
  assert(all(vec_lrge(1:card_min)>=thrs_low), 'The largest column entry is not well-aligned. The eigenvector alignment might be poor.')
  idx_curr = idx_lrge(1:card_min); 
  
  core_bin = false(1,m);    %logical, current indexing    
  core = NaN(1,m);          %integer, to store the additional steps
  num_save = numel(idx_curr);
  core_bin(idx_curr) = true;
  core(1:num_save) = idx_curr;  
  
  eXp_save = Inf(m,1);  
  vol_insd = sum(sum(adj(core_bin,core_bin)));
  vol_full = sum(sum(adj(core_bin,:)))+sum(sum(adj(:,core_bin)));
  wgt_full = sum(vw(core_bin));
  cut_curr = vol_full/2 - vol_insd;
  eXp_save(num_save) = cut_curr/wgt_full;   %if adj is K, and NCut is minimized    
  %eXp_save(num_save) = cut_curr;   %if adj is inv(M)*K instead of K, and "dynamic coupling" is minimized
  eXp_save(1:num_save-1) = linspace(1.1,1.025,num_save-1).*eXp_save(num_save);  %to discard those without messing up the plots
  
  idx_sweep = idx_lrge(card_min+1:m);
  vol_tril = sum(adj(idx_sweep,:), 2);  %out-degrees
  vol_triu = sum(adj(:,idx_sweep), 1)';  %in-degrees
  wgt_curr = vw(idx_sweep,:);  %wgt may not be the same as vol, thus two variables
  idx_seen = zeros(1,m);
  idx_seen(1:num_save) = 1:num_save;
  idx_skip = [];  % for dbg  
  for j = 1:1:numel(idx_sweep')
    if vecX(idx_sweep(j))<thrs_low
      if any(vecX(idx_sweep(j+1:end))>thrs_low)
        idx_skip(end+1) = idx_sweep(j);
      end      
      continue;
    end
    vol_insd = vol_insd + sum(adj(idx_sweep(j),core_bin)) + sum(adj(core_bin,idx_sweep(j)));    
    wgt_full = wgt_full + wgt_curr(j);        
    vol_full = vol_full + vol_tril(j) + vol_triu(j);
    cut_curr = vol_full/2 - vol_insd;
    num_save = num_save + 1;    
    eXp_save(num_save) = cut_curr/wgt_full;   %if adj is K, and NCut is minimized    
    %eXp_save(num_save) = cut_curr;   %if adj is inv(M)*K instead of K, and "dynamic coupling" is minimized
    core_bin(idx_sweep(j)) = true;   %update core0 for the next
    core(num_save) = idx_sweep(j);
    idx_seen(num_save) = card_min+j;
    if debug
      assert( abs(eXp_save(num_save) - (0.5*(sum(sum(adj(core_bin,:)))+sum(sum(adj(:,core_bin)))) - sum(sum(adj(core_bin,core_bin)))) / sum(vw(core_bin)) )<1e-10);  %if adj is K, and NCut is minimized       
      %assert( abs(eXp_save(num_save) - (0.5*(sum(sum(adj(core_bin,:)))+sum(sum(adj(:,core_bin)))) - sum(sum(adj(core_bin,core_bin)))) )<1e-10);  %if adj is inv(M)*K instead of K, and "dynamic coupling" is minimized
    end
  end    
  eXp_save = eXp_save(1:num_save);
  idx_seen = idx_seen(1:num_save);   
  
  [eXp_min, idx_core] = sort(eXp_save, 'ascend');
  iter = 1;
  iter_lim = 10;
  idx_rght = core;  % "above threshold"/legit nodes out of all considered  
  while iter <= iter_lim 
    core = idx_rght(1:1:idx_core(iter));
    adj_core = adj(core,core);
    [n_cc] = GraphUtils.alecconncomp(adj_core);
    if numel(core)<card_min   %case when initial core is smaller than card_min was excluded upstream once or even twice 
      iter = iter_lim+1; 
      n_cc = 2;
      break;
    end
    if n_cc==1
      eXp = eXp_save(idx_core(iter)); 
      idx_core = idx_core(iter);
      break;
    else
      iter = iter + 1;
    end
  end
  if iter>iter_lim && n_cc>1
    eXp = eXp_save(idx_core(1));
    idx_core = idx_core(1);
    core = idx_rght(1:1:idx_core);  % to later choose the largest connected component
  end
  
  if debug
    figure;
    subplot(3,1,1);
    plot(eXp_save, 'k.-', 'MarkerSize', 14);    
    hold on;
    plot(eXp_save(1:idx_core), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
    hold off;
    xl = xlim;
    % ylim([min(eXp_save)-0.001 min(min(eXp_save)+0.1, max(eXp_save))]);
    xlim(xl);
    % set(gca, 'Position', [0.1300    0.5838    0.6    0.2750])
    xlabel('Node number, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel(['$\phi(C_',num2str(col),')$, [-]'], 'Interpreter', 'Latex', 'FontSize', 13);
  end  
  if debug
    x = 1:1:num_save;
    subplot(3,1,2);
    plot(x, vec_lrge(idx_seen), 'k.-');
    xl = xlim;
    hold on;
    plot([xl(1),xl(2)], [vecX(core(end)), vecX(core(end))], 'k-.');
    plot(1:1:numel(core),vecX(core),'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
    for k = 1:num_save
        if iscell(lbl(idx_rght(k)))
          row_lbl = lbl{idx_rght(k)};
        else
          row_lbl = lbl(idx_rght(k));
        end
        text(k, vec_lrge(k), num2str(row_lbl), 'FontSize', 9, 'Color', 'k', 'FontWeight',...
            'bold', 'VerticalAlignment','top', 'HorizontalAlignment', 'center');
    end    
    hold off;
    % ylim([max(0.707,vec_lrge(num_save+1)) 1]);
    % set(gca, 'Position', [0.1300    0.1762    0.6    0.2750])
    xlabel('Node number, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel(['\boldmath$X^*_',num2str(col),'$, [-]'], 'Interpreter', 'Latex', 'FontSize', 12);
  end
  if debug
    subplot(3,1,3);
    plot(x, vecZ_lrge(idx_seen), 'k.-');
    xl = xlim;
    hold on;
    plot([xl(1),xl(2)], [vecZ(core(end)), vecZ(core(end))], 'k-.');
    plot(1:1:numel(core),sort(vecZ(core),'descend'),'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
    for k = 1:num_save
        if iscell(lbl(idx_rght(k)))
          row_lbl = lbl{idx_rght(k)};
        else
          row_lbl = lbl(idx_rght(k));
        end
        text(k, vecZ_lrge(k), num2str(row_lbl), 'FontSize', 9, 'Color', 'k', 'FontWeight',...
            'bold', 'VerticalAlignment','top', 'HorizontalAlignment', 'center');
    end    
    hold off;
    % ylim([max(0.707,vec_lrge(num_save+1)) 1]);
    % set(gca, 'Position', [0.1300    0.1762    0.6    0.2750])
    xlabel('Node number, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel(['\boldmath$Z^*_',num2str(col),'$, [-]'], 'Interpreter', 'Latex', 'FontSize', 12);
  end  
  
  % Here it could be simplified to leaving one connected component only
  adj_core = adj(core,core);
  [n_cc, idx_cc] = GraphUtils.alecconncomp(adj_core);    
  core = core(:);  
  while n_cc>1
    thrs = vecX(core(end));  % for curiosity...
    card_cc = tabulate(idx_cc);
    [~, card_ord] = sort(card_cc(:,2), 'descend');
    cc_pairs = [card_ord(1) card_ord(2)];
    vx2 = idx_cc==cc_pairs(2);
    vx2 = core(vx2);
    core = my_setdiff(core, vx2);
    eXp = ( sum(sum(adj(core,:))) - sum(sum(adj(core,core))) )/sum(vw(core));
    adj_core = adj(core,core);
    [n_cc, idx_cc] = GraphUtils.alecconncomp(adj_core);
  end  
  num_core = numel(core);  
  cut_off = vec_lrge(num_core);  
end


function [part_ind, flow] = improve(adj, core, others, int_infinity, debug)
m = size(adj, 1);

% Each flow problem add a fake source and sink as the m+1 and m+2 vertex
Aflow = adj;
Aflow(core,m+1) = int_infinity*ones(length(core),1);
Aflow(m+1,core) = int_infinity*ones(1,length(core)); % sparse(m+1,1);
Aflow(others,m+2) = int_infinity*ones(length(others),1);
Aflow(m+2,others) = int_infinity*ones(1,length(others)); % sparse(m+1,1);

% solve the max-flow problem
[flow, ci] = max_flow(Aflow, m+1, m+2);
if flow>0  % success with matlab-bgl
  part_ind = ci==sign(ci(m+1));  
else  % try the same with MATLAB (requires MATLAB newer than R2015b)
  [fr,to,wgt] = find(Aflow);
   Gflow = digraph(fr,to,wgt);
   [flow,~,cs,~] = maxflow(Gflow, m+1, m+2);
   part_ind = false(m+2,1);
   part_ind(cs) = true;
end
part_ind = part_ind(1:m);
if debug
  assert( flow == sum(sum(floor(adj(part_ind,~part_ind)))) );
end
end


