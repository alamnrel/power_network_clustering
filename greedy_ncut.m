function [T, eXp_curr] = greedy_ncut(adj, vw, cores, rest, debug)

assert(all(all(adj))>=0);
assert(~any(diag(adj)));
debug = true;

nc = numel(cores);
m = size(adj,1);
vol_insd = zeros(1,nc);
vol_full = zeros(1,nc);
wgt_full = zeros(1,nc);
eXp_curr = zeros(1,nc);

for k = 1:1:nc
  vol_insd(k) = sum(sum(adj(cores{k},cores{k})));
  vol_full(k) = sum(sum(adj(cores{k},:)))+sum(sum(adj(:,cores{k})));
  wgt_full(k) = sum(vw(cores{k}));
  cut_curr    = vol_full(k)/2 - vol_insd(k);
  eXp_curr(k) = cut_curr/wgt_full(k);  
end

%Binarize cores for faster access:
cores_bin = false(m,nc);
for k = 1:1:nc
  cores_bin(cores{k},k) = true;
end

T = zeros(m,1);
if isempty(rest)
  for k = 1:1:nc
    T(cores_bin(:,k)) = k;
  end
  return;
end

vol_tril = sum(adj, 2);   %out-degrees
vol_triu = sum(adj, 1)';  %in-degrees
while ~isempty(rest)  
  nr = numel(rest);
  eXp_redu = zeros(nr,nc);
  for j = 1:1:nr
    vx = rest(j);
    vx_adj = adj(vx,:);
    adj_vx = adj(:,vx);
    for k = 1:1:nc
      core = cores_bin(:,k);
      ins = vol_insd(k) + sum(vx_adj(core)) + sum(adj_vx(core));
      wgt = wgt_full(k) + vw(vx);
      vol = vol_full(k) + vol_tril(vx) + vol_triu(vx);
      cut = vol/2 - ins;
      eXp = cut/wgt;
      eXp_redu(j,k) = eXp_curr(k) - eXp;
    end
  end
  [max_rst,mov_cor] = max(eXp_redu,[],2);
  [max_all,idx_rst] = max(max_rst);
  mov_cor = mov_cor(idx_rst);  
  mov_rst = rest(idx_rst);
  vol_insd(mov_cor) = vol_insd(mov_cor) + sum(adj(mov_rst,cores_bin(:,mov_cor))) + sum(adj(cores_bin(:,mov_cor),mov_rst));
  vol_full(mov_cor) = vol_full(mov_cor) + vol_tril(mov_rst) + vol_triu(mov_rst);
  wgt_full(mov_cor) = wgt_full(mov_cor) + vw(mov_rst);
  cut_curr          = vol_full(mov_cor)/2 - vol_insd(mov_cor);
  assert( abs(cut_curr/wgt_full(mov_cor)-eXp_curr(mov_cor)+eXp_redu(idx_rst,mov_cor))<1e-6 ); 
  eXp_curr(mov_cor) = cut_curr/wgt_full(mov_cor);
  cores_bin(mov_rst,mov_cor) = true;
  cores{mov_cor} = [cores{mov_cor}; mov_rst];
  rest(idx_rst)  = [];
end

if debug
  for k = 1:1:nc
    assert(~any(cores_bin(setdiff(1:m,cores{k}),k)));
    assert(all(cores_bin(cores{k},k)));
    assert( abs(eXp_curr(k) - (0.5*(sum(sum(adj(cores{k},:)))+sum(sum(adj(:,cores{k})))) - sum(sum(adj(cores{k},cores{k})))) / sum(vw(cores{k})))<1e-10 );    
  end  
end

for k = 1:1:nc
  T(cores_bin(:,k)) = k;
end



