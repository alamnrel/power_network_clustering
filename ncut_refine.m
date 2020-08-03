function [T, eXp_curr] = ncut_refine(adj, vw, T, debug)
                              
assert(all(all(adj))>=0);
assert(~any(diag(adj)));
assert( all(diff(unique(T))==1) && min(T)==1 );
debug = true;

T = T(:);
nc = max(T);
m = size(adj,1);
[eXp] = GraphUtils.lbl_info(adj, vw, T);
ncut0 = mean(eXp(1,:));

% Current state
vol_insd = zeros(1,nc);
vol_full = zeros(1,nc);
wgt_full = zeros(1,nc);
eXp_curr = zeros(1,nc);
grp_card = zeros(1,nc);
for k = 1:1:nc
  core_bin = T==k;
  vol_insd(k) = sum(sum(adj(core_bin,core_bin)));
  vol_full(k) = sum(sum(adj(core_bin,:)))+sum(sum(adj(:,core_bin)));
  wgt_full(k) = sum(vw(core_bin));
  cut_curr    = vol_full(k)/2 - vol_insd(k);
  eXp_curr(k) = cut_curr/wgt_full(k);  
  grp_card(k) = nnz(core_bin);
end

% Binarize cores for faster access
cores_bin = false(m,nc);
for k = 1:1:nc
  cores_bin(T==k,k) = true;
end

% Perform (sequential) label propagation
vol_tril = sum(adj, 2);   %out-degrees
vol_triu = sum(adj, 1)';  %in-degrees
max_all = Inf;
while max_all>0
  eXp_rcv = zeros(m,nc);
  eXp_rem = zeros(m,nc); 
  for j = 1:1:m
    j_adj = adj(j,:);
    adj_j = adj(:,j);
    for k = 1:1:nc
      if k~=T(j) && grp_card(T(j))>1
        %Receiving group k:
        corN = cores_bin(:,k);
        insN = vol_insd(k) + sum(j_adj(corN)) + sum(adj_j(corN));
        wgtN = wgt_full(k) + vw(j);
        volN = vol_full(k) + vol_tril(j) + vol_triu(j);
        cutN = volN/2 - insN;
        eXpN = cutN/wgtN;
        eXp_rcv(j,k) = eXp_curr(k)-eXpN;
        %Sending group T(j):
        corO = cores_bin(:,T(j));
        insO = vol_insd(T(j)) - sum(j_adj(corO)) - sum(adj_j(corO));
        wgtO = wgt_full(T(j)) - vw(j);
        volO = vol_full(T(j)) - vol_tril(j) - vol_triu(j);
        cutO = volO/2 - insO;
        eXpO = cutO/wgtO;
        % All other groups keep their wgt_full, vol_full and vol_insd
        eXp_rem(j,k) = eXp_curr(T(j))-eXpO;
      end
    end
  end
  [max_grp,rcv_grp] = max(eXp_rcv+eXp_rem,[],2);
  [max_all,mov_nod] = max(max_grp);
  if max_all>0
    rcv_grp = rcv_grp(mov_nod);  %
    rem_grp = T(mov_nod);
    
    vol_insd(rcv_grp) = vol_insd(rcv_grp) + sum(adj(mov_nod,cores_bin(:,rcv_grp))) + sum(adj(cores_bin(:,rcv_grp),mov_nod));
    vol_full(rcv_grp) = vol_full(rcv_grp) + vol_tril(mov_nod) + vol_triu(mov_nod);
    wgt_full(rcv_grp) = wgt_full(rcv_grp) + vw(mov_nod);
    grp_card(rcv_grp) = grp_card(rcv_grp) + 1;
    cut_curr          = vol_full(rcv_grp)/2 - vol_insd(rcv_grp);
    assert( abs(cut_curr/wgt_full(rcv_grp)-eXp_curr(rcv_grp)+eXp_rcv(mov_nod,rcv_grp))<1e-6 );
    eXp_curr(rcv_grp) = cut_curr/wgt_full(rcv_grp);
    
    vol_insd(rem_grp) = vol_insd(rem_grp) - sum(adj(mov_nod,cores_bin(:,rem_grp))) - sum(adj(cores_bin(:,rem_grp),mov_nod));
    vol_full(rem_grp) = vol_full(rem_grp) - vol_tril(mov_nod) - vol_triu(mov_nod);
    wgt_full(rem_grp) = wgt_full(rem_grp) - vw(mov_nod);
    grp_card(rem_grp) = grp_card(rem_grp) - 1;    
    cut_curr          = vol_full(rem_grp)/2 - vol_insd(rem_grp);
    assert( abs(cut_curr/wgt_full(rem_grp)-eXp_curr(rem_grp)+eXp_rem(mov_nod,rcv_grp))<1e-6 );
    eXp_curr(rem_grp) = cut_curr/wgt_full(rem_grp);    
    
    cores_bin(mov_nod,rcv_grp) = true;
    cores_bin(mov_nod,rem_grp) = false;
    T(mov_nod) = rcv_grp;
  end
end

if debug
  for k = 1:1:nc
    assert( isequal(cores_bin,full(sparse(1:m,T,true,m,nc))) );
    assert( grp_card(k)==nnz(cores_bin(:,k)) );
    assert( abs(eXp_curr(k) - (0.5*(sum(sum(adj(cores_bin(:,k),:)))+sum(sum(adj(:,cores_bin(:,k))))) - sum(sum(adj(cores_bin(:,k),cores_bin(:,k))))) / sum(vw(cores_bin(:,k))))<1e-10 );    
  end  
  [eXpX] = GraphUtils.lbl_info(adj, vw, T);
  ncutX = mean(eXpX(1,:));  
  assert(ncutX<=ncut0);
end

