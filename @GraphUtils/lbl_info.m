function  [eXp, ixs, pcut] = lbl_info(adj, vw, T)
%
% This is similar to cc_info, but operating not on connected components,
% but on partitioning labels. That is, if partitions encoded in T are
% connected, the eXp and ixs outputs of cc_info and lbl_info will be the
% same except the column ordering. The objective to have the column
% ordering of eXp and ixs strictly enforced by the order of column labels
% in T is the main motivation and goal of this function (otherwise cc_info
% would be enough).
%

lbl = unique(T);
lbl = lbl(:)';
nc = numel(lbl);
eXp = zeros(6,nc);
ixs = cell(1,nc);
pcut = zeros(1,nc);
vol_t = sum(sum(adj,2));
for i = 1:1:nc
  cur = lbl(i);    
  idx = T==cur;  
  ixs{i} = find(idx);
  vol_i = sum(sum(adj(idx,idx)));  % internal volume
  vol_e = sum(vw(idx));            % total cluster node weight  
  vol_f = sum(sum(adj(:,idx)));  % "full" volume
  bnd_i = vol_f - vol_i;
  pcut(i) = bnd_i;
  eXp(1, i) = abs(bnd_i/vol_e);  % expansion (NCut) of ith connected component
  eXp(2, i) = abs(bnd_i/nnz(idx));  % RCut of ith connected component
  eXp(3, i) = nnz(idx);  % number of vertices in ith connected component
  eXp(4, i) = vol_i/vol_t - ((bnd_i+vol_i)/vol_t)^2; % similarity-based modularity
  eXp(5, i) = bnd_i;  % cut of ith connected component
  eXp(6, i) = vol_i;  % internal volume of ith connected component
end

end