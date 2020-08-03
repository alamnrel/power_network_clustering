function [LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = vcs_sc(S_GL, S_LL, param, ADJ, debug)

if iscell(S_GL.S)
  python = true;
  S_GL.S = DiGSI.unpack_mat(S_GL.S);
else
  python = false;
end
GEN = S_GL.buscol;
LOD = S_GL.busrow;
nb = numel(GEN) + numel(LOD);
bus = [GEN(:); LOD(:)];
igen = 1:1:numel(GEN);
ilod = numel(GEN)+1:1:nb;
epsilon = eps('double');

if nargin<3
  param.num_clu = 2:min(20,nb);
elseif iscell(param.num_clu)   % Python list/not a scalar
  param.num_clu = cell2mat(param.num_clu);
  param.num_clu(param.num_clu>nb) = [];
end
if nargin<5
  debug = false;
end

% Sensitivity
A_SC = zeros(nb, nb);
A_SC(igen,ilod) = S_GL.S';
A_SC(ilod,igen) = S_GL.S;
if ~isempty(S_LL)  
  if python
    S_LL.S = DiGSI.unpack_mat(S_LL.S);
  end  
  if ~issymmetric(S_LL.S)
    S_LL.S = (S_LL.S + S_LL.S')/2;
  end
  assert( isequal(S_LL.busrow, S_LL.buscol), ...
    [mfilename,':WrongInput'], '[%s] The load-to-load ',...
     'sensitivities should be between the same sets of buses.', mfilename );
  [~, ~, ib] = intersect(S_GL.busrow, S_LL.busrow, 'stable');  
  assert( isequal(S_LL.busrow(ib), S_GL.busrow) );
  S_LL.S = S_LL.S(ib,ib);   % rearrange S_LL to fit S_GL
  A_SC(ilod, ilod) = S_LL.S;
end
A_SC(1:nb+1:nb*nb) = 0;
A_SC(A_SC<epsilon) = 0;   % negative entries should normally not be present..

% Neglect minor connected components
A_SC = A_SC*1000;   % scale A_SC up (better numerics?)
[nc, cc] = GraphUtils.alecconncomp(A_SC);
pp_cc = tabulate(cc);
[sz_cc_max, ix_cc_max] = max(pp_cc(:,2));
cc_dom = pp_cc(ix_cc_max, 1);
A_SC = A_SC(cc == cc_dom, cc == cc_dom);
ibus = [igen, ilod];
ibus = ibus(cc == cc_dom);
busX = bus(cc==cc_dom);
m = size(A_SC,1);

% Create graph and apply spectral clustering
k_try = param.num_clu(1:end);
G_SC = PFgraph('adj', A_SC, 'bus', bus(cc == cc_dom));  % bus included to monitor which stations go where
L_SC = GraphUtils.createLaplacian(A_SC, 'sym');
L_SC(1:m+1:m*m) = 0;  L_SC = -L_SC;   % An for > speed and accuracy
[V0, v0] = GraphUtils.spClust(L_SC, 'la', max(k_try));
[VV, YS, tt] = GraphUtils.nc_ev_sym(V0, [1 1], [0 0 0 1], [], min(k_try), false);
VV = VV(k_try-min(k_try)+1);
YS = YS(k_try-min(k_try)+1);
CQ = ev_qt(VV);
%CQ = YS;
if debug && numel(CQ)>1
  figure;
  x_width = 1.5*241; y_width = 1.5*161;
  plot(1:numel(CQ), CQ, 'k.-', 'MarkerSize', 18, 'LineWidth', 1);
  xlabel('Nr. of eigenvectors, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  ylabel('$C$, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
  set(gcf, 'pos', [100 100 x_width y_width]);
  set(gca, 'XTick', 1:1:numel(CQ));
  xlim([1,numel(CQ)]);
  set(gca, 'XTickLabel', cellfun(@num2str, num2cell(2:1:numel(CQ)+1), 'unif', false));
  set(gca, 'xgrid', 'on');  
end

% Do SC (for vcs-clustering just choose the best cost among available)
[cq_srt, ord_try] = sort(CQ, 'descend');
nc = ord_try(1);
V = VV{nc};
if nargin<4 || isempty(ADJ)
  ADJ = G_SC.adj;
  lbl = 1:1:size(V,1);  % V directly stems from G_SC.adj -> direct correspondence
else
  lbl = G_SC.bus;
end
[core_nodes] = quick_cores(V, lbl, ADJ, 1, false);
% core_nodes = spec_cores(V, G_SC.adj, 1);
T_SC = G_SC.uckwaycut(core_nodes, 2, 1);
%[val, T_SC] = max(V, [], 2);  % option number 2
% T_SC = NaN(1,size(V,1));      % option number 3
% thrs = 0.866;
% for j = 1:1:size(V,2)
  % T_SC(V(:,j)>thrs) = j;
% end
T = NaN(nb,1);
T(ibus) = T_SC;

% Characterise load-to-load partition through expansion
 [val, T_SC] = max(V, [], 2);  % to reconsiliate possible NaNs from above
 T_SC = G_SC.contig(T_SC);
 vw = sum(A_SC,2);
[eXp, ixs, cut_tot] = GraphUtils.lbl_info(A_SC, vw, T_SC);
eXp(4,:) = cut_tot(:)';  % modularity is not needed

% Return buses and labels:
grp = unique(T_SC);
LBL = cell(numel(grp), 1);
bus = bus(:);
ind_gen = ismember(bus, GEN(:));
ind_lod = ismember(bus, LOD(:));
for k = 1:1:numel(grp)
  LBL{k} = bus(T==grp(k));
  LBL_GEN{k} = bus(T==grp(k) & ind_gen);
  LBL_LOD{k} = bus(T==grp(k) & ind_lod);
end
irst = isnan(T);
RST = bus(irst);
if python
  CQ = num2cell(CQ(:), 2);
  eXp = cellfun(@(x) num2cell(x,2), num2cell(eXp,1), 'un', 0);
end
% CQ = YS;
end

