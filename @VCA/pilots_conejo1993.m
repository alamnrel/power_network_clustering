function [plt_set, plt_qlt] = pilots_conejo1993(S_GL, S_LL, C_LL, Q_x, plt_ini, plt_max, debug)
  
% 
% Returns multiple sets of pilot buses, each next set having one more bus
% than the previous ones, and the corresponding quality of voltage control
% metric for each set. This metric is the basis of the algorithm and it is 
% from the paper [1] that is the basis of this script.
% 
% [1] A. Conejo, T. Gomez, and J. de la Fuente, "Pilot-bus selection for 
% secondary voltage control," International Transactions on Electrical 
% Energy Systems, vol. 3, no. 5, pp. 359-366, 1993
% [2] J.L. Sancha , J.L. Fernandez, A. Cortes, J.T. Abarca, "Secondary 
% voltage control: analysis, solutions and simulation results for the 
% Spanish transmission system," IEEE Trans. on Power Systems, vol. 11,
% no. 2, pp. 630-638, May 1996.
%

if iscell(S_GL.S)
  python = true;
  S_GL.S = DiGSI.unpack_mat(S_GL.S);
  S_GL.buscol = DiGSI.unpack_mat(S_GL.buscol);
  S_GL.busrow = DiGSI.unpack_mat(S_GL.busrow);
else
  python = false;
end
if python
  S_LL.S = DiGSI.unpack_mat(S_LL.S);
  S_LL.busrow = DiGSI.unpack_mat(S_LL.busrow);
  S_LL.buscol = DiGSI.unpack_mat(S_LL.buscol);
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
SET_PQG = S_GL.buscol;
SET_LOD = S_LL.buscol;
num_lod = numel(SET_LOD);
if nargin<5
  plt_ini = [];
else
  plt_ini = plt_ini(:);
end
if nargin<6 || ~isscalar(plt_max) || ~isnumeric(plt_max)
  plt_max = num_lod-1;
end
if nargin<7
  debug = false;
end
B = S_GL.S;
M = S_LL.S;

% Conejo-matrices not multiplied by C*...*C':
P_L = M*C_LL*M';

% Search for pilots
plt_qlt = [];
plt_cur = plt_ini; 
plt_set = cell(1, plt_max);
SET_CND = setdiff(SET_LOD, plt_cur);
for k = 1:1:plt_max
  J_CAND = zeros(1, numel(SET_CND));
  for j = 1:1:numel(SET_CND)
    C = form_C([plt_cur; SET_CND(j)], num_lod);
    J_CAND(j) = calc_J(B, P_L, Q_x, C);
  end
  [j_min, idx_min] = min(J_CAND);
  plt_qlt = [plt_qlt, j_min];
  plt_cur = [plt_cur; SET_CND(idx_min)];
  plt_set{k} = plt_cur;  
  if debug
    plot(J_CAND);
    set(gca, 'XTick', 1:numel(SET_CND))
    set(gca, 'XTickLabel', cellfun(@num2str,num2cell(SET_CND'), 'un', 0));
    set(gca, 'xgrid', 'on');
    xlim([min(SET_CND), numel(SET_CND)]);
    ylim([min(J_CAND), 2*min(J_CAND)])
  end
  SET_CND = setdiff(SET_CND, SET_CND(idx_min));
end

end


function C = form_C(idx_plt, num_lod)
% Forms binary pilot bus indicator matrix.

C = zeros(numel(idx_plt), num_lod);
for i = 1:1:numel(idx_plt)
  C(i, idx_plt(i)) = 1;
end
end


function J = calc_J(B, P_L, Q_x, C)
  X = C*B;
  F = X'/(X*X');
  n_lod = size(Q_x, 2);
  Y = eye(n_lod) - B*F*C;
  J = trace( Q_x*Y*P_L*Y' );  % see [2] for the derivation
end




