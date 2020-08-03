function [LBL, T] = vca_rest(adj, bus, LBL)
% 
% Attempts to assign unassigned buses to voltage control zones based on
% graph connectivity between assigned buses (non-NaN labels in T) and
% unassigned buses (NaN labels in T).
% If an unassigned bus is adjacent to several zones (several labels in
% T), it will remain labeled as NaN.
% 

m = size(adj, 1);
T =  NaN(1, m);
LBL = cellfun(@(x)x(:), LBL, 'un', false);
for i = 1:1:numel(LBL)
  T(ismember(bus, LBL{i})) = i;
end
if ~any(isnan(T))
  return;
end

% Assign generator nodes to the adjacent buses:
while any(isnan(T))
  i_rst = find(isnan(T));  % nonassigned buses
  n_rst = numel(i_rst);
  for i = 1:1:n_rst
    nbr_cur = find( adj(:,i_rst(i)) );
    lbl_cur = T(nbr_cur); %#ok<FNDSB>
    lbl_cur(isnan(lbl_cur)) = [];
    if isempty(lbl_cur)
      continue;
    end
    if all(lbl_cur==lbl_cur(1))
      T(i_rst(i)) = lbl_cur(1);
      LBL{lbl_cur(1)} = [LBL{lbl_cur(1)}; bus(i_rst(i))];
    else
      T(i_rst(i)) = Inf;   % seen, but unassigned
    end
  end
end
T(isinf(T)) = NaN;
end