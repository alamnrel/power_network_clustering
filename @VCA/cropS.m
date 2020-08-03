function S = cropS(S, SET_FROM, SET_TO)
idx_col = ismember(S.buscol, SET_FROM);
idx_row = ismember(S.busrow, SET_TO);
S.S = S.S(idx_row, idx_col);
S.buscol = S.buscol(idx_col);
S.busrow = S.busrow(idx_row);
end