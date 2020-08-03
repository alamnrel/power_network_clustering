function inc = adj2inc( adj )
%adj2inc This method takes an adjacency matrix and produces an incidence 
%  matrix out of it.
%
%  Author: Ilya Tyuryukanov
%  Date of first version: 30 May 2016
%  Last revision: 30 May 2016

assert(ismatrix(adj)&&isreal(adj)&&~any(diag(adj)),...
    '[%s] The adj-matrix should be a real matrix with zero diagonal', mfilename);
assert(isequal(logical(adj),logical(adj')),...
    '[%s] The adj-matrix should posess a symmetrical non-zero structure', mfilename);
n = nnz(adj)/2;
adj = tril(adj);
m = size(adj, 1);

[i_adj, j_adj] = find(adj);
lin_idx = [i_adj(:).'; j_adj(:).'];
i_inc = lin_idx(:);
j_inc = repmat(1:n, [2, 1]);
j_inc = j_inc(:);
v_inc = ones(1, 2*n);
inc = sparse(i_inc, j_inc, v_inc, m, n);

end