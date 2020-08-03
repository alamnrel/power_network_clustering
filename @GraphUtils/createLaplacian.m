function L = createLaplacian(adj, nm, vw)
%
% Return graph Laplacian from input (graph adjancency matrix, Laplacian
% type, vertex weights).
%
% Author: Ilya Tyuryukanov
% Date of first version: 18 August 2015
% Revision: 5 November 2015 (??)
% Revision: 30 March 2019 (added custom vertex weights as 3rd input and 
%   generalized the Laplacian computation accordingly).
% 

m = size(adj, 2);
d = sum(adj, 2);   % weighted degrees of vertices
if nargin<3
  vw = d;
end

if issparse(adj)
  D = spdiags(d, 0, m, m);
else
  D = diag(d);
end
L = D - adj;  %the classical graph Laplacian
    
% Compute normalized Laplacians 
switch nm
  case 'rwk'
    M = 1./vw;  
    if issparse(L)
      L = spdiags (M, 0, m, m) * L;  
    else
      L = bsxfun(@times, L, M); 
    end
  case 'sym'
    M = 1./sqrt(vw);
    [i, j, ~] = find(L);
    den = sparse(i, j, bsxfun(@times, M(i), M(j)), m, m);
    L = L.*den;  
    if ~issymmetric(L)
      L = (L+L.')/2;  % enforce symmetry to speed up eig() among other things
    end
  case 'off'
    return;
  otherwise
    error([mfilename,':WrongPublicInput'],...
      ['[%s] Invalid Laplacian normalisation method in GraphUtils.createLaplacian(g, nm). ',...
      'Use nm = ''off''/''rwk''/''sym''.'], mfilename);
end
end