function coh_grp = coh_loose(Vs, tol)
% Identify "loosely" coherent generators based on:
%   Vs - slow eigensubspace (electromechanical eigenvectors used for slow-cohrency)
%   tol - grouping tolerance
%
% References:
% [1] J.H. Chow, "New algorithms for slow coherency aggregation of large
%   power systems," in Systems and Control Theory for Power Systems, IMA
%   volumes in mathematics and its applications, vol. 64, ed. by J.H. Chow,
%   R.J. Thomas, P.V. Kokotovic? (Springer-Verlag, New York, 1994)
% [2] J.H. Chow and K.W. Cheung, ?A toolbox for power system dynamics and 
%   control engineering education and research,? IEEE Trans. Power Syst., 
%   vol. 7, no. 4, pp. 1559-1564, Nov. 1992.
%
% Ilya Tyuryukanov, 2016

n_gen = size(Vs, 1);

% Normalize the columns and rows of Vs to length 1
Vs(:,1) = ones(1,n_gen);  % 1st column is a constant real eigenvector
Vs = normalize_rows(Vs')';
Vs = normalize_rows(Vs);

% Matrix of cosines of angles between row vectors of Vs
dotprod = Vs * Vs';
c_matr = dotprod - tol*ones(n_gen, n_gen);

% Find (loosely) coherent generators
coh_grp = [1:1:n_gen; 1:1:n_gen];  % in the beginning each generator
                                   % belongs to "its own" coherent group
for i = 2:1:n_gen
  for j = 1:1:i-1
    if c_matr(i,j) > 0 
      if coh_grp(2,i) ~= coh_grp(2,j)
        grp = min(coh_grp(2,i), coh_grp(2,j));
        coh_grp(2,i) = grp;
        coh_grp(2,j) = grp;
      end
    end    
  end
end
[coh_grp(2,:), idx] = sort(coh_grp(2,:), 'ascend');
coh_grp(1,:) = coh_grp(1,idx);
offset = diff(coh_grp(2,:));
offset = offset - 1; 
offset(offset < 0) = 0;
offset = [offset(1), offset];
offset = cumsum(offset);
coh_grp(2,:) = coh_grp(2,:) - offset;

end

function V = normalize_rows(V)
V = bsxfun(@rdivide, V, sqrt(sum(V.^2,2)));
end