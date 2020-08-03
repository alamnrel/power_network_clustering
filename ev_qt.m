function Q = ev_qt(Vr)

if isa(Vr, 'cell')
  nc = numel(Vr);  
elseif isa(Vr, 'double')
  nc = 1;
  Vr = {Vr};
else
  error('');
end
Q = zeros(nc, 1);

% (Average cosine to the nearest axis)
m  = size(Vr{1}, 1);
for i = 1:1:nc
  V = Vr{i};
  maxsim = max(V, [], 2);
  Q(i) = sum(maxsim)/m;
end

end


%{
% Percent above threshold
P = nnz(abs(V)>thrs)/size(V,1);
X = discretize(V);

Q(1) = norm(X - Vd, 'fro')/m;
%Q(2) = sum(sum( abs(X - V) ));
%cossim = V*X';
%Q(3) = sum(diag(cossim))/m;
%}