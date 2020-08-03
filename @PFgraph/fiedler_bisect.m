function busind = fiedler_bisect(g)
%
% Bisect graph based on its Fiedler vector

if isempty(g.Loff)
  g = g.createLaplacian('off');
end

[X, lmbda] = eig(full(g.Loff));  % full eig is more stable...
lmbda = diag(lmbda);
[~, idx] = sort(lmbda, 'ascend');
X = real(X(:,idx));   
X = X(:,2);  % only 2nd eigenvector is important, real is "just in case"

if ~all(X==median(X))
  busind(1,:) = X >  median(X);
  busind(2,:) = X <= median(X);
else
  middle = ceil(numel(X)*0.5);
  busind = false(2,numel(X));
  busind(1, 1:middle) = true;
  busind(2, middle+1:end) = true;
end
end