function Q = loewdin(C)
[U, ~, V] = svd(C);
Q = U*V';
end