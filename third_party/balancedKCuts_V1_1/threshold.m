function clusters = threshold(F)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    [n, k] = size(F);
    [dummyx, clustering] = max(F, [], 2);
    clusters = zeros(n, k); for i=1:k, clusters(clustering==i, i) = 1; end
end
