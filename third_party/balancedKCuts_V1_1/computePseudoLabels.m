function [ixPL, clusterSizes] = computePseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance)
% ixPL are the indices of Pseudovertices
% ixPL(i, j) = t means the vertex t is a pseudo vertex with label is j
% Ordering between Pseudovertices: t1 = ixPL(i, j), t2 = ixPL(l, j) with i < l implies t1 > t2 (t1 is fixed before t2)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    nClasses = size(clusters, 2);
    [~, ~, ~, vertexOrderingPerClass] = computeVertexOrdering(W, clusters, vertexWeights, CheegerFlag, imbalance);                    
    clusterSizes = histc(vertexOrderingPerClass(:,2), 1:nClasses);
    minNPL = min(clusterSizes);
    ixPL = zeros(minNPL, nClasses);
    for l=1:nClasses
        dummy = vertexOrderingPerClass(vertexOrderingPerClass(:,2)==l, 1);      
        ixPL(:, l) = dummy(1:minNPL);
    end
    
end