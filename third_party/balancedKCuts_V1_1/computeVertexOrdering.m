function [vertexOrdering, scoresOrdered, bestAlternativeClassesOrdered, vertexOrderingPerClass] = computeVertexOrdering(W, clusters, vertexWeights, CheegerFlag, imbalance)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    
    [nVertices, nClasses] = size(clusters);
    classes = 1:nClasses;
    scoreMatrix = zeros(nVertices, nClasses-1);
    volV = sum(vertexWeights);    
    degree = sum(W,2);
    [~, ~, cuts, balances] = bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance);
    if size(balances,1) > 1
        balances = balances';
    end
    if CheegerFlag  % balances are not same as volumes in this case
        volumes = clusters'*vertexWeights;
    end

    for l=1:nClasses
        ixInClassL = find( clusters(:,l)==1 );
        %display(['cluster ', num2str(l), ': ', num2str(length(ixInClassL))]);
        remainingClasses = classes(~ismember(classes, l));
        if CheegerFlag            
            % dummy is matrix of size |ixInClassL| x |remainingClasses|: stores the
            % resulting volumes in location (i, j) where vertex i is added to
            % class j
            dummy = bsxfun(@plus, repmat(volumes(remainingClasses)', length(ixInClassL), 1), vertexWeights(ixInClassL)); 
            newBalances = min( dummy*imbalance, volV - dummy );
            dummy = volumes(l) - vertexWeights(ixInClassL);
            newBalanceOfL = min(dummy*imbalance, volV - dummy);
        else            
            newBalances = bsxfun(@plus, repmat(balances(remainingClasses), length(ixInClassL), 1), vertexWeights(ixInClassL)); 
            newBalanceOfL = balances(l) - vertexWeights(ixInClassL);
        end
        % Change in the mth ratio of balanced multi-cut
        deltaM = bsxfun(@plus, bsxfun(@plus, -2*(clusters(:, remainingClasses)'*W(:, ixInClassL))', cuts(remainingClasses)'), degree(ixInClassL))./newBalances ...
                                        - repmat(cuts(remainingClasses)'./balances(remainingClasses), length(ixInClassL), 1) ;         
                                    
        % New value for lth ratio is computed instead of change because the
        % old lth ratio is constant for all the remainingClasses.                                    
        scoreMatrix(ixInClassL, :)  = bsxfun(@plus, deltaM, (cuts(l) - degree(ixInClassL) + 2*(clusters(:, l)'*W(:, ixInClassL))')./(newBalanceOfL));
    end

    [scores, bestAlternativeClasses] = min(scoreMatrix, [], 2);     
    [scoresOrdered, vertexOrdering] = sort(scores, 'descend');
    bestAlternativeClassesOrdered = bestAlternativeClasses(vertexOrdering);
    
    [~, clustering] = max(clusters, [], 2);
    vertexOrderingPerClass = [vertexOrdering, clustering(vertexOrdering)];

end