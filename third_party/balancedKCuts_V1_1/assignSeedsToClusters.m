function newClusters = assignSeedsToClusters(W, clusters, vertexWeights, seeds, currentClass, targetClasses, CheegerFlag, imbalance)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    
    volV = sum(vertexWeights);
    degree = sum(W,2);
    [~, ~, cuts, balances] = bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance);

    volumes = clusters'*vertexWeights;
    if size(balances,1) > 1
        balances = balances';
    end
    if CheegerFlag
        % dummy is matrix of size |seeds| x |targetClasses|: stores the
        % resulting volumes in location (i, j) where seed i is added to
        % targetClass j
        dummy = bsxfun(@plus, repmat(volumes(targetClasses)', length(seeds), 1), vertexWeights(seeds)); 
        newBalances = min( dummy*imbalance, volV - dummy );
        dummy = volumes(currentClass) - vertexWeights(seeds);
        newBalanceOfL = min(dummy*imbalance, volV - dummy);
    else            
        newBalances = bsxfun(@plus, repmat(balances(targetClasses), length(seeds), 1), vertexWeights(seeds)); 
        newBalanceOfL = balances(currentClass) - vertexWeights(seeds);
    end

    % Change in the mth ratio of balanced multi-cut
    deltaM = bsxfun(@plus, bsxfun(@plus, -2*(clusters(:, targetClasses)'*W(:, seeds))', cuts(targetClasses)'), degree(seeds))./newBalances ...
                                    - repmat(cuts(targetClasses)'./balances(targetClasses), length(seeds), 1) ;         
    % New value for lth ratio is computed instead of change because the
    % old lth ratio is constant for all the targetClasses.
    scoreMatrix = bsxfun(@plus, deltaM, (cuts(currentClass) - degree(seeds) + 2*(clusters(:, currentClass)'*W(:, seeds))')./(newBalanceOfL));
        
    [scores, bestAlternativeClasses] = min(scoreMatrix, [], 2);
    [~, maxIx] = max(scores); % higher cut means the seed should stay in the currentClass
    newClusters = clusters;
    for i=1:maxIx-1
        newClusters(seeds(i), currentClass) = 0;
        newClusters(seeds(i), targetClasses(bestAlternativeClasses(i))) = 1;
    end
    for i=maxIx+1:length(seeds)
        newClusters(seeds(i), currentClass) = 0;
        newClusters(seeds(i), targetClasses(bestAlternativeClasses(i))) = 1;
    end        

end