function [feasibleF, permutation] = getFeasibleClusters(W, vertexWeights, F, imbalance, Labels, debugMode)        
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    if(~exist('debugMode', 'var')), debugMode = false; end    
    % Change F so that it has one given label (seed) per class
    [seeds, labeledClasses] = find(Labels);
    seedsPerCluster = sum(F(seeds, :));
    richClusters = find(seedsPerCluster > 1);
    while ~isempty(richClusters)
       for i=length(richClusters)
           seedsInThisCluster = seeds(F(seeds, richClusters(i))>0);
           if debugMode, assert(length(seedsInThisCluster) > 1); end
           orphanClusters = find(seedsPerCluster == 0);
           if debugMode, assert(~isempty(orphanClusters)); end
           % Keep the highest ordered seed in the current cluster and move the rest
           % to the best alternative cluster.
           F = assignSeedsToClusters(W, F, vertexWeights, seedsInThisCluster, richClusters(i), orphanClusters, 1, imbalance);
           seedsPerCluster = sum(F(seeds, :));           
       end
       richClusters = find(seedsPerCluster > 1); % each time all the rich clusters of the current iteration are neutralized; an orphan might become rich but will definitely have at least one less seed than the previous rich
    end
    % The following assertion fails if at least one class does not have a
    % label
    %if debugMode, assert( unique(sum(F(seeds, :)))==1 ); end

    % Reorder F so that it satisfies the label constraints:
    % clusterings maybe same but the labeling can be different
    feasibleF = F;
    [~, clustering] = max(F, [],2);
    permutation = clustering(seeds);
    nClasses = size(F,2);
    for l=1:length(labeledClasses)
        feasibleF(:, labeledClasses(l)) = F(:, permutation(l));
    end
    remColumns = find(~ismember(1:nClasses, permutation));
    unLabeledClasses = find(~ismember(1:nClasses, labeledClasses));
    for l=1:length(remColumns)
        feasibleF(:, unLabeledClasses(l)) = F(:, remColumns(l));
    end
    
    if debugMode
        %assert( norm(feasibleF(seeds, :)-eye(nClasses), 'fro') == 0 );
        [~, feasibleClustering] = max(feasibleF, [],2);
        assert(sum(feasibleClustering(seeds) ~= labeledClasses) == 0, 'Labels constraints are not satisfied')
        assert( abs( sum(bCutMultiClass(W, feasibleF, vertexWeights, 1, imbalance)) - sum(bCutMultiClass(W, F, vertexWeights, 1, imbalance)) ) < 1e-8 ); 
    end

end