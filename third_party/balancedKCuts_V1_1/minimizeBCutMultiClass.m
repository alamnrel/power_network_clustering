function [clustering, bcuts, clusters, F, lambdas] = minimizeBCutMultiClass(W, vertexWeights, nClasses, CheegerFlag, imbalance, F, initialNPseudoLabelsPerClass, minBalance, maxBalance, MAXITERS, iterFreq, simplexProjectionFlag, Labels, debugMode)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    if(~exist('debugMode', 'var')), debugMode = false; end
    yDual = 0;    
    nVertices = size(W,1);
    [LabeledVertices, ~] = find(Labels); % Labels are from transductive setting
    F(LabeledVertices, :) = 0;
    F = max(F, Labels);
    [lambdas, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance);
    clusters = threshold(F);       
    bcuts = bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance);
	if isnan(sum(bcuts))
        if debugMode
            display('Changing starting F to avoid NaN bcuts ... ');
        end
		[~,nanClustering] = max(F,[],2);
		nanIndices = setdiff(1:nClasses,unique(nanClustering));
		for nanIdx = nanIndices
			[~,ixmax] = max(F(:,nanIdx));
			F(ixmax,(F(ixmax,:)==1)) = 1-2*eps;
			F(ixmax,nanIdx)=1;
		end
		[lambdas, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance);
        clusters = threshold(F);       
		bcuts = bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance);
        if isnan(sum(bcuts)) % probably W is disconnected and eigvecs have to be modified.
            F = randn(nVertices, nClasses);
            F = F + abs(min(min(F))); % To avoid negative values in F and to keep the random distribution intact
            F = F/max(max(F));
            [lambdas, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance);
            clusters = threshold(F);       
            bcuts = bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance);
        end
	end
    %display(['Start : ', 'lambda = ', num2str(sum(lambdas), '%10.6f'), 'bcut = ', num2str(sum(bcuts), '%10.6f'), ', Ratios = ', num2str(bcuts', '%10.6f')]);
    if debugMode, display(['Start : ', 'lambda = ', num2str(sum(lambdas), '%10.6f'), ' bcut = ', num2str(sum(bcuts), '%10.6f')]); end
    
    %bestCutChanged = false;
    iter = 1; 
    CONVERGENCE = false; 
    finalBCut = inf; % used for restarting the iteration with discrete F
    currentNPseudoLabelsPerClass = initialNPseudoLabelsPerClass;   
    PL = sparse(nVertices, nClasses);                    
    degenerateFlag = 0;       
    bestBCutSoFar = sum(bcuts); % this is for tracking the best cut found so far
    bestFNewInTermsOfCut = clusters;        
    outerIteration = true;
    while ~CONVERGENCE
            
        betterLambdaFound = true;
        allPseudoLabelsEnforced = false;
        
        while betterLambdaFound || ~allPseudoLabelsEnforced
            
            % change this to get back to old  PL strategy: PL = sparse(nVertices, nClasses);            
            if currentNPseudoLabelsPerClass > 0
                if allPseudoLabelsEnforced
                    if debugMode
                        display('All Pseudo labels enforced... exiting the inner loop');
                    end
                    break; 
                end
                
                allPseudoLabelsEnforced = currentNPseudoLabelsPerClass >= min(clusterSizes); % one class is full of Pseudo lables; so exit in the next iteration.
                [~, jx, PseudoVertices] = find(ixPL(1:currentNPseudoLabelsPerClass, :));
                PseudoVerticesComplement = ~ismember(1:nVertices, PseudoVertices);
                if nVertices < nClasses*currentNPseudoLabelsPerClass + nClasses || sum(sum(W(PseudoVerticesComplement, PseudoVerticesComplement))) == 0
                    if debugMode
                        display('Corner cases... exiting the inner loop');
                    end
                    break;
                end
                
                PL = sparse(PseudoVertices, jx, 1, nVertices, nClasses);
                F(PseudoVertices, :) = 0;
                F = max(F, PL); % F is always in [0,1]
                clustersDummy = threshold(F);            
                bcutsDummy = bCutMultiClass(W, clustersDummy, vertexWeights, CheegerFlag, imbalance);                  
                [lambdas, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance);                                    
                if debugMode 
                    disp(['bcut (after thresholding with PL): ', num2str(sum(bcutsDummy)), ' lambda (after setting PL): ', num2str(sum(lambdas)), ' rel diff: ', num2str(abs(sum(bcutsDummy)-sum(lambdas))/sum(lambdas))]);                
                end
                if ~outerIteration && abs(sum(bcutsDummy)-sum(lambdas))/sum(lambdas) < 1e-4 && ~degenerateFlag %&& ~bestCutChanged
                    % first condition is for the case when we restrat with F=clusters
                    % second condition implies all (possibly except for PL
                    % vertices) are 0/1. The outer iteration will anyway
                    % check this case, so we can exit the inner loop
                    % third condition ensures that the last inner problem
                    % is properly solved without degeneration
                    if debugMode, display('stopping label enforcement as the continuous solution saturated'); end
                    break;
                else                    
                    if debugMode, display(['enforcing ', num2str(currentNPseudoLabelsPerClass), ' labels per class']); end
                    if sum(lambdas) > 2*sum(bcuts) % to cut-short the saturation of F
                        F = clusters;
                        lambdas = bcuts;
                    end
                end
            end                                                          
                                        
            if debugMode, display(['starting lambda: ', num2str(sum(lambdas)), ' and starting BCut: ', num2str(sum(bcuts))]); end
            timeForInnerProblem = tic;
            [FNew, lambdasNew, sgBalanceNew, degenerateFlag, ~, ~, ~, ~, ~, ~, yDual, bestFNewInTermsOfCutTemp] = ...
                innerProblem(W, F, sgBalance, double((PL+Labels)>0), minBalance, maxBalance, lambdas, MAXITERS, iterFreq, simplexProjectionFlag, vertexWeights, CheegerFlag, imbalance, yDual, debugMode);            
            timeForInnerProblem = toc(timeForInnerProblem);
            if debugMode, display(['Time for inner prlbem: ', num2str(timeForInnerProblem/60), ' minutes']); end
            
            tempBCut = sum(bCutMultiClass(W, bestFNewInTermsOfCutTemp, vertexWeights, CheegerFlag, imbalance));
            if bestBCutSoFar > tempBCut
                bestBCutSoFar = tempBCut;
                bestFNewInTermsOfCut = bestFNewInTermsOfCutTemp;        
            end
            
            clustersNew = threshold(FNew);            
            bcutsNew = bCutMultiClass(W, clustersNew, vertexWeights, CheegerFlag, imbalance);    
            nVerticesMoved = sum(sum(clusters~=clustersNew))/2;                        
                        
            betterBCutFound = false;
            if sum(bcutsNew) < sum(bcuts)
                clusters = clustersNew;
                bcuts = bcutsNew;      
                betterBCutFound = true;
                % change this to get back to old  PL strategy: currentNPseudoLabelsPerClass = initialNPseudoLabelsPerClass;                
                allPseudoLabelsEnforced = false;
                
                if currentNPseudoLabelsPerClass > initialNPseudoLabelsPerClass
                    %[ixPL1, clusterSizes1] = getPseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);      
                    [ixPL, clusterSizes] = computePseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);      
                    % no need to check if currentNPseudoLabelsPerClass goes beyond the limits: each cluster has a size at least equal to currentNPseudoLabelsPerClass
                    % and this increases only when we change the cut and at that time we anyway check for the limits.
                    if debugMode, display(['got the PL from the following cut: ', num2str(sum(bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance))) ' cluster sizes: ' num2str(clusterSizes')]); end
                end
            end            
            
            if sum(lambdasNew) < sum(lambdas) && abs(sum(lambdasNew) - sum(lambdas))/sum(lambdas) > 1e-4 % upto some numerical accuracy
                betterLambdaFound = true;
                F = FNew;
                lambdas = lambdasNew;
                sgBalance = sgBalanceNew;
            else
                betterLambdaFound = false;
            end
            
            if debugMode
                display(['iter ', num2str(iter), ', lambda = ', num2str(sum(lambdas), '%10.6f'), ', bcut = ', num2str(sum(bcuts), '%10.6f'), ', bestBcutSoFar = ', num2str(bestBCutSoFar, '%10.6f'), ...
                ', number of vertices moved to a different class : ', num2str(nVerticesMoved)]);
            else
                %display(['iter ', num2str(iter), ', lambda = ', num2str(sum(lambdas), '%10.6f'), ', bcut = ', num2str(sum(bcuts), '%10.6f')]);
               % ', number of vertices moved to a different class : ', num2str(nVerticesMoved)]);
            end
            iter = iter+1;            

            if (~betterLambdaFound || degenerateFlag) && ~betterBCutFound 
                if currentNPseudoLabelsPerClass == initialNPseudoLabelsPerClass
                    %[ixPL1, clusterSizes1] = getPseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);      
                    [ixPL, clusterSizes] = computePseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);      
                    if debugMode, display(['got the PL from the following cut: ', num2str(sum(bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance))) ' cluster sizes: ' num2str(clusterSizes')]); end
                    currentNPseudoLabelsPerClass = min([clusterSizes; max([2*initialNPseudoLabelsPerClass, 1])]);                    
                else                    
                    currentNPseudoLabelsPerClass = min([clusterSizes; currentNPseudoLabelsPerClass*2]); % to ensure nPL <= smallest cluster size                                            
                end
            end
            outerIteration = false;
                       
        end
        
        outerIteration = true;
        if debugMode
            if debugMode, display('----------------- Outer Iteration -----------------'); end
        end
        
        bestBCutsSoFar = bCutMultiClass(W, bestFNewInTermsOfCut, vertexWeights, CheegerFlag, imbalance);
        %bestCutChanged = true;
        if sum(bestBCutsSoFar) < sum(bcuts)
            F = bestFNewInTermsOfCut;
            clusters = threshold(F);
            %[ixPL1, clusterSizes1] = getPseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);           
            [ixPL, clusterSizes] = computePseudoLabels(W, clusters, vertexWeights, CheegerFlag, imbalance);      
            currentNPseudoLabelsPerClass = min([clusterSizes; currentNPseudoLabelsPerClass]); % to ensure nPL <= smallest cluster size (this is needed as the bcut is replaced by bestBcut) 
            if debugMode, display(['got the PL from the following cut: ', num2str(sum(bCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance))) ' cluster sizes: ' num2str(clusterSizes')]); end
            bcuts = bestBCutsSoFar; % lambdas will be computed below when we replace F by clusters.            
            %bestCutChanged = true;
        end

        if finalBCut <= sum(bcuts)
            if debugMode, display('No descent in the outer iteration... exiting...'); end
            break;
        end        
        
        % change this to get to PL strategy: currentNPseudoLabelsPerClass = initialNPseudoLabelsPerClass;
        if abs(sum(bcuts)-sum(lambdas))/sum(lambdas) < 1e-4 && finalBCut <= sum(bcuts) 
            CONVERGENCE = true;
        else
            F = clusters;
            [lambdas, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance); % bcuts stay the same.
            finalBCut = sum(bcuts);            
        end               
                       
    end    
    
    if debugMode
        display(['Number of Pseudo labels: ' num2str(currentNPseudoLabelsPerClass)]);
    end
    [~, clustering] = max(clusters,[],2);
        
end