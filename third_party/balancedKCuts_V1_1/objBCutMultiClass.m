function [lambdas, sgBalance, mtvs, balances] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, beta, degreePL, cutEdgesPL, eta, theta, thetaBar)
% evaluates the set function and computes its subgradient at F: eta + min{beta vol(C) + theta, vol(Cbar) + thetaBar}, for all C in the partition
% Conditions on theta and thetaBar to make the set function vanish on empty set: theta = 0 or vol(V) + thetaBar = 0
% The last five arguments are optional    
% Note that the subgradient for the set function is same as the one without the offset: eta
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    
    [n, k] = size(F);    
    if(~exist('degreePL', 'var')), degreePL = 0; end
    if(~exist('cutEdgesPL', 'var')), cutEdgesPL = 0; end    
    if(~exist('eta', 'var')), eta = zeros(k,1); end
    if(~exist('theta', 'var')), theta = zeros(k,1); end
    if(~exist('thetaBar', 'var')), thetaBar = zeros(k,1); end
    if(~exist('debugMode', 'var')), debugMode = false; end
        
    if CheegerFlag
        volV = sum(vertexWeights);
        	
        medianThreshold = (beta*volV + theta - thetaBar)/(beta+1);            
        [~, ixSorted] = sort(F, 1);
        vertexWeightsVec = vertexWeights(:);
        vertexWeightsSorted = reshape(vertexWeightsVec(ixSorted(:)), n, k);
        cumSumSortedVolumes = cumsum(vertexWeightsSorted);
                
        ixPlus  = cumSumSortedVolumes >= bsxfun(@plus, vertexWeightsSorted, medianThreshold');
        [~, jxTemp, valTemp] = find(ixSorted.*ixPlus);
        sgBalance = sparse(valTemp, jxTemp, beta*vertexWeightsSorted(ixPlus), n, k);
        
        ixMinus = cumSumSortedVolumes <= bsxfun(@plus, zeros(n,k), medianThreshold');
        [~, jxTemp, valTemp] = find(ixSorted.*ixMinus);
        sgBalance = sgBalance + sparse(valTemp, jxTemp, -vertexWeightsSorted(ixMinus), n, k);                
        
        ixZero = (ixPlus|ixMinus) == 0;        
        
		if sum(sum(ixZero))>0
			[~, jxTemp, valTemp] = find(ixSorted.*ixZero);
    	    
        	repetitions = sum(ixZero)';
	        [nonZeroIx, ~, nonZeroRepetitions] = find(repetitions);
    	    repeatIx = zeros(sum(nonZeroRepetitions), 1);
       	 	repeatIx(cumsum([1; nonZeroRepetitions(1:end-1)])) = 1; % starting indices for new value
	        repeatIx = cumsum(repeatIx); % repeat the indices                
    	    thetaRep = theta(nonZeroIx(repeatIx));
        	thetaBarRep = thetaBar(nonZeroIx(repeatIx));
        
	        sgBalance = sgBalance + sparse(valTemp, jxTemp, (beta+1)*cumSumSortedVolumes(ixZero) - beta*volV -vertexWeightsSorted(ixZero) + thetaBarRep(:) - thetaRep(:), n, k);       
		end
        
        % Sanity Check: 
        if debugMode
            vertexWeightsRep = repmat(vertexWeights, 1, k);
            if length(unique(F))==2   
                g = bsxfun(@gt, F, min(F)); 
            else
                g = bsxfun(@gt, F, mean(F));
            end

            % Note that the subgradient for the set function is same as the one without the offset: eta
            assert( max(abs( sum(sgBalance.*g) - min(beta*sum(vertexWeightsRep.*g)+theta', volV - sum(vertexWeightsRep.*g) + thetaBar') ))  <= 1e-6 ); 
        end        
    else
        sgBalance = repmat(vertexWeights, 1, k);
    end
    
    mtvs = mtv(W, F) + degreePL - sum(cutEdgesPL.*F)';
    balances = sum(sgBalance.*F) + eta';
    lambdas = mtvs./balances';
    
end