function [Fnew, lambdasNew, sgBalanceNew, degenerateFlag, convergenceFlag, nIters, degeneratedF, deltaPosNew, deltaNegNew, finalFNew, yDual, bestFNewInTermsOfCut] = ...
    innerProblem(W, F, sgBalance, PL, minBalance, maxBalance, lambdas, MAXITERS, iterFreq, simplexProjectionFlag, vertexWeights, CheegerFlag, imbalance, yDual, debugMode)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    if(~exist('debugMode', 'var')), debugMode = false; end    
    [n, k] = size(F);
    [ixPL, jxPL] = find(PL);
    offsetBalances = full(sum(sparse(ixPL, jxPL, vertexWeights(ixPL), n, k))');
    % if the number of labels is more than one per class, merge them.
    if max(sum(PL)) > 1
        Y = zeros(n, 1); % label vector for PL; it takes values from 1:k
        labels = 1:k;
        Y(ixPL) = labels(jxPL);
        [Wp, mapSparsify] = sparsify_vx(W, Y, debugMode);
        np = size(Wp, 1); 
        if np < k || sum(sum(Wp)) == 0
            if debugMode
                display('Corner cases... exiting the inner problem');
            end
            [Fnew, lambdasNew, sgBalanceNew, degenerateFlag, convergenceFlag, nIters, degeneratedF, deltaPosNew, deltaNegNew, finalFNew, yDual, bestFNewInTermsOfCut] ...
                = deal(F, lambdas, sgBalance, 0, 0, 0, F, 0, 0, F, yDual, F);
            return;
        end
        
        Fp(mapSparsify, :) = F; % compute F for the sparsified graph
        PLPrime(mapSparsify, :) = PL; % compute PL on the sparsified graph: now there is only one PL per class               
        vertexWeightsPrime = sparse(mapSparsify, 1, vertexWeights, np, 1);
        [Fnewp, lambdas1, degenerateFlag, convergenceFlag, nIters, degeneratedFp, deltaPosNew, deltaNegNew, finalFNewp, yDual, bestFNewInTermsOfCut] = ...
            PDHG(Wp, Fp, [], PLPrime, minBalance, maxBalance, lambdas, MAXITERS, iterFreq, simplexProjectionFlag, vertexWeightsPrime, CheegerFlag, imbalance, offsetBalances, yDual, debugMode); %sgBalance will be recomputed in this case by the function, PDHG
        % recover the values on the original graph
        Fnew = Fnewp(mapSparsify, :);
        degeneratedF = degeneratedFp;%(mapSparsify, :);
        finalFNew = finalFNewp(mapSparsify, :);
        bestFNewInTermsOfCut = bestFNewInTermsOfCut(mapSparsify, :);
        
    else
        [Fnew, lambdas1, degenerateFlag, convergenceFlag, nIters, degeneratedF, deltaPosNew, deltaNegNew, finalFNew, yDual, bestFNewInTermsOfCut] = PDHG(W, F, sgBalance, PL, minBalance, maxBalance, lambdas, MAXITERS, iterFreq, simplexProjectionFlag, vertexWeights, CheegerFlag, imbalance, offsetBalances, yDual, debugMode);
    end
    
    [lambdasNew, sgBalanceNew] = objBCutMultiClass(W, Fnew, vertexWeights, CheegerFlag, imbalance);
    if debugMode
        assert( norm(lambdasNew - lambdas1, 'inf') <= 1e-8 );
    end

end


function [Fnew, lambdasNew, degenerateFlag, convergenceFlag, nIters, degeneratedF, deltaPosNew, deltaNegNew, finalFnew, yDual, bestFNewInTermsOfCut] = ...
        PDHG(W, F, sgBalance, PL, minBalance, maxBalance, lambdas, MAXITERS, iterFreq, simplexProjectionFlag, vertexWeights, CheegerFlag, imbalance, offsetBalances, yDual, debugMode)

    [n, k] = size(F);
    Worig = W; 
    PLorig = PL;
    vertexWeightsOrig = vertexWeights;
    map = 1:n;
    [ixPL, jxPL] = find(PL); % Now only one Pseudo vertex per class (every thing else is merged into one)

    if sum(sum(PL)) <= 0
        degreePL = zeros(k, 1);
        cutEdgesPL = zeros(n, k);    
        theta = zeros(k,1); eta = zeros(k,1); thetaBar = zeros(k,1);
        clusters = threshold(F);
        initialBCut = sum(objBCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance));
    else    
        degree = sum(W,2); degreePL = zeros(k,1); degreePL(jxPL) = degree(ixPL);
        %map = setdiff(1:n, ixPL);
        map = ~ismember(1:n, ixPL);
        WMinusPL = W(map, map);
        FMinusPL = F(map, :);
        n = size(WMinusPL, 1);    
        jxPLComplement = find(~ismember(1:k, jxPL)); %setdiff(1:k, jxPL);
        if ~isempty(jxPLComplement)
            cutEdgesPL = sparse(n,k);
            cutEdgesPL(:, jxPLComplement) = repmat(-sum(W(map, ixPL), 2), 1, length(jxPLComplement));
        end
        cutEdgesPL(:, jxPL) = bsxfun(@minus, 2*W(map, ixPL), sum(W(map, ixPL), 2) );    
        
        if n < k || sum(sum(WMinusPL)) == 0
            if debugMode
                display('Corner cases... exiting PDHG');
            end
            [Fnew, lambdasNew, degenerateFlag, convergenceFlag, nIters, degeneratedF, deltaPosNew, deltaNegNew, finalFnew, yDual, bestFNewInTermsOfCut] ...
                = deal(F, lambdas, 0, 0, 0, F, 0, 0, F, yDual, F);
            return;
        end
        
        W = WMinusPL; F = FMinusPL;   % We solve the problem on the graph without PL (edges from them are hard-coded)
        vertexWeights = vertexWeights(map);        

        %theta = imbalance*nPLPerClass; eta = (k - 1 - imbalance)*nPLPerClass ;
        eta = min([imbalance*offsetBalances sum(vertexWeights)+sum(offsetBalances)-offsetBalances], [], 2);
        theta = imbalance*offsetBalances - eta; 
        thetaBar = sum(offsetBalances) - offsetBalances - eta; %sum(offsetBalances) - (imbalance+1)*offsetBalances;    
        clusters = threshold(F);
        initialBCut = sum(objBCutMultiClass(W, clusters, vertexWeights, CheegerFlag, imbalance, degreePL, cutEdgesPL, eta, theta, thetaBar)); 
        [lambdas1, sgBalance] = objBCutMultiClass(W, F, vertexWeights, CheegerFlag, imbalance, degreePL, cutEdgesPL, eta, theta, thetaBar); 
        % lambdas1 and lambdas passed should match!
        if debugMode
            assert( norm(lambdas - lambdas1, 'inf') <= 1e-8 );
        end
        PL = sparse(n, k); % Now PL have been hard-coded, hence restting it.
    end    
    
    cutEdgesPL = full(cutEdgesPL);
    degreePL = full(degreePL);
    PL = full(PL);
    F = full(F);
    sgBalance = full(sgBalance);
    gamma = 0; % this was used earlier for soft-enforcement of Pseudo Labels; now this is obselete.
    
    % Compute diagonal pre-conditioners
    m = minBalance; M = maxBalance;
    [ixt, jxt, wvalt] = find(triu(W));
    mhalf = full(sum(sum(W>0))/2);
    nDualVars = 2*k+n+2*mhalf*k + k*(~CheegerFlag); % if CheegerFlag is false, we can handle upper bound constraint on the balances
    Sigma = zeros(nDualVars, 1);
    counter = 1;
    % Descent constraint
    for l=1:k

        Sigma(counter) = 1/ (sum(wvalt)+sum(abs(-lambdas(l)*sgBalance(:,l) - cutEdgesPL(:,l) - gamma*PL(:,l)))+m+M);
        counter = counter+1;

    end

    % Lower bound constraint
    for l=1:k
        Sigma(counter) = 1/sum(abs(sgBalance(:,l)));
        counter = counter + 1;
    end

    % upper bound constraint
    if ~CheegerFlag    
        for l=1:k
            Sigma(counter+l-1) = 1/sum(abs(sgBalance(:,l)));
        end
    end

    % two sets of edges constraint
    index = 2*k + k*(~CheegerFlag);
    Sigma(index+1:index+mhalf*2*k) = 1/3; % three variables: alpha_ij >= |f_i - f_j|
    index = index+mhalf*2*k;
    % simplex constraint
    Sigma(index+1:index+n) = 1/k;

    nPrimalVars = mhalf*k+n*k+2*k;
    Tau = zeros(nPrimalVars, 1);

    for l=1:k
        Tau((l-1)*mhalf+1:l*mhalf) = 1./(wvalt+2);
    end

    nNeighbors = full(sum(W>0,2)); % from two sets of edge constraints
    for l=1:k
        dummy = (abs(-lambdas(l)*sgBalance(:,l) - cutEdgesPL(:,l) - gamma*PL(:,l)) + 1 + abs(sgBalance(:,l)) + abs(sgBalance(:,l))*(~CheegerFlag) + 2*nNeighbors);
        Tau(mhalf*k+(l-1)*n+1:mhalf*k+l*n) = 1./dummy;
    end

    Tau(end-2*k+1:end-k) = 1/m;
    Tau(end-k+1:end) = 1/M;

    % Primal initialization for the inner problem
    alpha = zeros(mhalf*k,1);
    for l=1:k    
        alpha( (l-1)*mhalf+1:l*mhalf ) = abs(F(ixt, l) - F(jxt, l));        
    end
    x = [alpha; F(:); zeros(2*k,1)];
    
    scale = min(max([wvalt+2; dummy; m; M]), 200);    
    if length(yDual) ~= nDualVars % previous iterate cannot be used because of the merging of Pseudo labels
        yDual = randn(nDualVars, 1);
    end
    [bestFDelta, finalFDelta, degenerateFlag, convergenceFlag, nIters, degeneratedF, yDual, bestFDeltaInTermsOfCut] = ...
        mexLpPDHG(triu(W), full(sgBalance(:)), m, M, full(lambdas), full(Tau), full(Sigma), full(x), full(yDual), n, mhalf, k, MAXITERS, imbalance, scale, full(cutEdgesPL(:)), full(degreePL), full(PL(:)), 0, iterFreq, full(eta), full(theta), full(thetaBar), full(vertexWeights), CheegerFlag==1, initialBCut, simplexProjectionFlag, debugMode);
    bestF = reshape(bestFDelta(1:n*k), n, k);
    bestFInTermsOfCut = reshape(bestFDeltaInTermsOfCut(1:n*k), n, k);
    finalF = reshape(finalFDelta(1:n*k), n, k);
    index = n*k;
    deltaPosNew = bestFDelta(index+1:index+k);
    index = index+k;
    deltaNegNew = bestFDelta(index+1:index+k);
    degeneratedF = reshape(degeneratedF, n, k);

    Fnew = zeros(size(Worig, 1),k);
    finalFnew = zeros(size(Worig, 1),k);
    bestFNewInTermsOfCut = zeros(size(Worig, 1),k);
    Fnew(map, :) = bestF; % bestF is on the graph without PL
    bestFNewInTermsOfCut(map, :) = bestFInTermsOfCut;
    finalFnew(map, :) = finalF;
    % Now update the Pseudovertices   
    if sum(sum(PLorig)) > 0
        Fnew = max(Fnew, PLorig);
        finalFnew = max(finalFnew, PLorig);
        bestFNewInTermsOfCut = max(bestFNewInTermsOfCut, PLorig);
    end
    
    Fnew = full(Fnew);
    finalFnew = full(finalFnew);
    bestFNewInTermsOfCut = full(bestFNewInTermsOfCut);
    lambdasNewp = objBCutMultiClass(W, bestF, vertexWeights, CheegerFlag, imbalance, degreePL, cutEdgesPL, eta, theta, thetaBar);
    lambdasNew = objBCutMultiClass(Worig, Fnew, vertexWeightsOrig, CheegerFlag, imbalance);

    if debugMode
        assert(abs(sum(lambdasNewp) - sum(lambdasNew)) <= 1e-8)
    end


end
