function [clustering, bCut, clusters] = balancedKCut(W, k, imbalance, vertexWeights, CheegerFlag, nInitializations, parForLoop, debugMode)
%
% balancedKCut minimizes balanced k-cut using the method described in  
% S. S. Rangapuram, P. K. Mudrakarta and M. Hein. 
% Tight Continuous relaxation of the Balanced k-Cut Problem. NIPS 2014.
% Cite the above paper if you are using this code.
%
% Usage: [clustering, bCut, clusters] = balancedKCut(W, k, imbalance, vertexWeights, CheegerFlag, nInitializations, parForLoop, debugMode)
%
% Inputs:
%   W:                  Symmetric (SPARSE) weight matrix with zero diagonal.
%   k:                  number of clusters.
%
% Optional Inputs:
%   imbalance:          imbalance paramter. Defalut value = k-1; for
%                       Normalized/Ratio Cheeger cut, use imbalance = 1.
%   vertexWeights:      vector of all ones for all types of Ratio cuts; i.e., balancing function is based on cardinalities (this is the defualt setting),
%                       degree vector for Normalized cut; one can also pass any generic vertex weights (strictly positive) instead of degree. 
%   CheegerFlag:        1 for Cheeger Cuts, 0 for Normalized/Ratio cut.
%                       In the latter case, imbalance must be equal to 1.
%                       Default value = 1.
%   nInitializations:   number of starting points for the method. 
%                       Default value = 12 (7 spectral clustering based initializations and 5 random initializations).
%                       For large datasets random initializations might take longer; 
%                       in this case set nInitializations to any value between 1
%                       and 7 (possibly compromising on quality).
%   parForLoop:         uses matlab parfor loop to to start the method from
%                       different starting points in parallel (this is default setting).
%                       Set this to false, if you do not want to use parfor.
%   debugMode:          prints additional debugging messages. Default value = false
%
% Outputs:
%   clustering:         result vector with labels 1, ..., k
%   bCut:               the balanced cut value
%   clusters:           cluster indicator matrix of size nVertices x k
%
% Use,
%   clustering = balancedKCut(W, k), for the asymmetric ratio Cheeger cut problem
%   clustering = balancedKCut(W, k, 1), for the ratio Cheeger cut problem
%   clustering = balancedKCut(W, k, 1, sum(W,2)), for the normalized Cheeger cut problem
%   clustering = balancedKCut(W, k, 1, ones(size(W,1),1), 0), for the ratio cut problem
%   clustering = balancedKCut(W, k, 1, sum(W,2), 0), for the normalized cut problem
%   clustering = balancedKCut(W, k, k-1, sum(W,2)), for the asymmetric normalized Cheeger cut problem
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    addpath('Ncut_9');
    addpath('GraphDemos/');    
    nVertices = size(W, 1);
    
    assert(sum(diag(W))==0, 'No self loops: diagonal of the weight matrix W must be zero.');
    assert(issparse(W) && isnumeric(W), 'The weight matrix W must be sparse and numeric');
    assert(issymmetric(W), 'The weight matrix W must be symmetric');
    assert(k >= 2, 'Wrong usage: number of clusters should at least be 2');
    if ~exist('imbalance', 'var'), imbalance = k-1; else assert(imbalance > 0, 'Wrong usgae: imbalance should be greater than zero'); end
    if ~exist('vertexWeights', 'var')
        vertexWeights = ones(nVertices, 1); 
    else
        assert(size(vertexWeights, 1) == nVertices & size(vertexWeights, 2) == 1, 'Wrong usgae: vertexWeights should be a vector of size nVertices x 1 (column vector)'); 
        assert(min(vertexWeights) > 0, 'Wrong usgae: vertexWeights should be strictly positive');
    end
    if ~exist('CheegerFlag', 'var'), CheegerFlag = 1; else assert(CheegerFlag == 1 || CheegerFlag ==0, 'CheegerFlag should either be zero or one'); end            
    if ~exist('nInitializations', 'var'), nInitializations = 12; else assert(nInitializations >= 1, 'Wrong usgae: number of starting points for the method should at least be 1'); end    
    if ~exist('parForLoop', 'var'), parForLoop = true; end            
    if ~exist('debugMode', 'var'), debugMode = false; end            
    
    if imbalance > 1
        assert(CheegerFlag == 1, 'Wrong usage: imbalance must be 1 for non-Cheeger balanced cuts')
        minBalance = imbalance*min(vertexWeights);
        maxBalance = imbalance*sum(vertexWeights)/k;        
    else
        minBalance = min(vertexWeights);
        maxBalance = sum(vertexWeights)/2;
    end
    
    nSpecStarts = min(2, nInitializations);
    nDiffusedStarts = min(5, ceil(0.5*(nInitializations-2)));
    nRandStarts = nInitializations-nSpecStarts-nDiffusedStarts;
        
    display('Computing Spectral clustering solution');    
    [ourSpecClusters, ~, ourSpecEigVecs, ~] = GD_SpectralClustering(W, k, vertexWeights); 
    display(['Spectral clustering solution: ' num2str(sum(bCutMultiClass(W, ourSpecClusters, vertexWeights, CheegerFlag, imbalance)))]);
    ourSpecF = zeros(nVertices, k); for i=1:k, ourSpecF(ourSpecClusters==i, i) = 1; end
    
    startingFs = zeros(nVertices, k, nInitializations);
    counter = 1;    
    startingFs(:, :, counter) = ourSpecF;    
    
    if nInitializations > 1
        ourSpecEigVecs = ourSpecEigVecs(:,2:k+1);
        ourSpecEigVecs = ourSpecEigVecs + abs(min(min(ourSpecEigVecs))); % To avoid negative values in F and to keep the random distribution intact
        ourSpecEigVecs = ourSpecEigVecs/max(max(ourSpecEigVecs));
        
        counter = counter+1;
        startingFs(:, :, counter) = ourSpecEigVecs;
        
        if nInitializations > 2
            Ncut9Clusters = ncutW(W,k); % nSeed = 0
            Ncut9SpecF = Ncut9Clusters; %zeros(nVertices, k); for i=1:k, Ncut9SpecF(Ncut9Clusters==i, i) = 1; end            
        end
        
        for diffusedStartIdx=1:nDiffusedStarts
            counter = counter+1;
            startingFs(:, :, counter) = getDiffusedStart(W,1,Ncut9SpecF,1,sparse(nVertices,k));                        
        end
               
        for randStartIdx=1:nRandStarts
            tempF = randn(nVertices,k);
            tempF = tempF + abs(min(min(tempF))); % To avoid negative values in F and to keep the random distribution intact
            tempF = tempF/max(max(tempF));
            counter = counter+1;
            startingFs(:, :, counter) = tempF;            
        end
                
    end
    
    if sum(isnan(ourSpecClusters)) > 0 % GD_spectral clustering failed, so use Ncut9
        if nInitializations <= 2 % Ncut9 is not computed yet
            Ncut9Clusters = ncutW(W,k); % nSeed = 0
            Ncut9SpecF = Ncut9Clusters; %zeros(nVertices, k); for i=1:k, Ncut9SpecF(Ncut9Clusters==i, i) = 1; end            
        end
        startingFs(:, :, 1) = Ncut9SpecF; 
    end
    
    bCuts_all = zeros(k, nInitializations);
    clusterings = zeros(nVertices, nInitializations);
    clusters_all = zeros(nVertices, k, nInitializations);
    
    if parForLoop && nInitializations <= 1
        parForLoop = false;
    end
    
    if parForLoop
        parfor i=1:nInitializations
            [clusterings(:, i), bCuts_all(:, i), clusters_all(:, :, i)] = minimizeBCutMultiClass(W, vertexWeights, k, 1, imbalance, startingFs(:, :, i), 0, minBalance, maxBalance, 3200, 400, 0, sparse(nVertices, k), debugMode);                
            if ~CheegerFlag
                bCuts_all(:, i) = bCutMultiClass(W, clusterings(:, i), vertexWeights, 0, 1);
            end
            display(['Solution obtained with the initialization ' num2str(i) ': '  num2str(sum(bCuts_all(:, i)))]);
        end 
    else
        for i=1:nInitializations
            [clusterings(:, i), bCuts_all(:, i), clusters_all(:, :, i)] = minimizeBCutMultiClass(W, vertexWeights, k, 1, imbalance, startingFs(:, :, i), 0, minBalance, maxBalance, 3200, 400, 0, sparse(nVertices, k), debugMode);                
            if ~CheegerFlag
                bCuts_all(:, i) = bCutMultiClass(W, clusterings(:, i), vertexWeights, 0, 1);
            end
            display(['Solution obtained with the initialization ' num2str(i) ': '  num2str(sum(bCuts_all(:, i)))]);
        end 
    end
    
    bCuts = sum(bCuts_all);
    [bCut, minIx] = min(bCuts);
    clustering = clusterings(:, minIx);
    clusters = clusters_all(:, :, minIx);
        
    display(['************* Final solution (unconstrained): ' num2str(bCut) ' ******************']);
    
    if debugMode
        assert(length(unique(clustering))==k, 'clustering does not have k labels');
        assert(length(clustering)==nVertices, 'clustering does not have n vertices');
        assert( abs(sum(bCutMultiClass(W, clustering, vertexWeights, CheegerFlag, imbalance)) - bCut) < 1e-8);        
    end
       
end