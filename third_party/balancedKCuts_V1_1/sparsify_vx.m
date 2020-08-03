function [Wp, map] = sparsify_vx(W, Y, debugMode)
% [Wp, map] = sparsify_vx(W, Y)
% Input: 
% W: Weight matrix 
% Y: a labeling of the vertices: same label (~=0) implies must-link. Label 0
%    means unlabeled vertex and they are not merged.
% Output:
% Wp: Weight matrix of the reduced graph
% map: the map from the original graph vertices to the vertices in the reduced graph
%
% To find the vertex weights of the reduced graph, use: sparse(map, 1, vertex_weights, size(Wp, 1), 1);
% To get the clustering of the reduced graph from that of the original
% graph, use: Yp(map) = Y; % the vector map has duplicate entires so the
% Y values of the mapped vertices are overwritten by the Y value of last vertex mapped
% To recover the original clustering from the clustering on the reduced
% graph, use: Y = Yp(map)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    if(~exist('debugMode', 'var')), debugMode = false; end
    n = size(W, 1);
    singleNodes = find(Y==0);
    nSingleNodes = length(singleNodes);
    labels = unique(Y(Y~=0));
    nSuperNodes = length(labels);
    np = nSingleNodes+nSuperNodes;
        
    % The reduced graph would have single nodes first and super nodes at the end
    map = zeros(np, 1);
    map(singleNodes) = 1:nSingleNodes;
    for i=1:length(labels)
        map(Y==labels(i)) = nSingleNodes+i;
    end
    
    Wp = sparse(np, np);    
    ixs_p = map(1:n); % we are using original weight matrix to sum the edges, so we need the map on the original graph    
    for i=1:nSuperNodes
        Wp(:, nSingleNodes+i) = sparse( ixs_p, 1, sum( W(:, Y==labels(i)), 2 ), np, 1); % It uses the fact that elements of vals with the duplicate (ix, jx) would be added!      
        % Here the column sum corresponds to merging edges from a super node and row sum corresponds to merging edges from different super nodes
    end    
    
    Wp = [W(singleNodes, singleNodes) Wp(1:nSingleNodes, nSingleNodes+1:end); Wp(:, nSingleNodes+1:end)'];       % Wp = [W(ix_single, ix_single) Wtemp(:, n_ixp+1:end); Wp(:, n_ixp+1:end)'];
    Wp(1:np+1:np*np) = 0;    %set the diagonal entries to zero 
    
    % Now sanity checks
    if debugMode
        Yp = unique([map, Y], 'rows'); Yp = Yp(:,2);
        assert( abs(tv(W, Y) - tv(Wp, Yp)) <= 1e-07 );
        assert( norm(Y - Yp(map)) == 0 );
    end
    
end