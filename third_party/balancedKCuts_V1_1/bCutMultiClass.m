function [bcuts, clustering, cuts, balances] = bCutMultiClass(W, F, vertex_weights, Cheeger_balance, imbalance)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    if ~exist('imbalance', 'var')
        imbalance = 1;
    end
    
    k = size(F, 2); % if F was the matrix nxk, then k should not change according to number of unique labels of clusters!
    if k > 1
        [~, clustering] = max(F,[],2);
        classes = 1:k;
    else
        clustering = F;
        classes = unique(clustering);
        k = length(classes);
    end

    bcuts = zeros(k,1);    
    cuts = bcuts; 
    balances = cuts;
    for i=1:k, 
        cuts(i) = sum(sum(W(clustering==classes(i),clustering~=classes(i))));
        if Cheeger_balance
            balances(i) = min(imbalance*sum(vertex_weights(clustering==classes(i))), sum(vertex_weights(clustering~=classes(i))));            
        else
            balances(i) = sum(vertex_weights(clustering==classes(i)));            
        end
        bcuts(i) = cuts(i)/balances(i); 
    end
    
end
