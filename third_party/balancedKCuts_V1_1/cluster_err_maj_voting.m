function [err, labels_predicted, misses] = cluster_err_maj_voting(truth, prediction)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    clusters_predicted = unique(prediction);
    k = length(clusters_predicted);
    labels_predicted = inf*ones(k, 1);
    misses = inf*ones(k,1);
    
    for i=1:k
        % for each predicted cluster assign a label according to the truth:
        % predicted cluster i receives the marjority ground truth label in that cluster
                
        truth_cluster_i = truth(prediction == clusters_predicted(i)); % it finds the truth on ith cluster
        labels_predicted(i) = mode(truth_cluster_i);  % majority label
        misses(i) = sum(truth_cluster_i ~= labels_predicted(i)); 
        
    end

    err = sum(misses)/length(prediction);

end