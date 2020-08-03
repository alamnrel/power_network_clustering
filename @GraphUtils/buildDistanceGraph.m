function [adjD, dm, nnInd, dst_srtd] = buildDistanceGraph(Y, adj, nm)
%
% Create graph with interconnections defined by adjacency matrix adj and 
% weights calculated from point cloud Y in Euclidian space with certain 
% geometric distance functions.
%
% Author: Ilya Tyuryukanov
% Date of first version: 18 November 2015
% Last revision: 18 November 2015

m = size(adj, 1);
[i, j, ~] = find(adj);

epsilon = eps('double');
if strcmp(nm, 'sym')
    dotprod = sum(Y(i,:).*Y(j,:), 2);
    sp_dst = real(acos(dotprod));  % spherical distances graph    
else
    sp_dst = sqrt(sum((Y(i,:)-Y(j,:)).^2, 2));
end
sp_dst(sp_dst<50*epsilon) = 50*epsilon;
adjD = sparse(i, j, sp_dst, m, m);  

if nargout > 1
    dm = GraphUtils.allShortestPaths(adjD, 'orig');
end

if nargout > 2
    [dst_srtd, nnInd] = sort(dm, 2);
    k = max(sum(logical(adj),2)) + 1; % max number of node's nearest neighbours including self   
    nnInd = nnInd(:, 1:k);
end

end