function adj_out = simadj2distadj( adj, mode )
% simadj2distadj This function takes (presumably) adjacency matrix of a 
% similarity graph adj and converts this graph to a distance graph
% (adjacency matrix adj_out) by using one of tranformations which are
% selected by the 'mode' input parameter.
%
% Author: Ilya Tyuryukanov
% Date: 30 May 2015
% 25 November 2015: ??
% 28 June 2016: Undirected graph is not enforced anymore via tril(). I.e. 
%   symmetric adj in --> symmetric adj_out (this is more logical & tril()
%   may be applied afterwards).

switch mode
  case 'inv1'
    adj_out = spfun(@(x) 1./x, adj);
  case 'inv2'  % according to Large-Scale Spectral Clustering on Graphs by Liu, Wang, Danilevsky
    maxA = max(max(adj));
    tmp = unique(adj(:), 'sorted');
    minD = -log10(tmp(end-1)/maxA);
    adj_out = spfun(@(x) -log10(x/maxA) + minD/1e3, adj);  % convert tight coupling (large weight) to "short distance"
  case 'orig'
    adj_out = adj; 
  otherwise
    error([mfilename,':WrongPublicInput'],...
      ['[%s] Invalid 2nd input in GraphUtils.simadj2distadj(adj, mode). ',...
      'Use mode = ''inv1''/''inv2''/''orig''.'], mfilename);
end


