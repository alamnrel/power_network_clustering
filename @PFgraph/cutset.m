function cutind = cutset(g, busind)
%
% This method calculates the branch index vector (into the incidence 
% matrix) given the bus index vectors of the islands (busind) and the 
% incidence matrix (g.inc).
%
% Author: Ilya Tyuryukanov
% Date of first version: 28 September 2015
% Last revision: 28 September 2015

inc_matr = g.inc;
cutind = GraphUtils.cutset(inc_matr, busind);
end

