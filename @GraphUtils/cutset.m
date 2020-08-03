function cutind = cutset(inc_matr, busind)
% 
% This method calculates the branch index vector (into the incidence 
% matrix) given the bus index vectors of the islands (busind) and the 
% incidence matrix (g.inc).
% 
% Author: Ilya Tyuryukanov
% Date of first version: 28 September 2015
% Last revision: 28 September 2015
% 

Nc = size(busind, 1);
n = size(inc_matr, 2);
cutind = zeros(Nc, n);  % indicator of cutset for each group of buses (cluster) in busInd
for i = 1:Nc
  j = double(busind(i,:));
  cut_i = j*inc_matr;
  % cut_i(cut_i==2) = 0;  % some functions might require this infos (internal lines of ith cluster)
  cutind(i,:) = cut_i;
end
end
