function ok = cmpIncAdj(I, A)
% Syntax: ok = cmpIncAdj(I, A)
%
% Purpose: 
%   cmpIncAdj(double I, double A) checks if the incidence matrix I and 
%   the adjacency matrix A are consistent with each other. The 
%   weights in matrix A are neglected for this task.
%
% Input: 
%   I - incidence matrix of a graph (mxn)
%   A - (weighted) adjacency matrix of a graph (mxm)
%        
% Output: 
%   ok - is false if matrices don't agree, true otherwise
% 
% Author: Ilya Tyuryukanov
% Date: 20 May 2015
% Last revision: 02 August 2015

siz = size(I); 
m = siz(1);
n = siz(2);
epsilon = eps('double');

a = I*I';
b = A;
b(b~=0) = 1;  % convert a weighted graph to simple graph
d = sum(b,2);
b = b + diag(d);
c = abs(a - b);

if sum(sum(c))>10*max(d)*epsilon
  ok = false;
else
  ok = true;
end
end