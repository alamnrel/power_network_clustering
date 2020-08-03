function [EigenVectors, R] = yushi_discrete(EigenVectors, R0, nbIterationsDiscretisationMax)
% 
% EigenvectorsDiscrete=yushi_discrete(EigenVectors)
% 
% Input: EigenVectors = continuous Ncut vector, size = ndata x nbEigenvectors 
% Output EigenvectorsDiscrete = discrete Ncut vector, size = ndata x nbEigenvectors
%
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

if nargin<3
  nbIterationsDiscretisationMax = 150;  %voir
end

[n,k]=size(EigenVectors);
i_rand = randi(n, min(20,n),1);
test = EigenVectors(i_rand,:);
test_length = sqrt(sum(test.^2, 2));

% if 20-random sample is unnormalised, normalise
if ~all(abs(test_length-1)<1e-10)
  vm = sqrt(sum(EigenVectors.*EigenVectors,2));
  EigenVectors = EigenVectors./repmat(vm,1,k);
end

if nargin<2
  R = eye(k);
else
  R = R0;
end
lastObjectiveValue = 0;
exitLoop = 0;
nbIterationsDiscretisation = 0;
epsilon = eps('double');
% nbIterationsDiscretisationMax = 20; %voir
while exitLoop== 0 
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;   
    EigenvectorsDiscrete = discretize(EigenVectors*R);
    [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors,0);    
    NcutValue=2*(n-trace(S));
    
    if abs(NcutValue-lastObjectiveValue) < epsilon || nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end
EigenVectors = EigenVectors*R;
end