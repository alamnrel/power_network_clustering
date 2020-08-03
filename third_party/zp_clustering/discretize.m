function Y = discretize(EigenVector)
% Y = discretisationEigenVectorData(EigenVector)
%
% discretizes previously rotated eigenvectors in discretisation
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

[n,k]=size(EigenVector);


[Maximum,J]=max(EigenVector,[],2);

Y = zeros(n,k);
idx = sub2ind([n,k], 1:n, J');
Y(idx) = 1;
end