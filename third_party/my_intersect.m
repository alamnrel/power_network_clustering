function C = my_intersect(A,B)
% MYINTERSECT Intersection of two sets of positive integers (much faster than built-in intersect)
% C = my_intersect(A,B)
% "Nice, although it returns repeats of B which may or may not be desired"
%
% By Nick (from MatlabCentral), probably alike to the work of Leon Peshkin

% 

%A = A(:)'; B = B(:)';
A = A(:); B = B(:);  % columns are more efficient!
if isempty(A) || isempty(B)
  C = [];
  return
else
  bits = zeros(max(max(A),max(B)), 1);  %bits = sparse(max(ma,mb),1);
  bits(A) = 1;
  C = B(logical(bits(B)));  
end