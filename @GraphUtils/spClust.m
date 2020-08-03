function [Y, v] = spClust(L, lo_up, nc)
%
% return first or last smallest nc eigenvalues and eigenvectors of L (sorted).
% Returning first or last depends on the LO_UP input (set to 'sm' or 'la' 
% respectively, or a numeric value about which to search). 
% 
% nm defines whether the eigenvectors are normalized to length 1 or left
% as they are as well as the type of the Laplacian property used for their
% calculation
%
% Author: Ilya Tyuryukanov
% Date of first version: 18 August 2015
%
% Revisions: 
% 1) 5 November 2015 (??)
% 2) 28 April 2019 (removed row normalization for nm='sym' and nm input itself, as it is irreversible and not always good)

opts.isreal = 1;
opts.disp = 0;
opts.maxit = 10000;
opts.tol = 2.2204e-16;


switch lo_up
  case 'lm'
    srt_ord = 'descend';
  case 'la'
    srt_ord = 'descend';
  case 'sm'
    srt_ord = 'ascend';
  case 'sa'    
    srt_ord = 'ascend';
  otherwise
    error(['The lower/upper eigenvalue switch is invalid. ',...
      'Choose among lm, la, sm, sa. See help of EIGS for more info.'])
end


% Construct R^nc (spectral nc-embedding)
v = NaN(5,1);  % dummy v!
while any(isnan(v)) && opts.maxit<1.1e6
  [X, v] = eigs(L, nc, lo_up, opts);
  X = real(X);  % just in case...
  if (size(v,1) > 1) && (size(v,2) > 1)
    v = diag(v);  % if condition is for robustness
  end
  v = real(v);  % just in case...
  if any(isnan(v))
    opts.maxit = opts.maxit*100;
  else    
    [v, idx] = sort(v, srt_ord);
    X = X(:,idx);
  end
end

if any(isnan(v))
  error('');
end

Y = X;
end