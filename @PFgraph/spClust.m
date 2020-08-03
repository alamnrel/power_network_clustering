function [Y, v] = spClust(g, nm, nc)
% 
% return first nc eigenvalues of L (sorted) and the corresponding first nc 
% eigenvectors of L if is false.
% nm defines the type of the Laplacian property used for the eigenvectors
% calculation
%
% Author: Ilya Tyuryukanov
% Date of first version: 18 August 2015
%
% Revisions: 
% 1) 5 November 2015 (??)
% 2) 28 April 2019 (removed the nm input to GraphUtils.spClust)

switch nm
  case 'off'
    L = g.Loff;
  case 'rwk'
    L = g.Lrwk;
  case 'sym'
    if ~isempty(g.Lsym)
      L = g.Lsym;
    else
      L = g.Lrwk;
    end
  otherwise
    error([mfilename,':WrongPublicInput'],...
      ['[%s] Invalid Laplacian normalisation method in PFgraph.spClust(g, nm, nc, out). ',...
      'Use nm = ''off''/''rwk''/''sym''.'], mfilename);
end

assert(~isempty(L), [ mfilename, ':InvalidProperty'],...
  ['[%s] The chosen Laplacian entry of the input graph object wasn''t ',...
  'initialized properly'], mfilename);

[Y, v] = GraphUtils.spClust(L, 'sm', nc);
end