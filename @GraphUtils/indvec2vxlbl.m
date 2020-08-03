function vxlbl = indvec2vxlbl(busind)
% 
% Converts the matrix of index vectors for each partition into a vector of
% labels for each vertex, where the label denotes the number of partition
% for that vertex.
% The index vectors assumed to be row vectors, i.e. the number of rows is 
% equal to the number of partitions, with each row beeing a binary index
% vector.
% As a side note, given vxlbl, the busind matrix of binary index vectors
% can be obtained with sparse(vxlbl, 1:n_vx, true). 
% In case the partitions in busind do not cover all graph vertices, the
% remaining vertices (not covered by any of the partitions) receive the
% label NaN.
%
% Author: Ilya Tyuryukanov
% Date: 15 September 2016   
% Last revision: 15 September 2016


% BY REMOVING THE ASSERTIONS WE ALLOW FOR SOME NaN VERTEX LABELS WHICH CAN 
% BE FILTERED OUT LATER.

% assert( all(sum(busind,1) == 1),...
%   [mfilename,':WrongPublicInput'],...
%   ['[%s] The input matrix is in the wrong format: it should contain binary ',...
%   'index vectors of partitions in its rows.'], mfilename);

n_part = size(busind, 1); 
n_vx = size(busind, 2);
vxlbl = NaN(1, n_vx);

for i = 1:n_part
  vxlbl(logical(busind(i, :))) = i;
end

% assert(~any(isnan(vxlbl)),...
%   [mfilename,':InconsistentResults'],...
%   ['[%s] The output contains NaNs. Not all vertices were assigned ',... 
%   'their partitions.']);

end