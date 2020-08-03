function g_redu0 = reduce_pml(g)
%
% Reduction of pairwise must-link constraints (e.g., transformer branches,
% critical transmission lines etc). 
%
% This step must be performed before searching for generator coherency
% cutset (e.g., with Packing Steiner Trees heuristics or MILP). Otherwise 
% cutsets found by generator coherency cutset algorithms may be invalid
% (e.g., MILP separates a transformer, or a Packing Steiner Trees heuristic
% constructs generator trees that have a connection through a transformer).
% 
% Author: Ilya Tyuryukanov
% Date of first version: 14 May 2019

% Aggregate pairwise must links 
if ~isempty(g.ml)
  ixs = GraphUtils.aggregate_mls(g.ml);  
  g_redu0 = merge_nodes(g, ixs);
  
  % filament is not used, as it merges "isolated/filament" generators 
  % unpredictably to the previous shortest paths by aggregate_mls (try 
  % case89pegase with 5 groups). The function behavior will change too.
  % So this is given as a theoretic possibility only.
  
  %{
  % filament = g.find_filaments();  
  %}
else
  g_redu0 = g;
end

end



