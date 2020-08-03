function g_out = extract_subgraph(g, busind)   
% 
% Syntax: g_out = extract_subgraph(g, busind)   
% 
% Purpose: Given the original power network graph (as a PFgraph object g) 
%   and the indicator vector of a connected component of this graph, create
%   a new connected graph representing the indicated connected component.
% 
% Input: 
% g: a PFgraph object (incidence and adjacency matrices + constraints) 
% busind: an indicator vector into a connected component of g
% 
% Output: 
% g_out: a new connected graph representing the indicated connected 
%   component  
% 
% Author: Ilya Tyuryukanov 
% Date of first version: 26 November 2016 

adj = g.adj;
m = size(adj,1);

% % The commented code below is just to produce a (superfluous!) warning
% busidx = busind + 1; 
% busidx = sparse(busidx, 1:m, true);
% cutind = cutset(g, busidx);
% cutind(cutind~=1) = 0;
% 
% adj_sep = cc_info(g, cutind);
% S = GraphUtils.alecconncomp(adj_sep);  % optional!
% if S>2
%   warning([ mfilename, ':BadResults'],...
%     '[%s] The input index vector represents more than one connected component',...
%     mfilename);
% end

busind = logical(busind);
adj_out = adj(busind,busind);

inc_out = g.inc;
inc_out = inc_out(busind, :);
idx = sum(inc_out,1)<2;  % remove branches that don't fully belong to the remaining graph
inc_out(:,idx) = [];

vw_out = g.vw(busind);
gen_out = g.gen(busind);
bus_out = find(busind);  % "old" indices of "remaining" buses
map = zeros(1,m);    
map(bus_out) = 1:1:numel(bus_out);  % mapping: new_idx = map(old_idx) and old_idx = new_idx(map)

coh_out = g.coh;
if ~isempty(coh_out)
  coh_out(1,:) = map(coh_out(1,:));
  coh_out(:,coh_out(1,:) == 0) = [];
end

ml_out = g.ml;
if ~isempty(ml_out)
  col1 = map(ml_out(:,1));
  col2 = map(ml_out(:,2));
  ml_out = [col1(:), col2(:)];
  idx = ~(ml_out(:,1) & ml_out(:,2));  % ml branches that don't fully belong to the remaining graph
  ml_out(idx, :) = [];
end

g_out = PFgraph('adj', adj_out, 'inc', inc_out, 'merge_map', map, 'vw',...
  vw_out, 'coh', coh_out, 'ml', ml_out, 'gen', gen_out);

end