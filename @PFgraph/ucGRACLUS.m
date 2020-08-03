function [T, n_stat] = ucGRACLUS(g, varargin)
% 
% Purpose:
% Uses the mex interface of graclus graph partitioning package to find
% clusters in g.adj (i.e. graph adjacency matrix without constraints).
% 
% Input:
% g: a (sensible) PFgraph object.
% 
% Note that all weights on graclus must be integers. Edge weights must be
% positive and vertex weights must be greater or equal to zero (better
% avoid using zero weights).
%
% Output:
% T : clustering assignment vector of the vertices of the final graph.    
% n_stat: some solution statistics...
%

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('cut_typ', 0, @(x) x == 0 || x == 1);
p.addParameter('n_ls', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('spec', 1, @(x) x == 0 || x == 1);
p.addParameter('n_rep', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
  'spelling. The allowed keys are: cut_typ, n_ls, spec, n_rep.'], mfilename);
cut_typ = varinp.cut_typ;
n_ls = varinp.n_ls;
spec = varinp.spec;
n_rep = varinp.n_rep;

np = max(g.coh(2,:));  % number of clusters (in this context)
assert(np >= 1, [mfilename,':WrongPrivateInput'],...
  ['[%s] The input number of clusters is less than two. Check the ''coh'' ',...
  'attribute of the input graph object (its 2nd row).'], mfilename);
m = size(g.adj, 2);
adj = g.adj;

all_bus = cell(n_rep, 1);
weXp = zeros(n_rep, 1);
nmin = false(n_rep, 1);
cntg = false(n_rep, 1);
cc_param = cell(n_rep, 1);
n_cntg = 0;
n_fail = 0;

for i = 1:1:n_rep   
  [map, obj] = graclus(adj, np, cut_typ, n_ls, spec);
   map = map - 1;  % to make consistent with METIS
  [cuts, bus, iscontig, eXp] = g.metis_info( map );  
  cntg(i) = logical(iscontig);
  [~, idx] = sort(eXp(3,:), 'descend');
  eXp = eXp(:,idx);  
  if (sum(eXp(3,np+1:end))~=0) && (sum(eXp(3,np+1:end))<1*eXp(3,np))
    [cuts, bus] = g.contig( cuts, 'SplitLargeMInor', true );
    [~, ~, ~, eXp] = g.cc_info( cuts );
    cntg(i) = true;
    n_cntg = n_cntg + 1;
  end
  weXp(i) = max(eXp(1,:));
  nmin(i) = min(eXp(3,:)) >= ceil(m/np*0.2);
  cc_param{i} = eXp;
  all_bus{i} = map + 1;
end

% Store the output of the best partition
[T, i_best] = GraphUtils.best_part(cntg&nmin, all_bus, weXp);
if ~any(cntg & nmin)
  n_fail = n_fail + 1;
end
n_stat = [n_cntg, n_fail];

end

