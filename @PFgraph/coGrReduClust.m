function [T] = coGrReduClust(g, trees, varargin)
% 
% Syntax: [T] = coGrReduClust(g, trees, varargin)
% 
% Purpose: performs Spectral Clustering on active power flow graph of power 
%   network model. The graph is reduced before Spectral Clustering is  
%   applied. In particular, must-link constraints are implemented by vertex
%   contractions. For that purpose, coherent generators should be first 
%   connected into subnetworks (the tree input). Cannot-link constraints are  
%   implemented by Voronoi-type  clustering around contracted coherent nodes 
%   (i.e. each coherent group is represented by a single node).
%          
% Input: 
% g: a (sensible) PFgraph object (with adjacency, incidence matrices and 
%   constraints)
% varargin: key-value pairs for this function (a struct or an isomorphic 
%   cell array). The keys and their possible values are listed below.
%                    
%   sc_alg: Laplacian normalization method for Spectral Clustering
%     'off': no normalization 
%     'rwk': random-walk normalized Laplacian
%     'sym': symmetric normalized Laplacian  
%   refine: If true and isunix is true, attempt label-propagation-based cut
%     refinement post-processing to improve the output partition.
%
% Output: 
% T: clustering assignment vector of the vertices of the final graph.
% 
% MAIN REFERENCES
% Many ideas about graph simplification are (initially) from: Guangyue Xu,
% Vijay Vittal (2010), "Slow Coherency Based Cutset Determination Algorithm
% for Large Power Systems", In IEEE Trans. Power Syst. 25 (2), pp. 877-884 
% 
% Author: Ilya Tyuryukanov
% Date of first version: 09 Sep 2015 
% Revision of 04 February 2016: almost totally rewritten as a PFgraph 
%    method, based on other PFgraph methods.     
% Revision somewhere down the road: adapted the function to generic "trees"
%    inputs that are no longer computed inside the function.
% Revision of 14 May 2019: Removed the reduction of pairwise must-links, as
%    it should happen BEFORE computing the trees that are reduced inside of
%    this function.

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('sc_alg', 'sym', @ischar);
p.addParameter('refine', false, @(x)isscalar(x)&&islogical(x));
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
  'spelling. The allowed keys are: sc_alg.'], mfilename);
nm = varinp.sc_alg;
rf = varinp.refine;

lnorm = {'off', 'rwk', 'sym'};  % Laplacian normalizations
assert(ismember({nm},lnorm), [ mfilename, ':WrongPublicInput'],...
  ['[%s] Only ''off'', ''rwk'' or ''sym'' Laplacian normalization codes are allowed ',...
  'as the ''sc_alg'' key-value input'], mfilename);
assert(~isempty(g.inc), [ mfilename, ':WrongPublicInput'],...
  ['[%s] The input graph must contain its incidence matrix as this is used to ',...
  'calculate the cutset from bus grouping information'], mfilename);

% The aggregation of pairwise must links (may be together with filaments),
% if any are present, must occur before computing the trees. The provided 
% input graph and the trees should already include the reduction of
% pairwise must-link branches to get the correct results.
assert(isempty(g.ml));

% Next, aggregate trees
trees = cellfun(@unique, trees, 'un', false);
g_redu = merge_nodes(g, trees);
n_coh = numel(trees);  
n_red = size(g_redu.adj, 1);
anchrs = n_red-n_coh+1:1:n_red;
[anch_i, anch_j] = find(g_redu.adj(anchrs,anchrs));
anch_i = anch_i + n_red - n_coh;
anch_j = anch_j + n_red - n_coh;
anch_p = sub2ind([n_red,n_red], anch_i, anch_j);
g_redu.adj(anch_p) = 0.1*min(abs(nonzeros(g_redu.adj)));   % reduce the weight between already merged "anchor nodes"

% Bus clustering (on the reduced graph)
g_redu = g_redu.createLaplacian(nm);
if strcmp(nm, 'sym')
  Asym = -g_redu.Lsym;
  m = size(Asym, 1);
  Asym(1:m+1:m*m) = 0;
  Yredu = GraphUtils.spClust(Asym, 'la', n_coh);
  Yredu = Utils.normalize_rows(Yredu);
else
  Yredu = g_redu.spClust(nm, n_coh);
end
emb_graph = GraphUtils.buildDistanceGraph(Yredu, g_redu.adj, nm);
dm = GraphUtils.nShortestPaths(emb_graph, anchrs, 'orig');  % num_lod x num_coh_ml distance matrix
[~, Tredu] = min(dm, [], 2);

cut0 = final_cutset(g_redu, Tredu);
[eXp0, pcut0] = g_redu.ici_info( cut0 );
if isunix && rf
  goals.weXp = Inf;
  goals.nbus = 1;
  goals.cut = pcut0;
  goals.Ncut = Inf;
  [T0, cntg0] = kahip_improve_v1(g_redu, Tredu, goals);
  if ~cntg0
    [T0, minor_cc] = contig(g_redu, T0);
  end
  cut1 = final_cutset(g_redu, T0);
  [eXp1, pcut1] = g_redu.ici_info( cut1 );
  if pcut1<pcut0
    Tredu = T0;
  end
end

% Go back to the original graph representation and form islands
T = Tredu(g_redu.merge_map);

% For plots (upon request)
%{
m_redu = size(g_redu.adj, 1);
row_gen_redu = m_redu-n_coh+1:m_redu;
row_lod_redu = 1:m_redu-n_coh;
%}


