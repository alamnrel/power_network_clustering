function [T] = ucSpEmbGrPart(g, varargin) 
% Syntax: 
% [T] = ucSpEmbGrPart(g, varargin) 
% 
% Purpose: 
% Performs spectral clustering on the active power flow graph of a power 
% network model. The graph should be supplied as a PFgraph object containing 
% among other things (weighted) adjacency and incidence matrices. Several
% Laplacian and postprocessing options can be chosen by specifying the key
% value arguments.
% 
% Input: 
% g:    a (sensible) PFgraph object  
% varargin: key-value pairs for this function (a struct or an isomorphic 
%   cell array). The keys and their possible values are listed below. 
%                   
%   alg(1:3): Laplacian normalization method
%     'off': no normalization 
%     'rwk': random-walk normalized Laplacian 
%     'sym': symmetric normalized Laplacian   
%   alg(4): Interpretation of the obtained points in spectral embedding
%     '+' means clustering of the embedded points in the Euclidian space
%     '-' means creation of an "embedded graph" and clustering of that graph
%   postpr: Way of clustering of points in spectral embedding
%     'pic' for PIC partitioning (GAC toolbox, Wei Zhang & Co) of the
%        embedded graph (in spectral embedding)
%     'kmeans': kmeans clustering as in MATLAB Statistics and Machine 
%        Learning Toolbox (only valid for alg(4) equal to '+')
%     'kmedoids': kmedoids clustering as in MATLAB Statistics and Machine 
%        Learning Toolbox (only valid for alg(4) equal to '+') 
%     OR hierarchical clustering linkage (a string as in the linkage()
%        function of the MATLAB Statistics and Machine Learning Toolbox)
% 
% Output: 
% T: clustering assignment vector of the vertices of the final graph.
% 
% MAIN REFERENCES
% The idea of using "embedded graphs" stems from the paper
% R. J. Sanchez-Garcia, M. Fennelly, S. Norris, N. Wright, G. Niblo,  
% J. Brodzki and J. W. Bialek, "Hierarchical Spectral Clustering of 
% Power Grids", In IEEE Trans. Power Syst. 29 (5), pp. 2229 2237
% 
% Author: Ilya Tyuryukanov
% Date of first version: 8 June 2015 
% Revision of 14 January 2016: rewritten as a PFgraph method + testing 

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('alg', 'sym-', @ischar);
p.addParameter('postpr', 'average', @ischar);
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
  'spelling. The allowed keys are: alg, postpr.'], mfilename);
alg = varinp.alg;
[nm, pts_repr] = parsealg(alg);
postpr = varinp.postpr;

hlink = {'single', 'weighted', 'average', 'complete', 'ward', 'centroid',...
    'median'};  % "hierarchical" linkages
slink = {'pic'};  % "similarity" linkages
plink = {'kmeans', 'kmedoids'};  % point clustering in Euclidian space
lnorm = {'off', 'rwk', 'sym'};  % Laplacian normalizations

assert(ismember({postpr},hlink) || ismember({postpr},slink) || ismember({postpr},plink),...
  [ mfilename, ':WrongPublicInput'],...
  '[%s] Only predefined ''postpr'' strings are allowed', mfilename);
assert(ismember({nm},lnorm), [ mfilename, ':WrongPublicInput'],...
  ['[%s] Only "off", "rwk" or "sym" Laplacian normalization codes are allowed ',...
  'as characters  1--3 of the "alg" input'], mfilename);
assert(strcmp(pts_repr,'eucld')&&(ismember({postpr},hlink)||ismember({postpr},plink))||...
  strcmp(pts_repr,'embgr')&&(ismember({postpr},hlink)||ismember({postpr},slink)),...
  [ mfilename, ':IncompatibleInput'],...
  ['[%s] Only point-cloud (e.g. kmeans) or hierarchical clustering linkages ',...
  'are allowed with alg="xxx+". Only hierarchical clustering or ',...
  'clustering-with-similarity-graphs linkages are allowed with alg="xxx-"'],...
  mfilename);

% Initialization
persistent run_count
if isempty(run_count), run_count = 1; end
adj = g.adj;
nc = max(g.coh(2,:));  % number of clusters (in this context)
assert(nc >= 1, [mfilename,':WrongPrivateInput'],...
  ['[%s] The input number of clusters is less than two. Check the ''coh'' ',...
  'attribute of the input graph object (its 2nd row).'], mfilename);

% Bus clustering
g = g.createLaplacian(nm);
if strcmp(nm, 'sym')
  Asym = -g.Lsym;
  m = size(Asym, 1);
  Asym(1:m+1:m*m) = 0;
  Y = GraphUtils.spClust(Asym, 'la', nc);
  Y = Utils.normalize_rows(Y);
else
  Y = g.spClust(nm, nc);
end
[T] = GraphUtils.sp_postpr(Y, postpr, alg, nc, adj);
end

function [nm, pts_repr] = parsealg(alg)
% Parse string with algorithm type and return its meaningful parts
% (normalization mode and points representation).
% 
% Author: Ilya Tyuryukanov
% Last revision: 14 January 2016

assert(numel(alg)==4, [ mfilename, ':WrongPublicInput'],...
  '[%s] The alg string can be only four characters long', mfilename);

nm = alg(1:3);
if strcmp(alg(4), '+')
  pts_repr = 'eucld';  % just as points in Euclidian space
elseif strcmp(alg(4), '-')
  pts_repr = 'embgr';  % as 'embedded graph'
else
  error('[parsealg] The 4th character of the alg string can be either "+" or "-"');
end
end


% Legacy storing graphic output..
%{
if ~strcmp(pics,'><?*')
  try
    plt = matfile(fullfile(pics, mfilename));
    plt.Properties.Writable = true;
    plt.for_plot(run_count, 1) = struct('Y', Y, 'busmap', g.bus,...
      'busind', T);
    run_count = run_count + 1;
  catch me
    warning(['[%s] Something went wrong when trying to save the ',...
      'graphical output. See the information below.'], mfilename);
    disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
  end
end
%}
