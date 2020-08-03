function graph_py = pfgraph2igraph( g, clust )
%
% Purpose: Converts a PFgraph object g to an igraph object graph_py for
%   the further usage with the Python igraph package (and with Python in
%   general).
%   The second input can model arbitratry node groupings (e.g. cluster 
%   assignments for every bus). By default (without 2nd input) generator 
%   coherency grouping in g.coh is passed as the node groping information
%   to graph_py. Thus, the format of the 2nd input is the same as the  
%   format of g.coh. Each group of nodes can get its own color for the 
%   visualization.
%
% Author: Ilya Tyuryukanov
% Date of first version: 1 October 2015
% Revision of 2 October 2015 (bug fixes?)
% Revision of 5 February 2016: conversion to a PFgraph method

adj = g.adj;
busmap = g.bus;
ml = g.ml;
if nargin == 1
  col = g.coh;
else
  col = clust;
end

% Note that all isolated vertices were removed in advance (no zero rows
% in adj), so the dimensions of 'adj' are preserved in 'i' and 'j' FOR
% INITIAL CONNECTED GRAPHS. The case of splitted networks should be tested
% however.
[i, j, w] = find(tril(adj));

% Convert everything to Python
% ('create -> append' is robuster, e.g. if 'i, j, w, ...' are scalars)
% -1 is because nodes will be indexed from 0 in Python
i_py = py.list(num2cell(uint32(i' - 1)));
j_py = py.list(num2cell(uint32(j' - 1)));
ws_py = py.list(num2cell(w'));
es_py = py.list(py.zip(i_py, j_py, ws_py));
bus_py = py.list(num2cell(uint32(busmap(1,:))));
if ~isempty(ml) && size(ml, 2) == 2
  ML1 = py.list(num2cell(uint32(ml(:, 1)' - 1)));
  ML2 = py.list(num2cell(uint32(ml(:, 2)' - 1)));
  % 2->1 order since we've took tril(adj) and ml(:,2)>ml(:,1) due to createLfGraph.m
  ml_py = py.list(py.zip(ML2, ML1));
else
  ml_py = py.None;
end

if ~isempty(col)
  col_py = py.list();
  if size(col,2)==3 && all(col(:)<=1) && all(col(:)>=0) % if rgb colors are provided
    for i = 1:1:size(col,1)
      grpi_py = py.tuple(num2cell(col(i,:)));
      col_py.append(grpi_py);      
    end
  elseif size(col,1)==2 && all(col(:)>=1) && Utils.isint(col)  % if node groupings are provided
    for i = unique(col(2,:))
      idx_grpi = col(2,:) == i;
      grpi = col(1, idx_grpi);
      grpi_py = py.tuple(num2cell(uint32(grpi - 1)));
      col_py.append(grpi_py);
    end
  else  % if something else (invalid) 
    error('')    
  end
else
  col_py = py.None;
end

% Change current directory to the directory of the method
old_pwd = pwd;
reset_pwd = onCleanup(@() cd(old_pwd));
file_dir = fileparts(mfilename('fullpath'));
cd(file_dir);

% Amend Python path (relative to the directory of the method)
py.importlib.import_module('sys');
pythonpath = py.sys.path;
if count(pythonpath, '..') == 0, 
  insert(pythonpath,int32(0), '..'); % one lvl above @PFgraph (where pfgraph2igraph.m resides)
end

% Reload the Python module each time (since it's in development...)
% (+ 'clear classes' from time to time, if reload fails to return good results)
mod = py.importlib.import_module('igraphmod');
py.importlib.reload(mod);
graph_py = py.igraphmod.mat2igraph(es_py, bus_py, ml_py, col_py);


