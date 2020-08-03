function [T, UBfactor, contig] = ucHMETIS(g,method,fix,imb)

% Purpose:
% Interface to call hmetis from the command line.
% Version used originally: 1.5.3, but for large networks (case300,
% case2383wp) resulted in the error "Out of netind memory!". This version
% is written for 32bit architectures and the architecture used was 64bit.
% Same issue in both windows 10 and Ubuntu 16.04.1 LTS.
% Version used: 2.0pre1
% For more information (no manual available): ./hmetis2.0pre1 -help
% See also the manual from 1.5.3 (many things are similar)
%
% Input:
% g: a (sensible) PFgraph object.
% method: parameter that specifies which method will be used to partition
%         the graph. method=0: multilevel recursive bisection.
%                    method=1: multilevel k-way partitioning.
% fix: parameter that sepcifies if the file that pre-assigns the generators
%      to certain partitions according to their coherency will be used.
%   fix=0: the file will not be used.
%   fix=1: the file will be used.
%
% Note that all weights on hMetis must be integers. Edge weights must be
% positive and vertex weights must be greater or equal to zero (better
% avoid using zero weights)
%
% Output:
% T: clustering assignment vector of the vertices of the final graph.
% UBfactor: parameter that specifies the imbalance between the partitions.
% contig: parameter that shows if the initial number of areas created are
% different than the number of connected components.
%   contig=0: bad results (partition is not valid).
%   contig=1: good results (partition is valid).

% Check the input parameters
assert(method==0 || method==1, '[%s] The specified method is not correct. Type 1 to use multilevel k-way partitioning or 0 to use multilevel recursive bisection.', mfilename);
assert(fix==0 || fix==1, '[%s] The specified fix number is not correct. Type 1 to use fixfile or 0 not to use it.', mfilename);

np = max(g.coh(2,:));  % number of clusters (in this context)
assert(np >= 2, [mfilename,':WrongPrivateInput'],...
  ['[%s] The input number of clusters is less than two. Check the ''coh'' ',...
  'attribute of the input graph object (its 2nd row).'], mfilename);
m = size(g.adj, 2);

vw = 1:1:m;
vw = vw(:);
GraphUtils.pfgraph2hmetis(g.adj, g.inc, vw, g.coh);

if nargin < 4
  all_imbalance=[9 12 15 18 21 24 27 30 33 37 40 43 47 50 55 60 65 70 75 100 150 300 500 750 1000];  % -> the same as in ucKaHIP
else
  all_imbalance = imb;
end
n_cases=length(all_imbalance);
all_bus = cell(n_cases, 1);
weXp = ones(n_cases, 1);
nmin = ones(n_cases, 1);
if fix==1
  outl_gen = size(g.coh,2)*ones(n_cases, 1);
  dif = m*ones(n_cases, 1);
end
cntg = false(n_cases, 1);
cc_param = cell(n_cases, 1);
n_cntg = 0;
n_fail = 0;

i = 1;
while i <= n_cases
  if fix==0
    if method==0
      system(sprintf('hmetis2.0pre1 test_hmetis.hgr %d -rtype=slow -ufactor=%d -nruns=20 -nvcycles=20  -kwayrefine -cmaxnet=150 -rmaxnet=150 >/dev/null',np,all_imbalance(i))); %-kwayrefine
    else
      system(sprintf('hmetis2.0pre1 test_hmetis.hgr %d -ptype=kway -rtype=kpslow -ufactor=%d -nruns=20 -nvcycles=20 -cmaxnet=150 -rmaxnet=150 >/dev/null',np,all_imbalance(i)));
    end
    map = transpose(load(sprintf('test.graph.part.%d',np)));
    [cuts, bus, iscontig, eXp] = g.metis_info( map );
    map = map + 1;
    cntg(i) = logical(iscontig);
    [~, idx] = sort(eXp(3,:), 'descend');
    eXp = eXp(:,idx);
    if (sum(eXp(3,np+1:end))~=0) && (sum(eXp(3,np+1:end))<1*eXp(3,np))
      [map] = g.contig( map, 'SplitLargeMinor', true );
      [cuts, bus] = final_cutset(g, map);
      [~, ~, ~, eXp] = g.cc_info( cuts );
      cntg(i) = true;
      n_cntg = n_cntg + 1;
    end
    weXp(i) = max(eXp(1,:));
    nmin(i) = min(eXp(3,:));
    cc_param{i} = eXp;
    all_bus{i} = map;
  else
    if method==0
      stat = system(sprintf('timeout 25 hmetis2.0pre1 test_hmetis.hgr %d -ptype=rb   -rtype=slow   -reconst -ufactor=%d -nruns=25 -nvcycles=25 -cmaxnet=3 -rmaxnet=3 -fixed=fixfile.txt >/dev/null',np,all_imbalance(i)));  % -kwayrefine  
    else
      stat = system(sprintf('timeout 25 hmetis2.0pre1 test_hmetis.hgr %d -ptype=kway -rtype=kpslow -ufactor=%d -nruns=25 -nvcycles=25 -cmaxnet=3 -rmaxnet=3 -fixed=fixfile.txt >/dev/null',np,all_imbalance(i)));
    end
    if stat == 0    
      map = transpose(load(sprintf('test_hmetis.hgr.part.%d', np)));
    else
      map = randi(m, 1, m) - 1;  % to result in a huge outl_gen(i)      
    end
    [cuts, bus, iscontig, eXp, pcut, n_outl_gen, difr] = g.metis_info(map);    
    map = map + 1;
    cntg(i) = logical(iscontig);
    % Applying contig.m in label-constrained setting is not reasonable, as
    % contig.m is not optimized to reassign labelled components (15.11.2017)
    weXp(i) = max(eXp(1,:));
    nmin(i) = min(eXp(3,:)); 
    cc_param{i} = eXp;
    all_bus{i} = map;
    outl_gen(i) = n_outl_gen;
    dif(i)=difr;
    i = i+1;    
  end
end

% Store the output of the best partition
if fix==0
  [T, i_best, contig] = GraphUtils.best_part(cntg,all_bus,weXp);
else
  [T, i_best, contig] = GraphUtils.best_part(cntg,all_bus,weXp,outl_gen,dif);
end
UBfactor = 1+all_imbalance(i_best)/100;
end

