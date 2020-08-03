function pfgraph2hmetis(adj, inc, vw, coh)

% Purpose:
% Creates the input files needed by hMETIS to perfrom the
% partitioning
%
% Input:
% adj: a valid  
%
%
% Output:
% The .hgr file for input to hMETIS
% The .txt fix file for input to hMETIS

if nargin<4
  fix = 0;
else
  fix = 1;
end
m=size(inc,1);  % number of vertices
n=size(inc,2);  % number of edges
if nargin < 3
  vw = ones(m, 1);
end
vw = vw(:);

% Create the hgr file for input in hmetis
file=fopen('test_hmetis.hgr','wt');
cleanup0 = onCleanup(@() fclose(file));
[~, ~, w] = find(adj);
if min(w) < 100  %  abs(min(w)-1)>100*epsilon
  mf = min(1e7,100/min(w));
  adj = ceil(adj*mf);  % because the weights need to be rounded, but avoiding zeros!
  w = ceil(w*mf);
end

if any(w~=1)
  % Create hypergraph with weighted edges
  fprintf(file,'%d %d 11\n',n,m);
  for i=1:1:n
    edge=find(inc(:,i));
    egdeweight=adj(edge(1,1),edge(2,1));
    fprintf(file,'%d %d ',ceil(abs( full(egdeweight) )), edge);  % weights on edges must be positive integers greater or equal to 1
    fprintf(file,'\n');
  end
  fprintf(file,'%d\n', vw);
else
  % Create unweighted hypergraph
  fprintf(file,'%d %d 10\n',n,m);
  for i=1:1:n
    edge=find(inc(:,i));
    fprintf(file,'%d ',edge);
    fprintf(file,'\n');
  end
  fprintf(file,'%d\n', vw);
end

% Create FixFile
if fix==1
  ff=-1*ones(m,1);
  for i=1:1:size(coh,2)
    ff(coh(1,i),1)=coh(2,i)-2;  % partition numbers start from 0, but 'free vertices' have index -1
  end
  fileff=fopen('fixfile.txt','wt');
  cleanup1 = onCleanup(@() fclose(fileff));
  fprintf(fileff,'%d\n',ff(:,1));
end
end