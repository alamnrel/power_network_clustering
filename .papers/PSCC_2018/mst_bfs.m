function [d dt pred childs] = mst_bfs(A,root)
% 
% MST_BFS Compute breadth first search distances, times, tree and children 
% for a rooted minimum spanning tree.
%
% [d dt pred] = bfs(A,u) returns the distance (d) and the discover time
% (dt) for each vertex in the spanning tree A in a breadth first search 
% starting from the tree ROOT (2nd argument).
%   d = dt(i) = -1 if vertex i is not reachable from ROOT
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from ROOT and i != ROOT.
%

% A more specialized version of bfs() by David F. Gleich
% Copyright, Stanford University, 2008-2009

% History:
% 2008-04-13: Initial coding
% October 2017: Adapted by me to return children of each node of a rooted 
%               spanning tree

[rp, ci]=sparse_to_csr(A);
num_nbr = diff(rp);
num_nbr = num_nbr-1;  % in a tree ONE neighbour is parent, others are children
num_nbr(root) = num_nbr(root)+1;  % but ROOT has no parents

n=length(rp)-1; 
d=-1*ones(n,1); dt=-1*ones(n,1); pred=zeros(1,n);
%childs = arrayfun(@(x) zeros(x,1), num_nbr, 'unif', false);
childs = cell(n,1);
for k = 1:1:n
  childs{k} = zeros(num_nbr(k),1);
end
curr = ones(n,1);
sq=zeros(n,1); sqt=0; sqh=0; % search queue and search queue tail/head

% start bfs at u
sqt=sqt+1; sq(sqt)=root; 
t=0;  
d(root)=0; dt(root)=t; t=t+1; pred(root)=root;
while sqt-sqh>0
    sqh=sqh+1; v=sq(sqh); % pop v off the head of the queue
    for ri=rp(v):rp(v+1)-1
        w=ci(ri);
        if d(w)<0            
            childs{v}(curr(v)) = w; curr(v) = curr(v)+1;
            sqt=sqt+1; sq(sqt)=w; 
            d(w)=d(v)+1; dt(w)=t; t=t+1; pred(w)=v; 
        end
    end
end
