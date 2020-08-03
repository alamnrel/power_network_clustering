function [LBL, avg_clu_dst] = vcs(S_GL, link, param)
% 
% This file implements "Clustering in Var Control Space" almost as in [1].
% (What is in [1] has a peculiar expression for the average zone distance).
% 
% The main input is the control-generator to loads sensitivity matrix S_GL
% obtained as S_GL = VCA.makeS(mpw, {'Q', 'Vm'}, SET_PQG, SET_LOD).
% 
% The way of obtainind the "average zone distance" is not well given in [1]
% (the results with formula in [1] are strange), thus here a more logical
% average cluster distance is used.
% 
% With this method (the original VCS [1]), the generators should be further 
% assigned to the nearest load clusters, and cluster connectivity issues 
% should be resolved.
% 
% MAIN REFERENCES:
% [1] H. Sun, Q. Guo, B. Zhang et al., "An adaptive zone-division-based 
%     automatic voltage control system with applications in China," In IEEE 
%     Trans. Power Syst., 28(2), pp. 1816-1828, May 2013
% [2] H. Ge, Q. Guo, H. Sun, B. Wang, B. Zhang, "Multivariate statistical 
%     analysis-based power-grid-partitioning method", IET 2017.
% 

if iscell(S_GL.S)
  python = true;
  S_GL.S = DiGSI.unpack_mat(S_GL.S);
  S_GL.buscol = DiGSI.unpack_mat(S_GL.buscol);
  S_GL.busrow = DiGSI.unpack_mat(S_GL.busrow);  
else
  python = false;
end

% Form the load-generator feature matrix X as in [1]
nld = size(S_GL.S, 1);
X = -log10(abs(S_GL.S)); 
bus = S_GL.busrow(:);

% Hierarchically cluster X 
if nargin<2
  link = 'ward';
end
D = pdist(X, 'euclidean');
[Z, avg_clu_dst] = linkageold(D, link(1:2));   % linkage tree
Z_tst = linkage(X, link);
assert( max(max(abs(Z-Z_tst)))<1e6 );
% dendrogram(Z)

% Test various cluster numbers
if nargin<3 || ~isfield(param, 'num_clu')
  clu_num = nld:-1:2;
  [d_max, i_max] = max(avg_clu_dst);
  nc = clu_num(i_max);
else
  nc = param.num_clu(end);
end
T_GL = cluster(Z, 'maxclust', nc);

% Return buses and labels
grp = unique(T_GL);
LBL = cell(numel(grp), 1);
bus = bus(:);
for k = 1:1:numel(grp)
  LBL{k} = bus(T_GL==grp(k));
end
if python
  LBL = cellfun(@(x) num2cell(x,2), LBL, 'un', 0);
  avg_clu_dst = num2cell(avg_clu_dst(:), 2);
end
end


function [Z, avg_dst] = linkageold(D, method)
%LINKAGEOLD Create hierarchical cluster tree using only MATLAB code.
% Copied from the MATLAB's Statistics Toolbox (linkage.m).

n = size(D,2);
m = round((1+sqrt(1+8*n))/2);
if isa(D,'single')
  Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
  Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if any(strcmp(method,{'ce' 'me' 'wa'}))
  D = D .* D;
end

avg_dst = NaN(n-1,1);
for s = 1:(n-1)
  if strcmp(method,'av')
    p = (m-1):-1:2;
    I = zeros(m*(m-1)/2,1);
    I(cumsum([1 p])) = 1;
    I = cumsum(I);
    J = ones(m*(m-1)/2,1);
    J(cumsum(p)+1) = 2-p;
    J(1)=2;
    J = cumsum(J);
    W = N(R(I)).*N(R(J));
    dst_avl = D./W;
    [v, k] = min(dst_avl);
    avg_dst(s) = mean(dst_avl);
  else
    [v, k] = min(D);
    if any(strcmp(method,{'ce' 'me' 'wa'}))
      avg_dst(s) = mean(sqrt(D));
    else
      avg_dst(s) = mean(D);
    end
  end
  
  i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
  j = k - (i-1)*(m-i/2)+i;
  
  Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
  
  % Update Y. In order to vectorize the computation, we need to compute
  % all the indices corresponding to cluster i and j in Y, denoted by I
  % and J.
  I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
  U = [I1 I2 I3];
  I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
  J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];
  
  switch method
    case 'si' % single linkage
      D(I) = min(D(I),D(J));
    case 'co' % complete linkage
      D(I) = max(D(I),D(J));
    case 'av' % average linkage
      D(I) = D(I) + D(J);
    case 'we' % weighted average linkage
      D(I) = (D(I) + D(J))/2;
    case 'ce' % centroid linkage
      K = N(R(i))+N(R(j));
      D(I) = (N(R(i)).*D(I)+N(R(j)).*D(J)-(N(R(i)).*N(R(j))*v)./K)./K;
    case 'me' % median linkage
      D(I) = (D(I) + D(J))/2 - v /4;
    case 'wa' % Ward's linkage
      D(I) = ((N(R(U))+N(R(i))).*D(I) + (N(R(U))+N(R(j))).*D(J) - ...
        N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
  end
  J = [J i*(m-(i+1)/2)-m+j];
  D(J) = []; % no need for the cluster information about j.
  
  % update m, N, R
  m = m-1;
  N(n+s) = N(R(i)) + N(R(j));
  R(i) = n+s;
  R(j:(n-1))=R((j+1):n);
end

if any(strcmp(method,{'ce' 'me' 'wa'}))
  Z(:,3) = sqrt(Z(:,3));
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);
end





