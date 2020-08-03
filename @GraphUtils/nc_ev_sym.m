function [Vr_out, Qt_out, t_run] = nc_ev_sym(V, rm, dm, batch_siz, k_min, debug)
% 
% Estimate a good number of groups by incrementally optimally rotating the
% eigenvectors (starting from k_min and up to k_max) and finding the number 
% of eigenvectors that corresponds to the best rotation.
% 

assert( sum(dm) <= 1 && size(dm,1) == 1  && size(dm,2) == 4 );
assert( size(rm,1) == 1  && size(rm,2) == 2 && sum(rm) <= 2 );

% Add path to Zelnik-Perona gradient descent code and kabsch
old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(which(mfilename))), 'third_party'));  % kabsch
addpath(fullfile(fileparts(fileparts(which(mfilename))), 'third_party', 'zp_clustering'));

k_max = size(V, 2);
if nargin < 5
  k_min = 2;
end
if nargin < 6
  debug = false;
end
nc = k_min:k_max;
m = size(V,1);

if nargin < 4 || isempty(batch_siz)
  batch_siz = round(linspace(round(m/10), floor(m/2), k_max-k_min+1));  % floor(m/2) is the maximal feasible batch size
else
  batch_siz = batch_siz*ones(1,k_max-k_min+1);
end
Vr_out = cell(1, k_max-k_min+1);
Qt_out = NaN(1, k_max-k_min+1);
Q_zm0 = eye(nc(1));

if dm(4) == 1
  cost_fcn = @(x) norm(discretize(x) - x, 'fro');
else
  cost_fcn = @(x) sum(sum(bsxfun(@rdivide, x.^2, max(x.^2,[],2))))/m;
end

if dm(1) == 1
  ZMP = @(x) evrot_NAG_batch(x, 20); 
elseif dm(3) == 1
  ZMP = @(x) evrot_orig(x, 1);   % ZMP = @(x) evrot_adam(x, 20, batch_siz(g));
end
t_run = zeros(1, k_max-k_min+1);

for g = 1:1:numel(nc)
  tic;
  dim = nc(g);    
  V_curr = Utils.normalize_rows(V(:,1:dim));
  %V_curr = V(:,1:dim);
  
  if rm(1)==1   % 'prev_rot'    
    V_curr = V_curr*Q_zm0;
    if dim == k_min   % good initialisation at dim==k_min for 'prev_rot'
      [Q_pr, cntrds] = randorthpts(V_curr, min(10*dim,m), cost_fcn, debug);  %  min(3*dim,m)
      V_curr = V_curr*Q_pr;
    else
      Q_pr = eye(dim);
    end
  end
  
  if rm(2)==1  % 'rand_orth'
    Qt1 = cost_fcn(V_curr);
    if ~(dim == k_min && rm(1)==1)
      [Q_pr, cntrds] = randorthpts(V_curr, min(10*dim,m), cost_fcn, debug);
    else
      Q_pr = eye(dim);
    end
    V_tst = V_curr*Q_pr;
    Qt2 = cost_fcn(V_tst);
    if Qt2 < Qt1
      V_curr = V_tst;
    else
      Q_pr = eye(dim); Qt2 = Qt1;
    end
  end
  
  % Extract centroids
  if ~exist('cntrds','var')
    n_cntr = 1;
    C = eye(dim);
  else
    n_cntr = size(cntrds, 1);
    if n_cntr<30
      C = zeros(dim, dim, n_cntr);
    else
      n_cntr=30;
      C = zeros(dim, dim, n_cntr);
    end
    for i = 1:1:n_cntr
      C(:,:,i) = V_curr(cntrds(i,:),:);
    end
    C = cat(3,C,eye(dim));
  end
  
  Qt3_fi = cost_fcn(V_curr);  
  Q_fi = eye(dim);
  if sum(dm(1:3)) == 0      
    for i = 1:1:size(C,3)
      R0 = loewdin(C(:,:,i)');
      [V_tst, Q2] = yushi_discrete(V_curr, R0, 150);
      if debug
        assert(isequal(V_tst, V_curr*Q2));
      end
      Qt_curr = cost_fcn(V_tst);
      if Qt_curr<Qt3_fi
        Qt3_fi = Qt_curr;
        Q_fi = Q2;
      end
    end  
    Qt3 = Qt3_fi;
    V_curr = V_curr*Q_fi;     
  else    
    % Minimise ZMP cost with gradient descend (do not use the Yu-Shi best result to get more or less independent ZMP results)
    if dm(2) == 1
      [Vr_ZM, Qt3_ZM, theta_fi] = evrot_adam(V_curr, 20, batch_siz(g));
    else
      [Vr_ZM, Qt3_ZM, theta_fi] = ZMP(V_curr);
    end
    
    % k-means initialization (as KP-heuristic in Chan-Schlag-Zien)
    T = kmeans( V_curr, dim, 'Distance', 'cosine', 'Replicates', 5,...
      'Start', cat(3, C(:,:,1:4), eye(dim)) );  % 'sqeuclidean'
    Eref = eye(dim);  % any axis would fit if reflections are allowed
    pair = Eref(T,:);
    Q1 = orth_map(V_curr', pair', false);
    if dm(2) == 1
      [Vr_kp, Qt3_kp, theta_kp] = evrot_adam(V_curr*Q1', 20, batch_siz(g));
    else
      [Vr_kp, Qt3_kp, theta_kp] = ZMP(V_curr*Q1');
    end
    
    Qt3_zm = Qt3_fi;
    Vr_zm = V_curr;
    Q_zm = eye(dim);    
    
    % Combine all three outcomes
    [~, idx_best] = min([Qt3_ZM, Qt3_kp, Qt3_zm]);
    
    switch idx_best
      case 1
        Qt3 = Qt3_ZM; Vr = Vr_ZM;
        Q_fi =  givens_matrix(theta_fi, dim);
      case 2
        Qt3 = Qt3_kp; Vr = Vr_kp;
        Q_fi =  Q1'*givens_matrix(theta_kp, dim);        
      case 3
        Qt3 = Qt3_zm; Vr = Vr_zm;
        Q_fi = Q_zm;
    end
    V_curr = V_curr*Q_fi;
    if debug % && (Qt3_kp<Qt3_ZM)
      assert( max(max(abs(V_curr-Vr))) < 1e-12 );
    end
  end
  
  Vr_out{g} = V_curr;
  Qt_out(g) = Qt3;
  
  if rm(1) == 1
    Q_zm0 = Q_zm0*Q_pr*Q_fi;
    Q_zm0 = [ Q_zm0, zeros(size(Q_zm0,1),1); zeros(1,size(Q_zm0,1)+1)]; Q_zm0(end,end) = 1;
    if debug && dim>2
      tst = [Utils.normalize_rows(V(:,1:dim)), ones(m,1)]*Q_zm0;   
      tst = tst(:,1:end-1);
      assert(max(max(abs(tst - Vr_out{g}))) < 1e-10);
    end
  end
  % fprintf('Number of clusters %d, iteration %d, Quality %f.\n', group_num(g), i, Qt);
  t_run(g) = toc;
end

end

function [Q, cost] = orth_map(A, B, rot)
M = B*A';
sz1 = size(A,1) ;      % dimension of space
[U, ~, V] = svd(M);
S = eye(sz1);
if rot
  sgn = sign(det(U*V'));
  S(end,end) = sgn;
end
Q = U*S*V';
if nargout > 1
  cost = norm(Q*A - B, 'fro');
end
end

function [Q_qr, R_qr] = randorthpts(V, iter_max, cost_fcn, debug)
m = size(V, 1);
ord = randperm(m);
ord = ord(:);
[Q_qr, R_qr] = greedy_basis(V(ord,:,:), iter_max, cost_fcn, debug);

if debug
  test = V(ord,:,:);
  for i = 1:1:size(R_qr,1)
    assert( isequal(test(R_qr(i,:),:,:), V(ord(R_qr(i,:)),:,:)) );
    R_qr(i,:) = ord(R_qr(i,:));
  end
else
  for i = 1:1:size(R_qr,1)
    R_qr(i,:) = ord(R_qr(i,:));
  end
end
end

function [Q_qr, R_qr] = greedy_basis(V, iter_max, cost_fcn, debug)

dim = size(V,2);
m = size(V,1);
Qt0 = cost_fcn(V);

R_qr = zeros(iter_max, dim);
Qt_best = Qt0*ones(iter_max, 1);
Q_qr = eye(dim);
iter = 1;
sim_cumul = zeros(m, 1);
del = false(m, 1);
sim_thr = 0.9986;  % 0.9994 -> 2 deg; 0.9975 -> 4 deg; 0.9986 -> 3 deg;
i_all = 1:1:m;
while ~all(del) && iter<=iter_max
  i_del = [find(del)];
  i_cand = my_setdiff(i_all, i_del);
  [dummy, sim_ord] = min(sim_cumul(i_cand));
  curr = i_cand(sim_ord);
  R0 = zeros(dim, dim);
  i_sel = zeros(1,dim);
  i_sel(1,1) = curr;
  vec_sel = V(curr,:)';
  R0(:,1) = vec_sel;
  sim_sel = V*vec_sel;  % does not consider prev. deletions
  similar = sim_sel>sim_thr;
  sim_cumul = sim_cumul + sim_sel;
  del = del | similar;
  c = abs(sim_sel);
  for j = 2:1:dim
    [min_sim, imin] = min(c);
    R0(:,j) = V(imin,:)';
    i_sel(1,j) = imin;
    assert(isequal(R0(:,any(R0,1)), V(i_sel(i_sel~=0),:)'));
    c = c + abs(V*R0(:,j));
  end
  
  Q_curr = loewdin(R0);
  V_test = V*Q_curr;
  Qt_curr = cost_fcn(V_test);
  if Qt_curr < min(Qt_best)
    Q_qr = Q_curr;
  end
  Qt_best(iter,:) = Qt_curr;
  R_qr(iter,:) = i_sel;
  iter = iter+1;
end
Qt_best = Qt_best(1:iter-1,:);
R_qr = R_qr(1:iter-1,:);
[Qt_best, qt_idx] = sort(Qt_best);
R_qr = R_qr(qt_idx,:);
if debug
  for i = 1:1:size(R_qr,1)
     x = (V(R_qr(i,:),:)*V(R_qr(i,:),:)'); x(1:dim+1:dim*dim) = 0; 
     RdotR = max(max(abs(x)));
     isequal(RdotR,Qt_best(i));
  end
end
end
