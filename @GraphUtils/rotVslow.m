function [Z_out, J_out, t_run] = rotVslow(V, mode, dm, batch_siz, k_min)

assert(sum(dm)<=1 && size(dm,1)==1  && size(dm,2)==3);

% Add path to Zelnik-Perona gradient descent code and kabsch
old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath(fullfile(fileparts(fileparts(which(mfilename))), 'third_party'));   %kabsch
addpath(fullfile(fileparts(fileparts(which(mfilename))), 'third_party', 'zp_clustering'));

k_max = size(V, 2);
if nargin < 4
  k_min = 2;
end
nc = k_min:k_max;
m = size(V,1);
Z_ys = cell(1, k_max-k_min+1);
J_ys = NaN(1, k_max-k_min+1);
t_ys = zeros(1, k_max-k_min+1);
cost_fcn = @(x) sum(sum(bsxfun(@rdivide, x.^2, max(x.^2,[],2))))/m;   % Zelnik-Manor cost function

if length(mode)>=5 && strcmp(mode(1:5),'YuShi')
  [X, YS, t_ys] = GraphUtils.nc_ev_sym(V, [1 1], [0 0 0 1], [], k_min);
  for g = 1:1:numel(nc)  
    tic;  
    dim = size(X{g},2);
    V_ini = Utils.normalize_rows(V(:,1:dim));
    R = V_ini(1:dim,:)\X{g}(1:dim,:);
    Z_ys{g} = V(:,1:dim)*R;      
    J0 = cost_fcn(Z_ys{g});   
    [Z_ys{g}, JJ, thetas] = evrot_NAG_batch(Z_ys{g});  %refine!       
    assert(JJ<=J0 || abs(JJ-J0)/min(JJ,J0)<0.0008); 
    J_ys(g) = cost_fcn(Z_ys{g});    
    t_ys(g) = t_ys(g) + toc;
    assert( abs(J_ys(g)-JJ)<1e-6 );
  end
end
if strcmp(mode,'YuShi')
  Z_out = Z_ys;
  J_out = J_ys;
  t_run = t_ys;
  return;
end

% Zelnik-Manor algorithm
Z_zm = cell(1, k_max-k_min+1);
J_zm = NaN(1, k_max-k_min+1);
t_zm = zeros(1, k_max-k_min+1);
if nargin < 3 || isempty(batch_siz)
  batch_siz = round(linspace(round(m/10), floor(m/2), k_max-k_min+1));   %floor(m/2) is the maximal feasible batch size
else
  batch_siz = batch_siz*ones(1,k_max-k_min+1);
end
if dm(1) == 1
  ZMP = @(x) evrot_NAG_batch(x, 20); 
elseif dm(3) == 1
  ZMP = @(x) evrot_orig(x, 1);   %ZMP = @(x) evrot_adam(x, 20, batch_siz(g));
end

for g = 1:1:numel(nc)
  tic;
  dim = nc(g); 
  if g==1
    V_curr = V(:,1:dim);
  else
    V_curr = [Z_zm{g-1}, V(:,dim)];
  end
  if dm(2) == 1
    [Z_zm{g}, J_zm(g), ~] = evrot_adam(V_curr, 20, batch_siz(g));
  else
    [Z_zm{g}, J_zm(g), ~] = ZMP(V_curr);
  end
  t_zm(g) = toc;
end

% Take the best of Yu-Shi and Zelnik-Perona
if length(mode)>=5 && strcmp(mode(1:5),'YuShi') && strcmp(mode(end-1:end),'ZM')
  for g = 1:1:numel(nc)
    if J_ys(g)<J_zm(g)
      Z_out{g} = Z_ys{g};
      J_out(g) = J_ys(g);
      t_run(g) = t_ys(g);
    else
      Z_out{g} = Z_zm{g};
      J_out(g) = J_zm(g);
      t_run(g) = t_zm(g);
    end
  end
else
  Z_out = Z_zm;
  J_out = J_zm;
  t_run = t_zm;
end

