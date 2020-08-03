function F_init = getDiffusedStart(W, nb_of_anchors_per_class, NcutDiscrete, nStarts, PL)
%
% We use the same procedure for generating the starting points as described in the paper,
% X. Bresson, T. Laurent, D. Uminsky and J.H. von Brecht, "Multiclass Total Variation Clustering", NIPS 2013.
% The following is taken from their code to generate these special starting points via the diffusion procedure.
% Note that we do not use their solution as a starting point to our method.
%

n = size(W,1);
D = sum(W,2);
Lap = spdiags(D,0,n,n) - W;
R = size(PL,2);
F_init = zeros(n, R, nStarts);
% Pick randomly labeled data points
if sum(sum(PL))>0
    indicesLabeledData = find(sum(PL,2)~=0);
    indictorFunctionLabeledData = full(PL(indicesLabeledData,:));
    %F_init = zeros(n,R);
    F_init(indicesLabeledData,:, 1) = indictorFunctionLabeledData;
    F_init(:,:, 1) = diffuse(Lap, F_init(:,:,1));
else    
    %nb_of_anchors_per_class = round(0.01* size(W,1)/ R); 
    if ~exist('NcutDiscrete','var')
        NcutDiscrete =ncutW(W,R);
    end
    [~,C_ncut] = max(NcutDiscrete,[],2);
    for t=1:nStarts
        [indicesLabeledData, indictorFunctionLabeledData] = pick_random_indices_from_classes(C_ncut, R, nb_of_anchors_per_class);
        %F_init = zeros(n,R);
        F_init(indicesLabeledData,:, t) = indictorFunctionLabeledData;
        F_init(:,:, t) = diffuse(Lap, F_init(:,:,t));
    end

end

end


function new_F = diffuse(Lap, F)
% Diffuse each column of F data by doing one step of implicit heat equation:
%                     (I+dt Lap) f_r^{new} = f_r    r=1,...,R

dt=1;
[n,R]=size(F); new_F=zeros(n,R);
md=@(x,type) x+dt*Lap*x; 
for i=1:R,
    f1 = F(:,i);
    [f1, ~] = pcg(md,f1,1e-6,50,[],[],f1);
    new_F(:,i) = f1;
end;

end

function  [idx, indicator_matrix]= pick_random_indices_from_classes(C,R,nb_per_class)
% idx is a column vectors of heights nb_per_class*R. 
% The first nb_per_class entries contains indices randomly picked from class 1, 
% The next nb_per_class  entries contains indices randomly picked from class 2, etc...
rng('default'); rng('shuffle');
idx=zeros(nb_per_class*R,1);
indicator_matrix=zeros(nb_per_class*R,R);

for k=1:R 
  idx_in_class_k= find(C==k);
  ClassSize=length(idx_in_class_k);
  if (ClassSize==0)
      display('one of the class is empty!')
  end
  v = randi(  [1 , ClassSize]  ,  nb_per_class , 1 );
  selected_idx = idx_in_class_k(v);
  idx(  (k-1)*nb_per_class +1 :  k*nb_per_class  ) = selected_idx;
  indicator_matrix((k-1)*nb_per_class +1 :  k*nb_per_class , k)=1;
end

end
