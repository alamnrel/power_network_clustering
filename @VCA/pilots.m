function [plt, plt_qt] = pilots(S_LL, LBL, mode, max_cand)

% 
% Given a load-to-load bus sensitivity matrix S_LL and the labels of each
% load T_LL, return the pilot bus indices (in S_LL) for each zone encoded  
% in T_LL. The mode parameter specifies the selection method.
% 
% If max_cand (4th input) is greater than 1, max_cand pilot candidates per
% zone will be returned, starting from highest to lowest quality (the  
% quality criterion is defined by mode, the 3rd input).
% 
% If 2nd output is requested, the value of quality criterion for each pilot 
% node in the first output is returned as well.
% 
% This file currently implements pilot buses from load-to-load voltage
% senstitivities. A generalisation that would additionaly consider the
% generator-to-load sensitivities  seems to be promising, although combining
% several important properties of pilot buses (pilots sensitive to the control
% generators of their own zones, while insensitive to the control generators
% of the other zones, and able to impose their voltage on the loads in their
% zones etc.) is not so straighforward.
%

if nargin<3
  mode = 'strongest_bus';
end
if nargin<4
  max_cand = 1;
end
bus = S_LL.busrow(:);
S_LL = S_LL.S;
nz = numel(LBL);
plt = cell(1,nz);
plt_qt = cell(1,nz);
m = size(S_LL, 1);  
% (no need to symmetrize S_LL here)

switch mode
  case 'strongest_bus'  % for this see the Sandro Corsi's book on voltage control
    sen_self = S_LL(1:m+1:m*m);
    for k = 1:1:nz
      cur = find(ismember(bus, LBL{k}));
      self = sen_self(cur);
      [val, idx] = sort(self, 'ascend');
      plt{k} = bus( cur( idx( 1:min(max_cand,numel(cur)) ) ) );
      plt_qt{k} = self( idx( 1:min(max_cand,numel(val)) ) );
    end
    
  case 'max_zone_sens'
    self = S_LL(1:m+1:m*m);
    nrm = 1./self(:)';
    S_LL = bsxfun(@mtimes, S_LL, nrm);  % dVi/dVj 
    S_LL(1:m+1:m*m) = 0;
    for k = 1:1:nz
      cur = find(ismember(bus, LBL{k}));
      S_CUR = S_LL(cur,cur);
      self = sum(S_CUR,1);   % cumulative impacts of dVj on kth zone voltages
      [val, idx] = sort(self, 'descend');
      plt{k} = bus( cur( idx( 1:min(max_cand,numel(cur)) ) ) );
      plt_qt{k} = self( idx( 1:min(max_cand,numel(val)) ) );
    end
    
  otherwise
    error('');
end

end