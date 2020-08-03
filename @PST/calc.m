function [delP,delQ,P,Q,conv_flag] = ...
                 calc(nbus,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol)
% Syntax:  [delP,delQ,P,Q,conv_flag] = 
%                calc(nbus,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol)
%
% Purpose: calculate power mismatch and check convergence
%
% Input:   nbus      - total number of buses
%          bus_type  - load_bus(3), gen_bus(2), swing_bus(1)
%          V         - magnitude of bus voltage
%          ang       - angle(rad) of bus voltage
%          Y         - admittance matrix
%          Pg        - real power of generation
%          Qg        - reactive power  of generation
%          Pl        - real power of load
%          Ql        - reactive power of load
%          tol       - a tolerance of computational error
%
% Output:  delP      - real power mismatch
%          delQ      - reactive power mismatch
%          P         - calculated real power
%          Q         - calculated reactive power
%          conv_flag - 0, converged
%                      1, not yet converged
%
% See also:  
%
% Calls:
%
% Call By:   loadflow

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ************************************************************
jay = sqrt(-1);
swing_bus = 1;
gen_bus = 2;
load_bus = 3;
% voltage in rectangular coordinate
V_rect = V'.*(cos(ang')+jay*sin(ang'));  
% bus current injection
cur_inj = Y*V_rect;
% power output based on voltages  
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg' - Pl' - P;
delQ = Qg' - Ql' - Q;

% zero out mismatches on swing bus and generation bus
for i = 1:nbus
  if bus_type(i) == swing_bus
      delP(i) = 0;
      delQ(i) = 0;
    elseif bus_type(i) == gen_bus
      delQ(i) = 0;
  end
end

%  total mismatch
mism = norm(delQ,'inf')+norm(delP,'inf');
if mism > tol,
    conv_flag = 1;
  else
    conv_flag = 0;
end
fprintf('mismatch is %g. \n',mism)
return

