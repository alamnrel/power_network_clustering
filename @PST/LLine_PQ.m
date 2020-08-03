function [S1,S2] = LLine_PQ(V1,V2,R,X,B,tap,phi)
% Syntax:   [S1,S2] = LLine_PQ(V1,V2,R,X,B,tap,phi) 
%
% Purpose:  Compute line flows. Inputs can be vectors.
%
% Input:    V1        - from bus complex voltage
%           V2        - to bus complex voltage
%           R         - line resistance
%           X         - line reactance
%           B         - line charging
%           tap       - tap ratio
%           phi       - phase shifter angle in degrees
% Output:   S1        - complex power injection at from bus
%           S2        - complex power injection at to bus
%
% See also:  
%
% Algorithm: 
%
% Calls:     
%
% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Authors:   Joe H. Chow
% Date:      March 1992
%
% ***********************************************************
jay = sqrt(-1);
nline = length(V1);
for i = 1:nline
  if tap(i) == 0
    tap(i) = 1;
  end
end
tps = tap.*exp(jay*phi*pi/180);

z = R + jay*X;
y = ones(nline,1)./z;

% This is now according to the standard MATPOWER/PST notation 
% for tps (v1->v1/tps->v2) (IT 23May2020):
V1p = V1./tps;  %from-vlt after ideal f->t PST (IT 23May2020)
S1 = V1p.*conj((V1p - V2).*y + jay*B/2.*V1p);  %IT 02Jan2020       
S2 = V2.*conj((V2 - V1p).*y + jay*B/2.*V2);  %IT 02Jan2020
