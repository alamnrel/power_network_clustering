function [bus_sol,line_flow] = ...
   LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,display,flag)
% Syntax:    [bus_sol,line_flow] = ...
% LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,display,flag)
% 
% Purpose: solve the load-flow equations of power systems
% 
% Input:    bus       - bus data
%           line      - line data
%           tol       - tolerance for convergence
%           iter_max  - maximum number of iterations
%           vmin      - voltage minimum limit
%           vmax      - voltage maximum limit
%           acc       - acceleration factor
%           display   - 'y', generate load-flow study report
%                        else, no load-flow study report
%           flag      - 1, form new Jacobian every iteration
%                       2, form new Jacobian every other iteration
% 
  
% Output:   bus_sol   - bus solution (see report for the solution format)                         
%           line_flow - line flow solution (see report)
% 
% See also:  
%
% Algorithm: Newton-Raphson method using the polar form of 
%   the equations for P(real power) and Q(reactive power).
%
% Calls:     PST.YYbus, PST.calc, PST.form_jac
%
% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Authors:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ***********************************************************
global bus_int

tt = clock;     % start the total time clock
jay = sqrt(-1);
load_bus = 3;
gen_bus = 2;
swing_bus = 1;
if exist('flag') == 0
  flag = 1;
end
%if flag <1 | flag > 2
if flag <1 
  error('LOADFLOW: flag not recognized')
end
nline = length(line(:,1));     % number of lines
nbus = length(bus(:,1));     % number of buses
% set maximum and minimum voltage
volt_min = vmin*ones(1,nbus);
volt_max = vmax*ones(1,nbus);
% build admittance matrix Y
[Y,nSW,nPV,nPQ,SB] = PST.YYbus(bus,line);  %IT 02Jan2020

% process bus data
bus_no = bus(:,1)';
V = bus(:,2)';
ang = bus(:,3)'*pi/180;
Pg = bus(:,4)';
Qg = bus(:,5)';
Pl = bus(:,6)';
Ql = bus(:,7)';
Gb = bus(:,8)';
Bb = bus(:,9)';
%cyb = (Gb + jay*Bb)';
bus_type = bus(:,10)';

% set up index for Jacobian calculation

%% form PQV_no and PQ_no
PQVptr = 1;     % PQV_no pointer
PQptr = 1;     % PQ_no pointer

for i = 1:nbus,
  if bus_type(i) == load_bus,
    PQV_no(PQVptr) = i;
    PQ_no(PQptr) = i;
    PQptr = PQptr + 1;
    PQVptr = PQVptr + 1;
  elseif bus_type(i) == gen_bus,
    PQV_no(PQVptr) = i;
    PQVptr = PQVptr + 1;
  end
end; %%

% construct angle reduction matrix
%ang_red = zeros(length(PQV_no),nbus);
%for i = 1:length(PQV_no)
%  ang_red(i,PQV_no(i)) = 1;
%end; %
il = length(PQV_no);
ii = [1:1:il];
ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);
% construct voltage reduction matrix
%volt_red = zeros(length(PQ_no),nbus);
%for i = 1:length(PQ_no)
%  volt_red(i,PQ_no(i)) = 1;
%end; %
il = length(PQ_no);
ii = [1:1:il];
volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);

iter = 0;     % initialize iteration counter
% calculate the power mismatch and check convergence
[delP,delQ,P,Q,conv_flag] = ...
             PST.calc(nbus,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol);
%keyboard
st = clock;     % start the iteration time clock
%% start iteration process
while (conv_flag == 1 & iter < iter_max)
iter = iter + 1;
% form the jacobian matrix; use partition Jacobian matrix
%  [Jac11,Jac12,Jac21,Jac22] = ...
%	     PST.form_jac(V,ang,Y,ang_red,volt_red);

% form the jacobian matrix; use full matrix formulation
%  t_s = clock;
   if flag >= 2
      if iter == flag*fix(iter/flag) + 1
       clear Jac
       Jac = PST.form_jac(V,ang,Y,PQV_no,PQ_no);
      end
     else
       clear Jac
       Jac = PST.form_jac(V,ang,Y,PQV_no,PQ_no);
   end

%  disp('Jacobian calculation')
%  etime(clock,t_s)
% reduced mismatch real and reactive power vectors
  red_delP = ang_red*delP;
  red_delQ = volt_red*delQ;
  clear delP delQ

% solve for voltage magnitude and phase angle increments
%  Jacinv = inv(Jac22);
%  clear J22
%  Jac = Jac11-Jac12*Jacinv*Jac21;
%  clear J11
%  red_del = red_delP-Jac12*Jacinv*red_delQ;
%  clear J12 red_delP
%  red_delAng = Jac\red_del;
%  red_delV = Jacinv*(red_delQ-Jac21*red_delAng);
%  clear J22 Jacinv red_delQ
%  t_s = clock;
  temp = Jac\[red_delP; red_delQ];
%  disp('network solution')
%  etime(clock,t_s)
%  clear Jac

% expand solution vectors to all buses
  delAng = temp(1:length(PQV_no),:)'*ang_red;
  delV = temp(length(PQV_no)+1:length(PQV_no)+length(PQ_no),:)'*volt_red;

% update voltage magnitude and phase angle
  V = V + acc*delV;
  V = max(V,volt_min);  % voltage higher than minimum
  V = min(V,volt_max);  % voltage lower than maximum
  ang = ang + acc*delAng;
% calculate the power mismatch and check convergence
%  t_s = clock;
  [delP,delQ,P,Q,conv_flag] =...
             PST.calc(nbus,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol);
%  etime(clock,t_s)
end; %%
ste = clock;   % end the iteration time clock

for i = 1:nbus
  if bus_type(i) == gen_bus,
    Pg(i) = P(i) + Pl(i);
    Qg(i) = Q(i) + Ql(i);
  elseif bus_type(i) == load_bus,
    Pl(i) = Pg(i) - P(i);
    Ql(i) = Qg(i) - Q(i);
  end
end
Pg(SB) = P(SB) + Pl(SB); Qg(SB) = Q(SB) + Ql(SB);

%VV(:) = V(:).*exp(jay*ang(:));  % solution voltage 
VV = V.*exp(jay*ang);   % solution voltage
% calculate the line flows and power losses
for i = 1:nline
  tap_ratio(i,1) = line(i,6);
  if tap_ratio(i,1) == 0,   % this line has no transformer
    tap_ratio(i,1) = 1;
  end
end
phase_shift(:,1) = line(:,7);
tps = tap_ratio.*exp(jay*phase_shift*pi/180);
from_bus = line(:,1);
from_int = bus_int(round(from_bus));
to_bus = line(:,2);
to_int = bus_int(round(to_bus));
r = line(:,3);
rx = line(:,4);
chrg = line(:,5);
z = r + jay*rx;
y = ones(nline,1)./z;

while(0)
  MW_bs = VV(:).*conj(VV(:)).*cyb(:);
  P_bs = real(MW_bs);   % active power sent out by from_bus to ground
  Q_bs = imag(MW_bs);   % reactive power sent out by from_bus to ground
end %while(0)

VV = VV(:); % Joe Chow 2013 1224
Vf = VV(from_int)./tps;  %IT 23May2020
Vt = VV(to_int);  %IT 02Jan2020

MW_s = Vf.*conj((Vf - Vt).*y + Vf.*(jay*chrg/2));  %IT 02Jan2020      
P_s = real(MW_s);     % active power sent out by from_bus to to_bus
Q_s = imag(MW_s);     % reactive power sent out by from_bus to to_bus

MW_r = Vt.*conj((Vt - Vf).*y + Vt.*(jay*chrg/2));  %IT 02Jan2020
P_r = real(MW_r);     % active power received by to_bus from from_bus
Q_r = imag(MW_r);     % reactive power received by to_bus from from_bus

for i = 1:nline
  line_flow(2*i-1:2*i,:) = ...
             [i from_bus(i) to_bus(i) P_s(i) Q_s(i)
		    i to_bus(i) from_bus(i) P_r(i) Q_r(i)];
end
% keyboard
P_loss = sum(P_s) + sum(P_r) ;
Q_loss = sum(Q_s) + sum(Q_r) ;
bus_sol=[bus_no'  V'  ang'*180/pi Pg' Qg' Pl' Ql' Gb' Bb' bus_type'];

% display results
if display == 'y',
  clc
  disp('                             LOAD-FLOW STUDY')
  disp('                    REPORT OF POWER FLOW CALCULATIONS ')
  disp(' ')
  disp(date)
  fprintf('SWING BUS                  : BUS %g \n', SB)
  fprintf('NUMBER OF ITERATIONS       : %g \n', iter)
  fprintf('SOLUTION TIME              : %g sec.\n',etime(ste,st))
  fprintf('TOTAL TIME                 : %g sec.\n',etime(clock,tt))
  fprintf('TOTAL REAL POWER LOSSES    : %g.\n',P_loss)
  fprintf('TOTAL REACTIVE POWER LOSSES: %g.\n\n',Q_loss)
  if conv_flag == 0,
    disp('                                      GENERATION             LOAD')
    disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
    disp(bus_sol(:,1:7))

    disp('                      LINE FLOWS                     ')
    disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
    disp(line_flow)
  end
end; %
if iter > iter_max,
  fprintf('Note: Solution did not converge in %g iterations.\n', iter_max)
end

return

