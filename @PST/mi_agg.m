function [n_bus,n_line,nmac_con] = ...
                       mi_agg(bus,line,area,nmach_a,basemva)
% Syntax   : [n_bus,n_line,nmac_con] = ...
%                      mi_agg(bus,line,area,nmach_a,basemva)
% 
% Purpose  : To aggregate coherent machines using the
%            inertia aggregation method. A solved loadflow
%            input is required.
%
% Input    : bus      - bus data
%            line     - line data
%            area     - matrix of coherent machines
%            nmach_a  - vector of number of machines in each coherent area                         
%            basemva  - base mva (optional)
% Output   : n_bus    - new system bus data
%            n_line   - new line data
%            nmac_con - aggregate generator data
% See also : podmore, s_coh3
%
% Calls    : 
%
% Call by  : 
%
% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 1.0
% Author   : Joe H. Chow (modified slightly by I. Tyuryukanov + extra correctness checks)
% Date     : March 11, 1992 (modified 05/94 to include exciter aggregation)
%                            
PST.pst_var;
jay = sqrt(-1);
if nargin == 4
  basemva = basmva;        % from global variable
end

n_bus  = bus;              % create new bus, line, 
n_line = line;             %  and machine data
n_mac  = mac_con;
[nrow,ncol] = size(mac_con);
if ncol <= 21,
  n_mac = [n_mac ones(nrow,23-ncol)];
end

% Solve input loadflow (this also reduces n_bus to 10 columns)
tol  = 1e-9;    % tolerance for convergence
iter_max = 300; % maximum number of iterations
vmin = 0.5;     % voltage minimum
vmax = 1.5;     % voltage maximum
acc  = 1.0;     % acceleration factor
[DUMMY,n_bus,lin_flw] = evalc('PST.LLoadflow(n_bus,n_line,tol,iter_max,vmin,vmax,acc,''n'',2);');
vlt0 = n_bus(:,1:3); nb0 = size(n_bus,1); 

% set up internal bus & generator index vectors
nbus = length(n_bus(:,1));
bus_int = zeros(round(max(n_bus(:,1))),1);
for i = 1:nbus
  bus_int(n_bus(i,1)) = i;
end
tot_mac  = length(mac_con(:,1)); % total number of machines
mac_int = zeros(round(max(mac_con(:,1))),1);
for i = 1:tot_mac
  mac_int(mac_con(i,1)) = i;
end

bus_vol  = n_bus(:,2);           % system bus voltages
bus_ang  = n_bus(:,3);           % system bus angles
num_area = length(area(:,1));    % number of coherent areas

nmac_con = [];
FLW_OUT = []; FLW_BCK = []; AGG_BUS = []; TRM_BUS = [];
Iagg = zeros(num_area,2);
for r=1:num_area              % cycle thru all areas
  if nmach_a(r) == 1          % single machine area
    nmac_con = [nmac_con; n_mac(mac_int(area(r,1)),:)];
    n_mac(mac_int(area(r,1)),1) = 0;  % set machine # to 0
  else                        % areas with more than 1 machine
    num_mach = nmach_a(r);    % # of coherent machines
    com_bus  = max(n_bus(:,1)) + 1;    % common bus
    AGG_BUS = [AGG_BUS; com_bus];
    mac_list  = mac_int(area(r,1:nmach_a(r)))';   %coherent machine numbers                                
    nlines   = length(n_line(:,1));    %# of lines
    bus_list = bus_int(mac_con(mac_list,2)); 
               % coherent machines bus numbers. It is okay
               % to have identical buses in bus_list
    bus_type = 2;
    if any( n_bus(bus_list,10)==1 )
      bus_type = 1;
    elseif all( n_bus(bus_list,10)==3 )  
      bus_type = 3;
    end
    % inertia weighted aggregate machine internal voltage and angle
    m = mac_con(mac_list,3).*mac_con(mac_list,16)/basemva;
    ma = sum(m);
    Pg = n_bus(bus_list,4).*mac_con(mac_list,22);
    Qg = n_bus(bus_list,5).*mac_con(mac_list,23);
    vb = n_bus(bus_list,2);
    angb = n_bus(bus_list,3)*pi/180;
    vbx = vb.*exp(jay*angb);
    xdp = basemva*mac_con(mac_list,7)./mac_con(mac_list,3);    
    %xdp = xdp*1.9;  %dbg    
    % current and voltage computations include possibilities of multiple
    % generator on the same bus
    int_cur = conj((Pg+jay*Qg)./vbx);
    int_vol = vbx + jay*xdp.*int_cur;
    Iagg(r,2) = sum(int_cur);
    %Total current preservation based aggregate voltage 
    int_pwr = int_vol.*conj(int_cur);    
    vlt_agg = sum(int_pwr)/conj(sum(int_cur));
    mag_cbus = abs(vlt_agg); 
    ang_cbus = angle(vlt_agg); 
    %Inertia-weighted aggregate voltage as in J.H. Chow (Iout won't match)
    %{
    mag_cbus = abs(int_vol)'*m/ma;    %OR abs(m'*int_vol/ma)
    ang_cbus = angle(int_vol)'*m/ma;  %OR angle(m'*int_vol/ma)    
    %}
    vg = abs(int_vol);
    delta = angle(int_vol);
    tapa = mag_cbus./vg;
    phia = ang_cbus*ones(size(angb))-delta;
    xdpa = xdp;
    v1 = mag_cbus*exp(jay*ang_cbus)*ones(size(vg));
    v2 = vb.*exp(jay*angb);
    [s1,s2] = PST.LLine_PQ(v1,v2,zeros(size(vg)),xdpa,zeros(size(vg)),tapa,phia*180/pi);
    assert( mean(abs(Pg+real(s2)))<1e-6 );  %"back" flow =!= generated power
    assert( mean(abs(Qg+imag(s2)))<1e-6 );  %"back" flow =!= generated power    
    FLW_OUT = [FLW_OUT;[com_bus*ones(num_mach,1),n_bus(bus_list,1),real(s1),imag(s1)]];
    FLW_BCK = [FLW_BCK;[n_bus(bus_list,1),com_bus*ones(num_mach,1),real(s2),imag(s2)]];
    xdpe = 0;  
    for i = 1:num_mach      
      xdpe = xdpe + 1/xdpa(i);
    end
    xdpe = 1/xdpe;    
    add_bus(1,1)   = com_bus;         % add common bus to tapa,phia*180/pi
    add_bus(1,2)   = mag_cbus;        % common bus voltage
    add_bus(1,3)   = ang_cbus*180/pi; % angle
    add_bus(1,4:9) = zeros(1,6);
    add_bus(1,10)  = 3;               % bus type
    n_bus = [n_bus; add_bus];
    n_bus(bus_list,4:5)= zeros(num_mach,2);   % remove P & Q 
                                      
    for j = 1:num_mach
      n_bus(bus_list(j),6:7)= n_bus(bus_list(j),6:7) - ...
        [Pg(j)+real(s2(j)) Qg(j)+imag(s2(j))];   % adjust load
    end
    n_bus(bus_list,10) = ones(size(bus_list))*3;   % change load type

    rl = nlines + num_mach;
    n_line(nlines+1:rl,1) = ones(size(bus_list))*com_bus;   % from bus                                          
    n_line(nlines+1:rl,2) = n_bus(bus_list,1);   % to bus
    n_line(nlines+1:rl,4) = xdpa;
    n_line(nlines+1:rl,6) = tapa; 
    n_line(nlines+1:rl,7) = phia*180/pi;
    xdeq = xdpe;  %xde0 = 1/sum(ones(num_mach,1)./xdp)
    term_bus = max(n_bus(:,1))+1;     % new terminal bus # 
    TRM_BUS = [TRM_BUS; term_bus];        
    com_vol = mag_cbus*exp(jay*ang_cbus);
    new_cur = conj(sum(s1)/com_vol);
    term_vol = com_vol - jay*xdeq*new_cur;    
    add_bus(1,1)   = term_bus;
    add_bus(1,2)   = abs(term_vol);
    add_bus(1,3)   = angle(term_vol)*180/pi;
    add_bus(1,4)   = real(term_vol * conj(new_cur));
    add_bus(1,5)   = imag(term_vol * conj(new_cur)); 
    add_bus(1,6:9) = zeros(1,4); 
    add_bus(1,10)  = bus_type;
    n_bus     = [n_bus; add_bus];
    Iagg(r,1) = term_bus;
    n_line(rl+1,1)      = term_bus;
    n_line(rl+1,2)      = com_bus;
    n_line(rl+1,3)      = 0.0;
    n_line(rl+1,4)      = -xdeq;
    n_line(rl+1,5:7)    = zeros(1,3);

    %  perform machine aggregation
    agg_mac = PST.eqgen(n_mac,mac_list,basemva, term_bus,area(r,1)); 
    if agg_mac(1,8) ~= 0
      agg_mac(1,8) = xdpe;
    else
      agg_mac(1,7) = xdpe;
    end
    nmac_con = [nmac_con; agg_mac];
    n_mac(mac_list,1) = zeros(num_mach,1);  %set machine # to zero
  end
end
row_del = sum(Iagg,2)==0;
Iagg(row_del,:) = [];

%  organize machine data
for i = 1:tot_mac
  if n_mac(i,1) ~= 0
    nmac_con = [nmac_con; n_mac(i,:)];
  end
end



%a. Check loadflow voltages
[DUMMY,n_bus,lin_flw] = evalc('PST.LLoadflow(n_bus,n_line,tol,iter_max,vmin,vmax,acc,''n'',2);');
assert( norm(vlt0(1:nb0,1)-n_bus(1:nb0,1),'fro')<1e-9 );
assert( norm(vlt0(1:nb0,2)-n_bus(1:nb0,2),'fro')/nb0<1e-7 );
assert( norm(vlt0(1:nb0,3)-n_bus(1:nb0,3),'fro')/nb0<1e-5 );
%b. Check aggregated power flows
FLW_OUT = sortrows(FLW_OUT,[1,2]);
FLW_BCK = sortrows(FLW_BCK,[2,1]);
NEW_OUT = NaN(size(FLW_OUT));
NEW_BCK = NaN(size(FLW_BCK));
pos = 1;
for i = 1:1:numel(AGG_BUS)  
  row = lin_flw(:,2) == AGG_BUS(i) & lin_flw(:,3) < AGG_BUS(i);
  NEW_OUT(pos:pos+nnz(row)-1,:) = lin_flw(row,2:end);
  row = lin_flw(:,3) == AGG_BUS(i) & lin_flw(:,2) < AGG_BUS(i);
  NEW_BCK(pos:pos+nnz(row)-1,:) = lin_flw(row,2:end);    
  pos = pos+nnz(row);
end
NEW_OUT = sortrows(NEW_OUT,[1,2]);
NEW_BCK = sortrows(NEW_BCK,[2,1]);
assert(norm(FLW_OUT(:,1:2)-NEW_OUT(:,1:2),'fro')/nb0<1e-12);
assert(norm(FLW_BCK(:,1:2)-NEW_BCK(:,1:2),'fro')/nb0<1e-12);
assert(norm(FLW_OUT(:,3:4)-NEW_OUT(:,3:4),'fro')/nb0<1e-6);
assert(norm(FLW_BCK(:,3:4)-NEW_BCK(:,3:4),'fro')/nb0<1e-6);
%c. Check aggregated currents (Ie = sum(Ik))
Iagg = sortrows(Iagg,1);
ng = numel(TRM_BUS);
Iout = zeros(ng,2);
for k = 1:1:ng
  row = find(n_bus(:,1)==TRM_BUS(k));
  Pout = n_bus(row,4);
  Qout = n_bus(row,5);
  Vout = n_bus(row,2)*exp(1i*n_bus(row,3)*pi/180);
  Iout(k,:) = [n_bus(row,1), conj((Pout + 1i*Qout)/Vout)];  
end
Iout = sortrows(Iout,1);
assert(norm(Iout(:,1)-Iagg(:,1),'fro')/ng<1e-12);
assert(norm(Iout(:,2)-Iagg(:,2),'fro')/ng<1e-6);
