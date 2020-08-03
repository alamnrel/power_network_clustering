function [mac_ord, tim_sim, out_ang, out_spd] = SIM_PST(pst_str, mac_ref, bus_flt, bus_adj, t_flt, t_clr, t_end)

jay = sqrt(-1);

% Save & set up global variables
pststate = PST.savepst();
clearvars -global;
PST.pst_var 	     % set up global variable 
mac_con = pst_str.gen;
line = pst_str.lin;
bus = pst_str.bus;
basrad = 2*pi*pst_str.f0;   % system frequency is 60 Hz
basmva = pst_str.basmva;    % 100 MVA base
m = size(mac_con, 1);
line = line(:,1:7);

% Solve loadflow
tol = 1e-9;     % tolerance for convergence
iter_max = 30;  % maximum number of iterations
vmin = 0.5;     % voltage minimum
vmax = 1.5;     % voltage maximum
acc = 1.0;      % acceleration factor
[DUMMY, bus, line_flw] = evalc('PST.LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,''n'',2);');
syn_ref = mac_ref;    % synchronous reference frame

% simulation
t_switch(1) = 0;      % all time in second+s, start time
t_switch(2) = t_flt;  % time to apply fault
t_switch(3) = t_clr;  % time to clear fault, 6 cycles
t_switch(4) = t_end;  % end time
h = 0.01;             % integration stepsize
stepsize = h;
k_switch(1) = round((t_switch(2)-t_switch(1))/stepsize)+1;
k_switch(2) = round((t_switch(3)-t_switch(1))/stepsize)+1;
k_switch(3) = round((t_switch(4)-t_switch(1))/stepsize)+1;

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
[Y_red,V_rec] = PST.RED_YYbus(bus,line);   % pre-fault admittance matrix

% create bus matrix with short circuit on bus 29
bus_f = bus;
row_f = bus_f(:,1)==bus_flt;
bus_f(row_f,6) = 100000.;
[Y_red_f,V_rec_f] = PST.RED_YYbus(bus_f,line);   % fault-on admittance matrix

% remove line bus_flt-bus_adj
line_pf = line;
row_flt = ismember(line(:,1),[bus_flt;bus_adj])&ismember(line(:,2),[bus_flt;bus_adj]);
line_pf(end+1,:) = [ bus_flt bus_adj -line(row_flt,3:5) line(row_flt,6:7) ];   % remove line bus_flt-bus_adj
[Y_red_pf,V_rec_pf] = PST.RED_YYbus(bus,line_pf);   % post-fault admittance matrix

% step 2: initialization
flag = 0;
f = PST.mac_em(0,1,bus,flag);   % all machine electromechanical model

% step 3: perform a predictor-corrector integration
for k = 1:k_switch(3)
  % step 3a: network solution
  %mach_ref(k) = 0;
  mach_ref(k) = mac_ang(syn_ref,k);
  flag = 1;
  f = PST.mac_em(0,k,bus,flag);   % network-machine interface
  psi = psi_re(:,k) + jay*psi_im(:,k);
  if k >= k_switch(2)
    cur = Y_red_pf*psi; % network solution
    %    bus_v(:,k) = V_rec_pf*psi; % bus voltage reconstruction
  elseif k >= k_switch(1)
    cur = Y_red_f*psi; % network solution
    %    bus_v(:,k) = V_rec_f*psi;   % bus voltage reconstruction
  else
    cur = Y_red*psi; % network solution
    %    bus_v(:,k) = V_rec*psi;   % bus voltage reconstruction
  end
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
  % step 3b: compute dynamics and integrate
  flag = 2;
  pmech(:,k) = pmech(:,1); % constant mechanical input power
  f = PST.mac_em(0,k,bus,flag); % dynamics calculation
  if k ~=k_switch(3)
    j = k+1;
    % following statements are predictor steps
    mac_ang(:,j) = mac_ang(:,k) + h*dmac_ang(:,k);
    mac_spd(:,j) = mac_spd(:,k) + h*dmac_spd(:,k);
    edprime(:,j) = edprime(:,k) + h*dedprime(:,k);
    eqprime(:,j) = eqprime(:,k) + h*deqprime(:,k);
    flag = 1;
    
    %mach_ref(j) = 0;
    mach_ref(j) = mac_ang(syn_ref,j);
    f = PST.mac_em(0,j,bus,flag);
    psi = psi_re(:,j) + jay*psi_im(:,j);
    if k >= k_switch(2)
      cur = Y_red_pf*psi;
    elseif k >= k_switch(1)
      cur = Y_red_f*psi;
    else
      cur = Y_red*psi;
    end
    cur_re(:,j) = real(cur); cur_im(:,j) = imag(cur);
    pmech(:,j) = pmech(:,k);
    flag = 2;
    f = PST.mac_em(0,j,bus,flag);
    % following statements are corrector steps
    mac_ang(:,j) = mac_ang(:,k) + h*(dmac_ang(:,k)+dmac_ang(:,j))/2.;      
    mac_spd(:,j) = mac_spd(:,k) + h*(dmac_spd(:,k)+dmac_spd(:,j))/2.;      
    edprime(:,j) = edprime(:,k) + h*(dedprime(:,k)+dedprime(:,j))/2.;      
    eqprime(:,j) = eqprime(:,k) + h*(deqprime(:,k)+deqprime(:,j))/2.;      
  end
end

mac_ord = mac_con(:,1);
t = [0:stepsize:t_switch(4)]';  % time
tim_sim = t;

%Restore the original PST state
out_ang = mac_ang;
out_spd = mac_spd;
PST.loadpst(pststate);


    



