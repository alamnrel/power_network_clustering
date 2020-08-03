function [dVldVg_out, dQgdVg_out, dVldQl_out, dVldQg_out, JAC_out, Ybus_out] = csvc_sen(mpw, debug)

% 
% This function returns the following matrices needed to simulate CSVC:
% 1) dVload/dVgen
% 2) dQgen/dVgen 
% 3) dVload/dQload
% The first of this matrices (dVload/dVgen) can also be used to characterize
% the coupling between loads and generators.
% The remaining output matrices JAC_out and Ybus_out can be requested for
% completeness.
% These sensitivies are analytically correct and match with power flow
% simulations... However, may be the most practical way to obtain
% sensitivities is to "perturb and observe" the load flows.
% 
% The buses in SET_FROM and SET_TO should be expressed with the MATPOWER   
% internal indexing (see ext2int(), int2ext()).
% 
% Returned busrow and buscol fields of S_out use the MATPOWER internal  
% indexing (i.e, as if buses were network nodes consequently numbered from
% one, with isolated buses excluded). To convert to the original bus 
% indices, use the mapping given in mpw.order.bus.i2e.
% 

epsilon = eps('double');
if ~isfield(mpw, 'order')
  input_state = 'e';
else
  input_state = mpw.order.state;
end
if input_state=='e'
  mpw = ext2int(mpw);
  i2e = mpw.order.bus.i2e;
end
if nargin<2
  debug = false;
end
m = size(mpw.bus, 1);
SET_PV = find(mpw.bus(:,2)==2);  % row index..
SET_PQ = find(mpw.bus(:,2)==1);
SET_SW = find(mpw.bus(:,2)==3);  

% Get dQ/dV
[JAC, Ybus] = makeJac(mpw, true);
SW_EXL = sort([SET_PV;SET_PQ], 'ascend'); 
PQ_EXL = sort([SET_PV;SET_SW], 'ascend');
H = JAC(SW_EXL, SW_EXL);  % excl. p_sw and theta_sw
K = JAC(m+1:end,SW_EXL);  % excl. theta_sw
N = JAC(SW_EXL,m+1:end);  % excl. p_sw
L = JAC(m+1:end,m+1:end);
dQdV = L - K/H*N;

% Get dVload/dQload
Bll = dQdV(SET_PQ,SET_PQ);
dVldQl = full(inv(Bll)); 

% Get dVload/dVgen
dQldVg = dQdV(SET_PQ,PQ_EXL);
dVldVg = -dVldQl*dQldVg; 

% Get dVgen/dQgen
Bgg = dQdV(PQ_EXL,PQ_EXL);
Bgl = dQdV(PQ_EXL,SET_PQ);
Blg = dQdV(SET_PQ,PQ_EXL);
dQgdVg = (Bgg - Bgl*dVldQl*Blg);

% % Optionally enforce symmetry 
% dVldQl = (dVldQl + dVldQl')/2;

% Get dVlod/dQgen (Var control space)
SET_PQ0 = SET_PQ;
dVldQg = zeros(numel(SET_PQ),numel(PQ_EXL));
for j = 1:1:numel(PQ_EXL)  
  SET_PQ = [SET_PQ0;PQ_EXL(j)];
  dQxdVx = dQdV(SET_PQ,SET_PQ);
  dVxdQx = inv(dQxdVx);
  dVldQg(:,j) = dVxdQx(1:numel(SET_PQ0),numel(SET_PQ0)+1);
end
SET_PQ = SET_PQ0;

% Buses in output matrices (what 1. rows and 2. columns are)
dVldVg_out.S = dVldVg;
dVldVg_out.buscol = mpw.bus(PQ_EXL, 1);
dVldVg_out.busrow = mpw.bus(SET_PQ, 1); 
dQgdVg_out.S = dQgdVg;
dQgdVg_out.buscol = mpw.bus(PQ_EXL, 1); 
dQgdVg_out.busrow = mpw.bus(PQ_EXL, 1);
dVldQl_out.S = dVldQl;
dVldQl_out.buscol = mpw.bus(SET_PQ, 1);
dVldQl_out.busrow = mpw.bus(SET_PQ, 1);
dVldQg_out.S = dVldQg;
dVldQg_out.buscol = mpw.bus(PQ_EXL, 1);
dVldQg_out.busrow = mpw.bus(SET_PQ, 1);
JAC_out.S = JAC;
JAC_out.buscol = mpw.bus(:,1);
JAC_out.busrow = mpw.bus(:,1);
Ybus_out.S = Ybus;
Ybus_out.buscol = mpw.bus(:,1);
Ybus_out.busrow = mpw.bus(:,1);

% Run some clean-ups of small values:
i_dqdv = abs(dQgdVg*mpw.baseMVA)>0.5;   % for 1 p.u.(!) gen. voltage rise, just 0.5Mvar power change
dQgdVg = dQgdVg.*i_dqdv;
i_dvdv = abs(dVldVg)>1e-3;   % for 1 p.u.(!) gen. voltage rise, just 0.001 p.u. lod. voltage change
dVldVg = dVldVg.*i_dvdv;

% Run some tests
if debug  
  opt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1,...
    'pf.tol', 1e-6, 'pf.nr.max_it', 15);  
  bus0 = mpw.bus;  
  gen0 = mpw.gen;
  mpw0 = mpw;
  
  % Check dVldQg
  [S_GL] = VCA.makeS(mpw, {'G'}, 1); 
  [S_GL.busrow,idx] = sort(S_GL.busrow);
   S_GL.S = S_GL.S(idx,:);
  [S_GL.buscol,idx] = sort(S_GL.buscol);
   S_GL.S = S_GL.S(:,idx);  
  [dVldQg_out.busrow,idx] = sort(dVldQg_out.busrow);
   dVldQg_out.S = dVldQg_out.S(idx,:);
  [dVldQg_out.buscol,idx] = sort(dVldQg_out.buscol);
   dVldQg_out.S = dVldQg_out.S(:,idx);     
  assert( abs(mean(mean(dVldQg_out.S - S_GL.S)))<1e-4 );
    
  % Check dQgdVg (small voltage setpoint increments are critical for convergence!)  
  dQgdVg_X = dQgdVg_out;
  for k = 1:1:numel(PQ_EXL)
    mpw = mpw0;
    gen = gen0;
    row = find(mpw.gen(:,1)==mpw.bus(PQ_EXL(k),1));
    gen(row,6) = gen(row,6) + 0.01;  % change voltage setpoint
    mpw.gen = gen;
    [mpw, success] = runpf(mpw, opt);
    gen = mpw.gen;
    if ~success, continue; end
    dqdv = (gen(:,3) - gen0(:,3))/mpw.baseMVA;  
    i_dqdv = abs(dqdv*100*mpw.baseMVA)>0.5;   % for 1 p.u.(!) voltage rise, just 0.5Mvar power change
    dqdv = dqdv.*i_dqdv;
    col = dQgdVg_out.buscol==mpw0.bus(PQ_EXL(k),1);  % mpw goes to ext. indexing after runpf()
    err = abs(dqdv - dQgdVg(:,col)/100);
    xxx = [dqdv, dQgdVg(:,col)/100];
    dQgdVg_X.S(:,col) = dqdv*100;
    if mean(err/max(abs(dQgdVg(:,col)/100))) > 1e-2
      cprintf('r', 'dQgdVg of generator at bus %d is far from the observed one, mean(err)=%.3f.\n', PQ_EXL(k), mean(err));
    end
  end
  dQgdVg_out = dQgdVg_X;
    
  % Check dVldVg - OK  
  dVldVg_X = dVldVg_out;
  for k = 1:1:numel(PQ_EXL)
    mpw = mpw0;
    gen = gen0;
    row = find(gen(:,1)==mpw.bus(PQ_EXL(k),1));
    gen(row,6) = gen(row,6) + 0.01;
    mpw.gen = gen;
    [mpw, success] = runpf(mpw, opt);
    bus = mpw.bus;
    if ~success, continue; end
    dv = bus(SET_PQ,8) - bus0(SET_PQ,8);
    i_dvdv = abs(dv*100)>1e-3;   % for 1 p.u.(!) gen. voltage rise, just 0.001 p.u. lod. voltage change
    dv = dv.*i_dvdv;
    col = dVldVg_out.buscol==mpw0.bus(PQ_EXL(k),1);
    err = abs(dVldVg(:,col)/100 - dv);
    xxx = [dv, dVldVg(:,col)/100];
    dVldVg_X.S(:,col) = dv*100;
    if mean(err/max(abs(dVldVg(:,col)/100))) > 1e-2
      cprintf('r', 'dVldVg of generator at bus %d is far from the observed one, mean(err)=%.3f\n', mpw.bus(PQ_EXL(k)), mean(err));
    end
  end
  dVldVg_out = dVldVg_X;
    
  % Check dVldQl - OK  
  for k = 1:1:numel(SET_PQ)
    mpw = mpw0;
    bus = bus0;
    bus(SET_PQ(k),4) = bus(SET_PQ(k),4) + mpw.baseMVA;  % + 1 p.u.
    mpw.bus = bus;
    [mpw, success] = runpf(mpw, opt);
    bus = mpw.bus;
    if ~success, continue; end
    dv = bus0(SET_PQ,8) - bus(SET_PQ,8);  % [dVdQ_LL(:,SET_LOD(k))'; dv']
    col = dVldQl_out.buscol==mpw0.bus(SET_PQ(k),1);
    xx = dVldQl(:,col) - dv;
    if mean(abs(xx))>3e-2
      cprintf('r', 'dVldQl of load at bus %d is far from the observed one, mean(err)=%.3f.\n', mpw.bus(SET_PQ(k)), mean(abs(xx)) );
    end
  end
end

if input_state=='e'
  dVldVg_out.buscol = i2e(dVldVg_out.buscol);
  dVldVg_out.busrow = i2e(dVldVg_out.busrow);
  dQgdVg_out.buscol = i2e(dQgdVg_out.buscol);
  dQgdVg_out.busrow = i2e(dQgdVg_out.busrow);
  dVldQl_out.buscol = i2e(dVldQl_out.buscol);
  dVldQl_out.busrow = i2e(dVldQl_out.busrow);  
  dVldQg_out.buscol = i2e(dVldQg_out.buscol);
  dVldQg_out.busrow = i2e(dVldQg_out.busrow);  
  JAC_out.buscol = i2e(JAC_out.buscol);
  JAC_out.busrow = i2e(JAC_out.busrow);
  Ybus_out.buscol = i2e(Ybus_out.buscol);
  Ybus_out.busrow = i2e(Ybus_out.busrow);
end

end


