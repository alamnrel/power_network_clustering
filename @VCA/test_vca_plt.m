function [mpw_out, dv_mean, qc_lvl, resnorm] = test_vca_plt(mpw0, mpw, SET_PQG, SET_PLT, param, debug)
% 
% Here mpw0 and mpw are assumed to have the same topology (or at least the 
% same buses and gens), which is for simplicity to save time.
% 
% SET_PQG, SET_PLT should be in MATPOWER external indexing, same as mpw and
% mpw0 (for mpw and mpw0 checks can be implemented, but not for SET_PQG or
% SET_PLT).
% 

if isfield(mpw0,'order') && mpw0.order.state=='i', error(''); end
if isfield(mpw,'order') && mpw.order.state=='i', error(''); end
assert( isequal(mpw0.bus(:,1), mpw.bus(:,1)) );
assert( isequal(mpw0.gen(:,[1,4,5]), mpw.gen(:,[1,4,5])));
assert(issorted(mpw0.gen(:,1)) && all(diff(mpw0.gen(:,1))>0));
assert(issorted(mpw0.bus(:,1)) && all(diff(mpw0.bus(:,1))>0));
assert(issorted(mpw.gen(:,1)) && all(diff(mpw.gen(:,1))>0));
assert(issorted(mpw.bus(:,1)) && all(diff(mpw.bus(:,1))>0));
% assert( isequal(mpw0.branch(:,1:2), mpw.branch(:,1:2)));  % branches in test and reference may be different 

% Initial setup:
if ~isfield(param,'step')
  freq = 0.1;
else
  freq = param.step;
end
if ~isfield(param,'balcost')
  objbal = 1e-4;
else
  objbal = param.balcost;
end
if ~isfield(param, 'dvmetr')
  dvmetr = 'MSE';
else
  dvmetr = param.dvmetr;
end

% Get the sensitivities. Derive all sets from sensitivity sets:
[dVgdVl, dQgdVg, dVldQl] = VCA.csvc_sen(mpw, true);
SET_PQ  = dVgdVl.busrow;
SET_PV  = dVgdVl.buscol;  
idx_lod = ~ismember(SET_PQ, SET_PQG);
idx_plt =  ismember(SET_PQ, SET_PLT);
assert(nnz(idx_plt)==numel(SET_PLT));
idx_pqg =  ismember(SET_PV, SET_PQG);
assert(nnz(idx_pqg)==numel(SET_PQG));
idx_avr = ~ismember(SET_PV, SET_PQG);
SET_PQG = SET_PV(idx_pqg);  % redefine these sets
SET_AVR = SET_PV(idx_avr);
SET_PLT = SET_PQ(idx_plt);
SET_LOD = SET_PQ(idx_lod & ~idx_plt);  % "other loads"
gen_pqg = ismember( mpw.gen(:,1), SET_PQG );
bus_pqg = ismember( mpw.bus(:,1), SET_PQG );
gen_avr = ismember( mpw.gen(:,1), SET_AVR );
bus_avr = ismember( mpw.bus(:,1), SET_AVR );
bus_lod = ismember( mpw.bus(:,1), SET_LOD );
bus_plt = ismember( mpw.bus(:,1), SET_PLT );

% Set other constants:
np = numel(SET_PLT);
nl = numel(SET_LOD);
nc = numel(SET_PQG);

% Set the references and pre-control infos:
v_ref = mpw0.bus(:,8);
p_ref = v_ref(bus_plt);
v_ad = mpw.bus(:,8);

% Set solver options
lsopt.Display = 'off';
lsopt.MaxIterations = 300;
lsopt.OptimalityTolerance = 1e-9;
mpopt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1,...
  'pf.tol', 1e-6, 'pf.nr.max_it', 15);

% Start the control:
iter =  0;
v_cur = v_ad;
p_cur = v_ad(bus_plt);
q_lvl = mpw.gen(gen_pqg,3);
p_err = p_ref - p_cur;
q_inc = zeros(nc,1);
err_avg = mean(abs(p_err(:,end)));
Qpu = mpw.gen(gen_pqg,4)/mpw.baseMVA;
Qpl = mpw.gen(gen_pqg,5)/mpw.baseMVA;
[dVgdVl, dQgdVg, dVldQl] = VCA.csvc_sen(mpw, true);
dVgdVl = dVgdVl.S;
dQgdVg = dQgdVg.S;

% To ensure valid initial conditions: v_lo < v0+ < v_hi
v_lo = min(0.9, 0.99*min( mpw.bus(:,8) ));  
v_hi = max(1.1, 1.01*max( mpw.bus(:,8) )); 

% Make negative reactive power levels comparable with the positive ones (by 
% relating to the same basis). This effectively disables the balancing of 
% negative reactive power levels w.r.t. to the actual (often small) lower
% reactive power limits, but these bottom limits are still taken care of in
% the constraints.
Qpl = -Qpu;

while err_avg(end)>1e-4
  
  % Set up the objective (names as in Conejo SVC):
  Cpqv = dVgdVl(idx_plt,idx_pqg);
  Cp = dQgdVg(idx_pqg,idx_pqg);
  Qp = mpw.gen(gen_pqg,3)/mpw.baseMVA;
  Qx = zeros(nc, 1);
  Qx(Qp>=0) = Qpu(Qp>=0);
  Qx(Qp<0) = Qpl(Qp<0);  % any(Qp<0)
  Cp1 = Cp./repmat(Qx,1,nc);
  Cp2 = circshift(Cp1, [1 0]);
  Qp1 = Qp./Qx;  % reactive power level..
  Qp2 = circshift(Qp1, [1 0]);
  C1 = Cpqv;
  d1 = freq*(p_ref - mpw.bus(bus_plt,8));
  C2 = Cp1 - Cp2;
  d2 = Qp2 - Qp1;  
  
  % PQG gen voltage limits:
  A1u = eye(nc);
  b1u = v_hi*ones(nc,1) - mpw.bus(bus_pqg,8);
  A1l = -eye(nc);
  b1l = mpw.bus(bus_pqg,8) - v_lo*ones(nc,1); 
  
  % Non-pilot loads voltage limits:
  Cpq = dVgdVl(idx_lod&~idx_plt, idx_pqg); 
  A2u = Cpq;
  b2u = v_hi*ones(nl,1) - mpw.bus(bus_lod,8);
  A2l = -Cpq;
  b2l = mpw.bus(bus_lod,8) - v_lo*ones(nl,1);  
  
  % Non-PQG (AVR) gen reactive power limits:
  Cpv = dQgdVg(idx_avr, idx_pqg);  
  A3u = Cpv;
  b3u = mpw.gen(gen_avr,4)/mpw.baseMVA - mpw.gen(gen_avr,3)/mpw.baseMVA;
  A3l = -Cpv;
  b3l = mpw.gen(gen_avr,3)/mpw.baseMVA - mpw.gen(gen_avr,5)/mpw.baseMVA;
  
  % PQG gen reactive power limits:
  A4u = Cp;
  b4u = mpw.gen(gen_pqg,4)/mpw.baseMVA - Qp; 
  A4l = -Cp;
  b4l = Qp - mpw.gen(gen_pqg,5)/mpw.baseMVA;   
  
  % Solve it!
  C = [C1; sqrt(objbal)*C2];
  d = [d1; sqrt(objbal)*d2];
  A = [A1u; A1l; A2u; A2l; A3u; A3l; A4u; A4l];
  b = [b1u; b1l; b2u; b2l; b3u; b3l; b4u; b4l];
  lb = -ones(nc,1);  % "remove" these bounds, regulate control step size by PARAM.STEP aka FREQ
  ub =  ones(nc,1);  % if PARAM.STEP=1, loose lb and ub enable to solve in one step
  % First attempt ordinary least squares, then check constraints
  dVc = C\d;  % for better linear algebra & basic solution direction
  resnorm = (C*dVc-d)'*(C*dVc-d);
  residual = C*dVc-d;
  % BUG IN LSQLIN: http://mathforum.org/kb/message.jspa?messageID=5719189
  if any(dVc<lb)||any(dVc>ub)||any(A*dVc>b)
    resnorm0 = resnorm;
    [dVc, resnorm, residual, exitflag] = lsqlin(C, d, A, b, [], [], lb, ub, zeros(nc,1), lsopt);
    assert(resnorm >= resnorm0);
  end
  mpw.gen(gen_pqg,6) = mpw.gen(gen_pqg,6) + dVc;
  [mpw, success] = runpf(mpw, mpopt);
  if ~success, error(''); end
  
  iter =  [iter,iter(end)+1];
  v_cur = [v_cur, mpw.bus(:,8)];
  p_cur = [p_cur, mpw.bus(bus_plt,8)];
  q_lvl = [q_lvl, mpw.gen(gen_pqg,3)];
  q_rel = bsxfun(@rdivide, q_lvl, Qx*mpw.baseMVA);
  p_err = [p_err, p_ref - p_cur(:,end)];  
  q_inc = [q_inc, q_lvl(:,end) - q_lvl(:,end-1)];
  err_avg = [err_avg, mean(abs(p_err(:,end)))];
  if (numel(err_avg)>4 && abs(err_avg(end)-mean(err_avg(end-4:end)))<1e-7)
    break;
  end
  if iter(end)>70
    break,
  end
end
mpw_out = mpw;
v_ac = mpw_out.bus(:,8);
if debug
  clrs = {rgb('violet'),rgb('Turquoise'),rgb('orange'),rgb('blue'),rgb('yellow'),...
      rgb('gray'),rgb('red'),rgb('DarkGreen'),rgb('magenta'),rgb('LimeGreen'),...
      rgb('maroon'), rgb('azure')};
  tstart = 1; tend = size(p_cur,2);
  vlt_plt = [];
  % figure;
  for i = 1:1:numel(SET_PLT)    
    h = plot([p_ref(i),p_cur(i,:)], '.-', 'Color', clrs{i}, 'LineWidth', 1.5, 'MarkerSize', 5);
    vlt_plt = [vlt_plt,h];
    hold on;    
    plot([tstart,tend+1], [p_ref(i),p_ref(i)], '--', 'Color', clrs{i}, 'LineWidth', 1);
  end
  hold off;
  xlim([tstart,tend+1]);
  set(gca, 'xtick', tstart+1:4:tend+1);
  set(gca, 'xticklabel', cellfun(@num2str, num2cell(tstart:4:tend+1), 'unif', false));
  gridLegend(vlt_plt,2,cellfun(@num2str,num2cell(SET_PLT),'un',0), 'location', 'north');  % , 'Orientation', 'horizontal'
  xlabel('Control steps, [-]', 'FontSize', 14);
  ylabel('Voltage magnitude, [-]', 'FontSize', 14);
end

dv_lod = v_ac(bus_lod|bus_plt)-v_ref(bus_lod|bus_plt);
dv_gen = v_ac(bus_pqg|bus_avr)-v_ref(bus_pqg|bus_avr);
dv_lod_avr = v_ad(bus_lod|bus_plt)-v_ref(bus_lod|bus_plt);

switch dvmetr
  case 'MAE'
    dv_mean(1) = mean( abs(dv_lod) );  
    dv_mean(2) = mean( abs(dv_gen) );  % generator deviation
    dv_mean(3) = mean( abs(dv_lod_avr) );  % no CSVC 
  case 'MSE'
    dv_mean(1) = mean( dv_lod'*dv_lod );  
    dv_mean(2) = mean( dv_gen'*dv_gen );  % generator deviation
    dv_mean(3) = mean( dv_lod_avr'*dv_lod_avr );  % no CSVC     
  case 'RMSE'
    dv_mean(1) = rms(dv_lod);  
    dv_mean(2) = rms(dv_gen);  % generator deviation
    dv_mean(3) = rms(dv_lod_avr);  % no CSVC     
  otherwise
    error('');
end

qc_lvl = mpw_out.gen(gen_pqg,3)./mpw_out.gen(gen_pqg,4);
% qc_perf = std(qc_lvl);
% qc_mean = mean(qc_lvl);
% qc_perf = mean( ((qc_lvl - qc_mean).^2)/qc_mean^2 );

if debug
  x_width = 2*241; y_width = 1.5*161;    
  figure;
  trgt = find(bus_lod|bus_plt);
  trgt0 = trgt;
  [err, ord] = sort(abs(v_ref(trgt)-v_ad(trgt)), 'descend');
  trgt = trgt(ord);
  trgt = trgt(1:29);
  plot(v_ref(trgt), '.-', 'Color', 'r', 'MarkerSize', 8, 'Linewidth', 1.5); hold on;
  plot(v_ad(trgt), '.-', 'Color', 'b', 'MarkerSize', 8, 'Linewidth', 1.5);
  data = [v_ref(trgt);v_ad(trgt)];
  set( gca, 'XTick', 1:nnz(trgt) );
  set( gca, 'XTickLabel', cellfun(@num2str,num2cell(trgt)','unif',false));
  ylim( [0.995*min(data), 1.005*max(data)] );
  xlim([1,nnz(trgt)]);
  xlabel('Bus ID, [-]', 'FontSize', 14);
  ylabel('Bus voltage, [p.u.]', 'FontSize', 14);
  set(gca, 'xgrid', 'on');
  set(gcf, 'pos', [100 100 x_width y_width]);    
  legend('VOLTAGE REFERENCE', 'VOLTAGE AFTER LOAD DISTURBANCE');
   
  figure;
  plot(v_ref(trgt), '.-', 'Color', 'r', 'MarkerSize', 8, 'Linewidth', 1.5); hold on;
  plot(v_ac(trgt), '.-', 'Color', rgb('ForestGreen'), 'MarkerSize', 8, 'Linewidth', 1.5);
  set( gca, 'XTick', 1:nnz(trgt) );
  set( gca, 'XTickLabel', cellfun(@num2str,num2cell(trgt)','unif',false));
  ylim( [0.995*min(data), 1.005*max(data)] );
  xlim([1,nnz(trgt)]);
  xlabel('Bus ID, [-]', 'FontSize', 14);
  ylabel('Bus voltage, [p.u.]', 'FontSize', 14);
  set(gca, 'xgrid', 'on');  
  set(gcf, 'pos', [100 100 x_width y_width]);    
  legend('VOLTAGE REFERENCE', 'VOLTAGE AFTER SVC');  
  
  figure;
  plot(abs(v_ref(trgt0)-v_ad(trgt0)), '.-', 'Color', 'b', 'MarkerSize', 8, 'Linewidth', 1.5); hold on;
  plot(abs(v_ref(trgt0)-v_ac(trgt0)), '.-', 'Color', rgb('ForestGreen'), 'MarkerSize', 8, 'Linewidth', 1.5);
  xtick = linspace(1,nnz(trgt0),nnz(trgt0));
  set( gca, 'XTick', xtick );
  xlim([1,nnz(trgt0)]);
  xlabel('Bus ID, [-]', 'FontSize', 14);
  ylabel('\DeltaV, [p.u.]', 'FontSize', 14);
  set(gca, 'xgrid', 'on'); 
  set(gcf, 'pos', [100 100 x_width y_width]);    
  legend('\DeltaV_{ad}', '\DeltaV_{ac}');    
  close;
  close;  
  close;  
end

end

