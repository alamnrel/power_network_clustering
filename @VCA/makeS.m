function S_out = makeS( mpw, mode, incr)

% 
% This function returns power flow sensitivities of voltage angles and 
% magnitudes in SET_TO with respect to reactive power injections from buses 
% in SET_FROM. The sensitivity matrix is computed by perturbing the powers 
% of SET_FROM buses by incr MVar disturbance and rerunning the load flow. 
% 
% The MODE input also defines SET_FROM and SET_TO. SET_FROM can be the set 
% of all generator buses or the set of load buses. SET_TO can only be set 
% of all load buses.
%
% The sets SET_FROM and SET_TO can be filtered with certain criteria (e.g.
% exclude all TO buses below certain voltage or nominal power, or all FROM
% buses below certain generation reactive power limit), but this is not yet
% implemented. 
% Filters are introduced to compute less power flows and reduce the run
% time, but they are not strictly needed.
% 
% Returned busrow and buscol fields of S_out use the indexing of the input
% MPW structure (i.e., if MPW is in internal indexing, buscol and busrow 
% will be in internal indexing too to be consistent with MPW.bus(:,1)).
%

if ~isfield(mpw,'order')
  input_state = 'e';
else
  input_state = mpw.order.state;
end
if input_state=='e'
  mpw = ext2int(mpw);
  i2e = mpw.order.bus.i2e;
end
if nargin<2
  mode = {'G'};
end
switch mode{1}
  case 'G'
    SET_FROM = mpw.gen(:,1);  % just any generator (max. possible set), as we assess possible impacts, and to narrow down use FILT_FROM
  case 'L'
    SET_FROM = mpw.bus(mpw.bus(:,2)==1, 1);  % just any PQ-bus (max. possible set), to reduce it use FILT_TO
  otherwise
    error('');
end
if mode{1}=='L'
  SET_TO = SET_FROM;  % as it should be square matrix
else
  SET_TO = setdiff( mpw.bus(:,1), SET_FROM ); 
end
if nargin<3 
  incr = mpw.baseMVA;   % [MW/MVar]...
end   

% Decide on input variables
col_in = 4;  % Qlod is 4 in mpw.bus
col_out = 8;  % Vm is 9 in mpw.bus     
bus0 = mpw.bus;
rows = mpw.bus(:,1);
m_in = numel(SET_FROM);
m = size(mpw.bus, 1);

% Make S by perturbed power flows (this also considers actual operating 
% conditions more precisely than the FDLF-assumption dV->dQ)
mpOpt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1,...
  'pf.tol', 1e-10, 'pf.nr.max_it', 25);
S = NaN(m,m_in);
out0 = bus0(:,col_out); 
SET_FROM = sort(SET_FROM, 'ascend');
for k = 1:1:m_in
  row = find( rows==SET_FROM(k) );  % which row in mpw.bus this bus index is
  bus = bus0;
  if bus(row,2)~=1   % if PV or SLACK
    bus(row,2) = 1;   % set as PQ bus
  end  
  bus(row, col_in) = bus(row, col_in) - incr;   % add Qgen for default-positive dV/dQ
  mpw.bus = bus;
  [~, bus, ~, ~, success] = runpf(mpw, mpOpt);
  if ~success
    error([mfilename,':PowerFlowNotConverged'],...
      '[%s] Power flow not converged, try a smaller power disturbance.',...
      mfilename);
  end
  out = bus(:,col_out);
  dout = out - out0;
  S(:,k) = dout/incr*mpw.baseMVA; 
end

idx_to = ismember(rows, SET_TO);  
S = S(idx_to,:); 
S_out.S = S;
S_out.buscol = SET_FROM;
S_out.busrow = rows(idx_to);

if input_state=='e'
  S_out.buscol = i2e(S_out.buscol);
  S_out.busrow = i2e(S_out.busrow);
end

end



%{
% An unjustified move to convert PV surrounded by other PV to PQ to avoid
% zero voltage sensitivities for such buses:
nbr_cur = find( adj(:,row) );   % doing the same for loads
is_nbr_pv = ismember(nbr_cur,  mpw.bus(mpw.bus(:,2)==2 | mpw.bus(:,2)==3, 1)  ); -> the VOLTAGE_REGULATED buses (not just gens)
if all(is_nbr_pv)
  bus(nbr_cur,2) = 1; % turn the nbr PV buses to PQ to avoid all-0 sensitivities
end
%}