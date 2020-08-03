function pst = struct_pst(bus, lin, gen, basmva, coh, f0)

% The below "power flow insertion" is not good because it can only insert 
% an AC power flow solution of PST. Better get a PST structure from here,  
% then transform it to a MATPOWER structure through pst2mp(), then solve 
% either AC or DC power flow (or OPF) in MATPOWER, and then transform back 
% to PST with mp2pst.

%{
tol = 2e-13; iter_max = 30; vmin = 0.5; vmax = 1.5; acc = 1.0;
bus(:,13) = 1;      % dummy nominal voltage
[DUMMY, bus, lin_flw] = evalc('loadflow(bus,lin,tol,iter_max,vmin,vmax,acc,''n'',2);');

lin_flw(:,1) = [];
num_arcs = size(lin_flw,1);
lin_flw = [lin_flw(1:2:num_arcs-1,:),lin_flw(2:2:num_arcs,:)];
lin_flw(:,3) = 0.5*(abs(lin_flw(:,3)) + abs(lin_flw(:,7)));
lin_flw(:,4:end) = [];
lin_flw(:,1:2) = sort(lin_flw(:,1:2), 2, 'ascend');
[frto,idx] = sortrows(lin_flw(:,1:2), [1,2]);
lin_flw = lin_flw(idx,:);   % standartize the order of lines
assert(isequal(lin_flw(:,1:2),lin(:,1:2)));
lin(:,12) = lin_flw(:,3);   % actual line flow
%}

if isempty(bus)
  pst = struct('bus', bus, 'lin', lin, 'gen', gen, 'basmva', basmva, 'coh', coh, 'f0', f0);
  return;
end
if size(bus,2)<13
  bus(:,13) = 1;   %dummy nominal voltage
end
if size(lin,2)<12
  lin(:,12) = eps;   %dummy line flow
end
pst = struct('bus', bus, 'lin', lin, 'gen', gen, 'basmva', basmva, 'coh', coh, 'f0', f0);
end