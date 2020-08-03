function [red_mod,map] = zhukovpp(pst,grp,nm,trmbusagg,modefit)
% 
% This file implements Zhukov generator aggregation at generator internal 
% buses (or generator terminal buses, if the third input is set to true).
% The input data should be in the Power System Toolbox (Chow-Rogers) format
% and generators should be represented by the classical second order model.
% 
% Input: 
% pst - power system data structure in the PST format
% grp - groups of machines to be agregated (in rows of a 2D matrix)
% nm  - number of slow modes to fit during the aggregation process
% trmbusagg - if true aggregate machines at terminal buses (Zhukov, Podmore) 
%           - if false aggregate machines at internal modes (JH Chow et al)
% modefit   - if true try to fit nm slow system modes during aggregation
%           - if false just conventionally aggregate the machine buses 
%
% A SIDE NOTE: Aggregating generator internal nodes excludes any chance of 
% connections between  "buses to be aggregated",  incl. those belonging to 
% distinct groups of buses. And it is indeed usually more accurate.
% 
if nargin<4
  modefit = true;
end
nbus = pst.bus; nlin = pst.lin; ngen = pst.gen; basmva = pst.basmva; 
gen0 = pst.gen; lin0 = pst.lin; f0 = pst.f0; w0 = 2*pi*f0;
%Old data files may have round-off issues
nbus(:,1) = round(nbus(:,1)); nbus(:,10) = round(nbus(:,10));
ngen(:,1) = round(ngen(:,1)); ngen(:,2) = round(ngen(:,2)); 
gen0(:,1) = round(gen0(:,1)); gen0(:,2) = round(gen0(:,2)); 
nlin(:,1:2) = round(nlin(:,1:2)); lin0(:,1:2) = round(lin0(:,1:2));  
%Adjust matrix sizes
nlin = nlin(:,1:7); lin0 = lin0(:,1:7);
if size(ngen,2)<22, ngen(:,22:23) = 1; end 
if size(gen0,2)<22, gen0(:,22:23) = 1; end
ngen(:,19) = round(ngen(:,19));
%Initialize generator mapping
m = size(gen0,1);
map = (1:m)';
%Remove single machine groups from aggregation
grp_del = sum(logical(grp),2)<=1;
if ~any(grp(~grp_del,:))
  red_mod.lin = nlin;
  red_mod.bus  = nbus;
  red_mod.gen = ngen;  
  red_mod.basmva = pst.basmva; 
  red_mod.f0 = pst.f0; 
end

%Solve input loadflow
tol  = 1e-9;    % tolerance for convergence
iter_max = 300; % maximum number of iterations
vmin = 0.5;     % voltage minimum
vmax = 1.5;     % voltage maximum
acc  = 1.0;     % acceleration factor
[DUMMY,nbus,lin_flw] = evalc('PST.LLoadflow(nbus,nlin,tol,iter_max,vmin,vmax,acc,''n'',2);');
vlt0 = nbus(:,1:3); nb0 = size(nbus,1); 

%Obtain the original inv(M)*K system matrix
[~, lam0, MKo, Ko, Mo] = PST.pstVslow(pst, size(gen0,1));
assert(max(lam0)<0.0001);
frq0 = sort(sqrt(abs(lam0))/2/pi);

%Finally remove single machine areas
grp(grp_del,:) = [];
ng = size(grp,1);

%Ensure that classical model generators are aggregated
gen_agg = nonzeros(grp);
ind_m = ismember(ngen(:,1),gen_agg);
assert(all(sum(ngen(ind_m,[4:6,8:15,18,20:21]),2)==0));

% Extend bus/line data to internal nodes (to aggregate INTERNAL NODES with Zhukov)
oldgenbus = [];
for k = 1:1:ng
  grp_k = nonzeros(grp(k,:));
  if numel(grp_k) <= 1, continue; end  %single gen area -> no aggregation
  idx_m = find(ismember(ngen(:,1),grp_k(:)));
  agg_e = ngen(idx_m,2);
  for i = 1:1:numel(idx_m)
    im = idx_m(i);
    ib = find(nbus(:,1)==agg_e(i));
    vlt = nbus(ib,2); ang = nbus(ib,3)*pi/180;
    Pgen = nbus(ib,4)*ngen(im,22);
    Qgen = nbus(ib,5)*ngen(im,23);
    vph = vlt*exp(1i*ang);
    if ~trmbusagg
      xdp = basmva*ngen(im,7)/ngen(im,3);
    else
      xdp = 0.001*basmva*ngen(im,7)/ngen(im,3);  %aggregate only slightly "deeper" than the terminal bus
    end
    cur = conj((Pgen+1i*Qgen)./vph);
    vbk = vph + 1i*xdp*cur;
    Sgen = vbk*conj(cur);
    
    add_bus(1,1) = max(nbus(:,1))+1;
    add_bus(1,2) = abs(vbk); add_bus(1,3) = angle(vbk)*180/pi;
    add_bus(1,4) = real(Sgen);
    add_bus(1,5) = imag(Sgen);
    %add_bus(1,6) = nbus(ib,6); add_bus(1,7) = nbus(ib,7);  %dbg!
    %add_bus(1,8) = nbus(ib,8); add_bus(1,9) = nbus(ib,9);  %dbg!
    add_bus(1,6:9) = zeros(1,4);
    add_bus(1,10) = nbus(ib,10);
    oldgenbus = [oldgenbus; ib];
    nbus = [nbus;add_bus];
    
    nl = size(nlin,1);
    nlin(nl+1,1) = add_bus(1,1);%from bus
    nlin(nl+1,2) = nbus(ib,1);  %to bus
    nlin(nl+1,4) = xdp;         %xdp
    lin0(nl+1,1) = add_bus(1,1);%from bus
    lin0(nl+1,2) = nbus(ib,1);  %to bus
    lin0(nl+1,4) = xdp;         %xdp
    ngen(im,2)   = add_bus(1,1);%new genbus
    ngen(im,7)   = ngen(im,7) - xdp; %is changed later to Xdpeq anyway
    ngen(im,19)  = add_bus(1,1);%new genbus    
    ngen(im,22) = 1; ngen(im,23) = 1;
  end
end
nbus(oldgenbus,4) = 0; nbus(oldgenbus,5) = 0; 
%nbus(oldgenbus,6) = 0; nbus(oldgenbus,7) = 0;  %dbg!
%nbus(oldgenbus,8) = 0; nbus(oldgenbus,9) = 0;  %dbg! 
nbus(oldgenbus,10) = 3;  %PV/SW->PQ
%pstx = struct_pst(nbus, nlin, ngen, basmva, [], 60); 
%gx = pst2graph(pstx);
%part_viz(gx, ones(size(nbus,1),1), 'ext_graph');

if exist('pst_old\loadflow.m','file')==2
  [DUMMY,tst,~] = evalc('loadflow(nbus,nlin,tol,iter_max,vmin,vmax,acc,''n'',2);');   %(dbg!)
  assert( norm(vlt0(:,1)-tst(1:nb0,1),'fro')<1e-9 );
  assert( norm(vlt0(:,2)-tst(1:nb0,2),'fro')/nb0<1e-7 );
  assert( norm(vlt0(:,3)-tst(1:nb0,3),'fro')/nb0<1e-5 );
end
[DUMMY,nbus,lin_flw] = evalc('PST.LLoadflow(nbus,nlin,tol,iter_max,vmin,vmax,acc,''n'',2);');
bus0 = nbus;
assert( norm(vlt0(:,1)-nbus(1:nb0,1),'fro')<1e-9 );
assert( norm(vlt0(:,2)-nbus(1:nb0,2),'fro')/nb0<1e-7 );
assert( norm(vlt0(:,3)-nbus(1:nb0,3),'fro')/nb0<1e-5 );

% Retained buses go first (incl. internal bus indices). This should
% simplify the indexing of aggregated buses later on.
gen_a = nonzeros(grp); 
idx_m = ismember(ngen(:,1),gen_a);
agg_e = ngen(idx_m,2);
rtn_e = setdiff(nbus(:,1),agg_e); %retained-externalindexing
agg_i = find(ismember(nbus(:,1),agg_e));
rtn_i = find(ismember(nbus(:,1),rtn_e));
nbus_a = nbus(agg_i,:); %#ok<FNDSB>
nbus_r = nbus(rtn_i,:); %#ok<FNDSB>
[~,idx] = sort(nbus_a(:,1),'ascend');
nbus_a = nbus_a(idx,:);
[~,idx] = sort(nbus_r(:,1),'ascend');
nbus_r = nbus_r(idx,:);
nbus = [nbus_r;nbus_a];
nr   = size(nbus_r,1);

% Create ext2int and int2ext
bmx = max(nbus(:,1));
nb0 = size(nbus,1);
ext2int = zeros(bmx,1);
for i = 1:1:nb0
  ext2int(nbus(i,1)) = i;
end
assert(all(ext2int(nbus(:,1))==(1:1:nb0)'));
int2ext = zeros(1,nb0);
for i = 1:1:nb0
  [int2ext(i),~,x] = find(ext2int==i); 
end
assert(isequal(nbus(:,1)',int2ext(1:1:nb0))); %nbus(:,1) maps to 1:1:nb0

% Unreduced Ybus and initial reduced Ybus
Y = PST.YYbus(nbus,nlin);
Yjj = sum(Y,2);
Yii = sum(Y,1);
assert(norm(Yii(:)-Yjj(:),'fro')/nb0<1e-6);
rtn_i = find(ismember(nbus(:,1),rtn_e));
Yrr = Y(rtn_i,rtn_i);
Yagg = [[Yrr,zeros(nr,ng)];[zeros(ng,nr),zeros(ng,ng)]];
del_bus = [];
new_bus = [];
FLW_OUT = [];
FLW_BCK = [];
AGG_BUS = [];
SOUT = [];
IDXa = false(nb0,1);
Vout = zeros(nb0,1);
Iagg = zeros(ng,2);
for k = 1:1:ng
  grp_k = nonzeros(grp(k,:));  
  idx_m = find(ismember(ngen(:,1),grp_k(:)));
  agg_e = ngen(idx_m,2);
  agg_i = ext2int(agg_e);
  %Rearrange admittance matrix
  rtn_i = setdiff(1:1:nb0,agg_i);
  Yrc = Y(rtn_i,agg_i); Ycr = Y(agg_i,rtn_i); Ycc = Y(agg_i,agg_i);
  %Compute voltages, currents, and powers at "to be aggregated" buses
  Poutk = nbus(agg_i,4); Qoutk = nbus(agg_i,5);
  vlt_k = nbus(agg_i,2); ang_k = nbus(agg_i,3)*pi/180;
  vph_k = vlt_k.*exp(1i*ang_k); cur_k = conj((Poutk+1i*Qoutk)./vph_k);
  Saggk = vph_k.*conj(cur_k);
  assert(abs(sum(Poutk+1i*Qoutk-Saggk))<1e-6);  %reassure terminal bus usage
  %Compute aggregate voltage based on power balance at aggregate bus
  tot_cur = sum(cur_k); totsagg = sum(Saggk);
  vlt_agg = totsagg/conj(tot_cur);
  mag_agg = abs(vlt_agg); phi_agg = angle(vlt_agg);
  %Save initial electric output at the adjacent retained buses (for debugging)
  na = numel(agg_e); int_agg = zeros(na,1);
  adj_e = zeros(na,1); adj_i = zeros(na,1);
  flw_out = zeros(na,4); flw_bck = zeros(na,4);
  for i = 1:1:na
    row = find(lin_flw(:,2)==agg_e(i)); assert(nnz(row)==1);
    flw_out(i,:) = lin_flw(row,2:end);
    adj_e(i) = flw_out(i,2); adj_i(i) = ext2int(adj_e(i));
    int_agg(i,1) = ext2int(flw_out(i,1));
    flw_out(i,1) = agg_e(1);
    row = find(lin_flw(:,3)==agg_e(i)); assert(nnz(row)==1);
    flw_bck(i,:) = lin_flw(row,2:end); flw_bck(i,2) = agg_e(1);
  end
  FLW_OUT = [FLW_OUT;flw_out]; FLW_BCK = [FLW_BCK;flw_bck];
  Vout(adj_i) = nbus(adj_i,2).*exp(1i*nbus(adj_i,3)*pi/180);
  IDXa(adj_i) = true;
  SOUT = [SOUT; [agg_e(1)*ones(numel(adj_i),1),nbus(int_agg,4),nbus(int_agg,5)]];
  %Inertia-weighted aggregate voltage as in J.H. Chow (Iout won't match)
  %{
  if ~trmbusagg  %inertial aggregation
    H_i = ngen(idx_m,3).*ngen(idx_m,16)/basmva; H_s = sum(H_i);
    mag_agg = abs(vph_k)'*H_i/H_s;  %as in JChow&Co and Machowski
    phi_agg = angle(vph_k)'*H_i/H_s;  %as in JChow&Co and Machowski
  end
  %}
  %Phase shifter ratios TH for Yagg (in PST/MATPOWER modeling TPS = 1/TH)
  vlt_agg = mag_agg*exp(1i*phi_agg); th = vph_k/vlt_agg;
  th_abs = abs(vph_k)/mag_agg; th_ang = angle(vph_k)-phi_agg;
  %Aggregate admittances
  Yra = Yrc*th;
  Yar = conj(th.')*Ycr;
  Yaa = conj(th.')*Ycc*th;
  Yagg(nr+k,nr+k) = Yaa;
  Yagg(1:nr,nr+k) = Yra(1:nr);
  assert(abs(sum(Yra(nr+1:end)))<1e-6);
  Yagg(nr+k,1:nr) = Yar(1:nr);
  assert(abs(sum(Yar(nr+1:end)))<1e-6);
  %%%Add the aggregate generator bus
  AGG_BUS = [AGG_BUS; agg_e(1)];
  Iagg(k,:) = [agg_e(1),tot_cur];
  add_bus(1,1) = agg_e(1);
  add_bus(1,2) = mag_agg; add_bus(1,3) = phi_agg*180/pi;
  add_bus(1,4) = sum(nbus(agg_i,4)); add_bus(1,5) = sum(nbus(agg_i,5));
  add_bus(1,6) = sum(nbus(agg_i,6)); add_bus(1,7) = sum(nbus(agg_i,7));
  add_bus(1,8) = sum(nbus(agg_i,8));
  add_bus(1,9) = sum(nbus(agg_i,9));
  if any(nbus(agg_i,10)==1)
    bus_type = 1;
  else
    bus_type = 2;
  end
  add_bus(1,10) = bus_type;
  new_bus = [new_bus;add_bus];
  del_bus = [del_bus;agg_i];  %remove nbus(agg_i,:) with all Pgen/Qgen/Plod/Qlod
  %Correct R->A branches with ideal phase shifters
  na = numel(agg_e);
  nlnx = lin0;
  for i = 1:1:na
    f = agg_i(i);
    row = find(sum(nlin(:,1:2)==nbus(f,1),2)==1);
    rwx = find(sum(nlnx(:,1:2)==nbus(f,1),2)==1);
    assert(nnz(row)==1); assert(nnz(rwx)==1);
    to = adj_e(i); t = adj_i(i);
    nlin(row,1) = agg_e(1); nlin(row,2) = to;
    nlin(row,6) = 1/th_abs(i); nlin(row,7) = -th_ang(i)*180/pi; %TPS = 1/TH  
    nlnx(rwx,1) = agg_e(1); nlnx(rwx,2) = to;
    nlnx(rwx,6) = 1/th_abs(i); nlnx(rwx,7) = -th_ang(i)*180/pi; %TPS = 1/TH      
    assert(abs(-Y(f,t)-1/(nlin(row,3)+1i*nlin(row,4)))<1e-6);
    assert(abs(-Y(f,t)-1/(nlnx(rwx,3)+1i*nlnx(rwx,4)))<1e-6);
    assert(all(abs(nonzeros(Y(agg_i,t))-nonzeros(Yrc(t,:)))<1e-6));
    nbus(t,8) = nbus(t,8);
    nbus(t,9) = nbus(t,9);
  end
  %Aggregate generators
  xdpe = 0;
  idx0m = find(ismember(gen0(:,1),grp_k(:)));
  xdpa = basmva*gen0(idx0m,7)./gen0(idx0m,3);
  for i = 1:1:numel(xdpa)
    xdpe = xdpe + 1/xdpa(i);  % Ourari-Dessaint-Do: th_abs(i)*cos(th_ang(i))^2/xdpa(i)
  end
  xdpe = 1/xdpe;
  agg_mac = PST.eqgen(ngen,idx_m,basmva,agg_e(1),ngen(idx_m(1),1));
  agg_mac(1,7) = xdpe;  
  ngen(idx_m,:) = 0;  
  ngen = [ngen;agg_mac];  
  %Try to fit xd' from 0 to xdmax to preserve slow modes
  if modefit
    ngnx = gen0; nbsx = bus0;
    ngnx(idx_m,:) = []; nbsx(agg_i,:) = [];
    ngnx = [ngnx;agg_mac]; nbsx = [nbsx;add_bus];
    errs = Inf(1,75);
    xdpv = sort([xdpe,linspace(1e-6,max(0.05,30*xdpe),numel(errs)-1)]);
    for j = 1:1:numel(errs)
      ngnx(end,7) = xdpv(j);
      pstx = struct_pst(nbsx,nlnx,ngnx,basmva,[],pst.f0);
      [Vs, lamX, MKx, Kx, Mx] = PST.pstVslow(pstx,nm);
      mx = size(MKx,1);
      Amat = zeros(2*mx,2*mx);
      H = ngnx(:,16).*ngnx(:,3)/basmva;
      D = ngnx(:,17).*ngnx(:,3)/basmva;
      Amat(1:mx,mx+1:2*mx) = diag(w0*ones(mx,1));
      Amat(mx+1:2*mx,1:mx) = MKx/w0;
      Amat(mx+1:2*mx,mx+1:2*mx) = -0.5*diag(H)\diag(D);
      LAM = eig(Amat);
      [~, idx] = sort(abs(imag(LAM)),'ascend');
      LAM = LAM(idx(1:2:end-1));
      %assert(norm(abs(LAM(2:nm))-sqrt(abs(lamX(2:end))))/(nm-1)<1e-3)  %dbg            
      if any( real(LAM(real(LAM)>0))>0.005 )
        break;  
      end
      frqX = sort(abs(imag(LAM))/2/pi);
      errs(j) = sum(abs(frq0(2:nm)-frqX(2:nm))./frq0(2:nm));
      %plot(frq0(1:nm),'r.','MarkerSize',10); hold on;  %dbg
      %plot(frqX(1:nm),'b.','MarkerSize',10); hold off;  %dbg
      %xlim([1,nm]); set(gca,'xgrid','on'); set(gca,'xtick',1:nm);  %dbg      
    end
    [errmin,idxmin] = min(errs);
    ngen(end,7) = xdpv(idxmin);
  else
    if trmbusagg
      ngen(end,7) = xdpe;
    else
      ngen(end,7) = 1e-6;      
    end
  end
end
del_gen = sum(ngen,2)==0;
ngen(del_gen,:) = [];
nbus(del_bus,:) = [];
nbus = [nbus;new_bus];

% Check the aggregation result
%a. Check the agrregated admittance matrix
[DUMMY,nbus,lin_flw] = evalc('PST.LLoadflow(nbus,nlin,tol,iter_max,vmin,vmax,acc,''n'',2);');
Ya = PST.YYbus(nbus,nlin);
assert(norm(Ya(1:nr,1:nr)-Y(1:nr,1:nr),'fro')<1e-6);
assert(norm(Ya(1:nr,1:nr)-Yagg(1:nr,1:nr),'fro')<1e-6);
for k = 1:1:ng
  assert(norm(Ya(1:nr,nr+k)-Yagg(1:nr,nr+k),'fro')/nnz(Ya(1:nr,nr+k))<1e-6);
  assert(norm(Ya(nr+k,1:nr)-Yagg(nr+k,1:nr),'fro')/nnz(Ya(1:nr,nr+k))<1e-6);
end
assert(norm(Ya(nr+1:nr+ng,nr+1:nr+ng)-Yagg(nr+1:nr+ng,nr+1:nr+ng),'fro')<1e-6);
%b. Check loadflow voltages
vmag = nbus(:,2); vang = nbus(:,3)*pi/180;
vphr = vmag.*exp(1i*vang);
assert(max(abs(vphr(IDXa)-Vout(IDXa)))<1e-6);
assert( norm(vlt0(1:nr,1)-nbus(1:nr,1),'fro')<1e-9 );
assert( norm(vlt0(1:nr,2)-nbus(1:nr,2),'fro')/nb0<1e-7 );
assert( norm(vlt0(1:nr,3)-nbus(1:nr,3),'fro')/nb0<1e-5 );
%c. Check aggregated power flows
FLW_OUT = sortrows(FLW_OUT,[1,2]);
FLW_BCK = sortrows(FLW_BCK,[2,1]);
NEW_OUT = NaN(size(FLW_OUT));
NEW_BCK = NaN(size(FLW_BCK));
pos = 1;
for i = 1:1:numel(AGG_BUS)  
  row = lin_flw(:,2) == AGG_BUS(i);
  NEW_OUT(pos:pos+nnz(row)-1,:) = lin_flw(row,2:end);
  row = lin_flw(:,3) == AGG_BUS(i);
  NEW_BCK(pos:pos+nnz(row)-1,:) = lin_flw(row,2:end);    
  pos = pos+nnz(row);
end
NEW_OUT = sortrows(NEW_OUT,[1,2]);
NEW_BCK = sortrows(NEW_BCK,[2,1]);
assert(norm(FLW_OUT(:,1:2)-NEW_OUT(:,1:2),'fro')<1e-12);
assert(norm(FLW_BCK(:,1:2)-NEW_BCK(:,1:2),'fro')<1e-12);
assert(norm(FLW_OUT(:,3:4)-NEW_OUT(:,3:4),'fro')<1e-6);
assert(norm(FLW_BCK(:,3:4)-NEW_BCK(:,3:4),'fro')<1e-6);
%d. Check aggregated currents (Ie = sum(Ik))
Iagg = sortrows(Iagg,1);
Iout = zeros(ng,2);
for k = 1:1:ng
  row = find(nbus(:,1)==AGG_BUS(k));
  Pout = nbus(row,4);
  Qout = nbus(row,5);
  Vout = nbus(row,2)*exp(1i*nbus(row,3)*pi/180);
  Iout(k,:) = [nbus(row,1), conj((Pout + 1i*Qout)/Vout)];  
end
Iout = sortrows(Iout,1);
assert(norm(Iout(:,1)-Iagg(:,1),'fro')/ng<1e-12);
assert(norm(Iout(:,2)-Iagg(:,2),'fro')/ng<1e-6);
%[mat2cell(lin_flw(:,2:3),ones(size(lin_flw,1),1),ones(1,2)),mat2cell(lin_flw(:,4:5),ones(size(lin_flw,1),1),ones(1,2))] -> shows linflows
red_mod.lin = nlin;
red_mod.bus  = nbus;
red_mod.gen = ngen;
red_mod.basmva = pst.basmva; 
red_mod.f0 = pst.f0; 
end



