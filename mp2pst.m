function pst = mp2pst( bus, varargin) 
%
% Syntax: pst = mp2pst( bus, varargin )
% 
% Purpose: 
% convert the data in MATPOWER format to PST (Power System Toolbox by  
% J. Chow, G. Rogers & L. Vanfretti) format; transform a MATPOWER study 
% case to a PST study case with some dynamic data added.
%
% Input: 
% obj: a MatpowerIn object or a MATPOWER case structure.
%
% Output: 
% obj: the updated input MatpowerIn object. It is updated by the bus, line
%   and gen data in the PST format for the ith power flow case. This added 
%   data also includes the solution of each load flow and baseMVA. The load
%   flow solution can be either for AC or DC load flow.
% 
% Remarks:
% By MATPOWER implementation, there can be only one slack generation bus 
% in a connected network. Therefore,  we neglect the "distributed" slack 
% output of initScalPQ.
% 
% Author: Ilya Tyuryukanov
% Date: 14 July 2015
% Revision: 15 January 2016
% Last revision: 16 December 2017 (not a MatpowerIn method anymore)
% 
% Correspondence between MATPOWER and PST formats:
%---------------------------------------------------------------------------------------------------|
% bus:                                                                                              |
%           PST              |          MATPOWER                                                    |
% col1 number                | 1  bus number (positive integer)                                     |
% col2 voltage mag(pu)       | 8  Vm, voltage magnitude (p.u.)                                      |
% col3 voltage angle(deg)    | 9  Va, voltage angle (degrees)                                       |
% col4 p_gen(pu)             | 2G Pg, real power output (MW)                                        |
% col5 q_gen(pu),            | 3G Qg, reactive power output (MVAr)                                  |
% col6 p_load(pu)            | 3  Pd, real power demand (MW)                                        |
% col7 q_load(pu)            | 4  Qd, reactive power demand (MVAr)                                  |
% col8 G shunt(pu)           | 5  Gs, shunt conductance (MW demanded at V = 1.0 p.u.)               |
% col9 B shunt(pu)           | 6  Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)             |
% col10 bus_type             | 2  bus type                                                          |
%       1: swing bus         |        1: PQ bus                                                     |
%       2: PV bus            |        2: PV bus                                                     |
%       3: PQ bus            |        3: reference bus                                              |
% col11 q_gen_max(pu)        | 4G Qmax, maximum reactive power output (MVAr)                        |
% col12 q_gen_min(pu)        | 5G Qmin, minimum reactive power output (MVAr)                        |
% col13 v_rated (kV)         | 10 baseKV, base voltage (kV)                                         |
% col14 v_max  pu            | 12 maxVm, maximum voltage magnitude (p.u.)                           |
% col15 v_min  pu            | 13 minVm, minimum voltage magnitude (p.u.)                           |
%---------------------------------------------------------------------------------------------------|
% line:                                                                                             |
%           PST              |          MATPOWER                                                    |
% 1  from bus                | 1  f, from bus number                                                |
% 2  to bus                  | 2  t, to bus number                                                  |
% 3  resistance(pu)          | 3  r, resistance (p.u.)                                              |
% 4  reactance(pu)           | 4  x, reactance (p.u.)                                               |
% 5  line charging(pu)       | 5  b, total line charging susceptance (p.u.)                         |
% 6  tap ratio               | 9  ratio, transformer off nominal turns ratio ( = 0 fmpwor lines )   |
% 7  tap phase               | 10 angle, transformer phase shift angle (degrees), positive -> delay |
% 8  tapmax                  | ...                                                                  |
% 9  tapmin                  | ...                                                                  |
% 10 tapstep                 | ...                                                                  |
% 11 MUST-LINK (NEW!)        | ...                                                                  |
% 12 linflow (NEW!)          | ...                                                                  |
%---------------------------------------------------------------------------------------------------|
% gen:                                                                                              |
%           PST              |          MATPOWER                                                    |
% 1  machine number          |                                                                      |
% 2  bus number              | 1  bus number                                                        |
% 3  base mva (machine)      | 7  mBase, total MVA base of this machine, defaults to baseMVA        |
%                            |           (i.e. to the network base power)                           |
% 4  x_l(pu)                 |                                                                      |
% 5  r_a(pu)                 |                                                                      |
% 6  x_d(pu)                 |                                                                      |
% 7  x'_d(pu)                |                                                                      |
% 8  x"_d(pu)                |                                                                      |
% 9  T'_do(sec)              |                                                                      |
% 10 T"_do(sec)              |                                                                      |
% 11 x_q(pu)                 |                                                                      |
% 12 x'_q(pu)                |                                                                      |
% 13 x"_q(pu)                |                                                                      |
% 14 T'_qo(sec)              |                                                                      |
% 15 T"_qo(sec)              |                                                                      |
% 16 H(sec)                  |                                                                      |
% 17 d_o(pu)                 |                                                                      |
% 18 d_1(pu)                 |                                                                      |
% 19 bus number              | 1  bus number                                                        |
% 20 saturat. factor S(1.0)  |                                                                      |
% 21 saturat. factor S(1.2)  |                                                                      |
% 22 active power fraction   |                                                                      |
% 23 reactive power fraction |                                                                      |
%---------------------------------------------------------------------------------------------------|

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addOptional('branch', [], @(x)isnumeric(x)&&ismatrix(x) );
p.addOptional('gen', [], @(x)isnumeric(x)&&ismatrix(x) );
p.addOptional('basmva', [], @(x)isnumeric(x)&&ismatrix(x) );
p.addOptional('f0', [], @(x)isnumeric(x)&&isscalar(x)&&x>0 );
p.addParameter('acdc', 'ac', @(x)ischar(x)&&length(x)==2 );
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
  'spelling. The allowed keys are: acdc.'], mfilename);
branch = varinp.branch;
gen = varinp.gen;
basmva = varinp.basmva;
f0 = varinp.f0;
acdc = varinp.acdc;

% Analyse inputs
if isa(bus, 'struct')
  mpw = bus;  
else
  mpw = struct_mp(bus, branch, gen, basmva);
  mpw.f0 = f0;
end
if isempty(mpw.bus)
  pst = struct_pst([], [], [], 0, [], 0);
  return;
end

% Solve power flow for the whole network
mpOpt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1);
if strcmp(acdc,'ac')
  mpw = runpf(mpw, mpOpt);
else
  mpw = rundcpf(mpw, mpOpt);
end
if ~mpw.success
  warning([mfilename,':PowerFlowNotConverged'],...
    ['[%s] For the network %s and the config number %d, no power flow ',...
    'solution could be found'], mfilename, obj.caseid, obj.curr-1);
  mpw = struct_mp([], [], [], 0, 0);
  return;  %isempty(mpw.bus)==true <=> no loadflow for ith load scaling
end

% Delete out-of-service (also adds order field to mpw...)
mpw = mp_del_oos(mpw);

% Set up power flow case
bus = mpw.bus;
branch = mpw.branch;
gen = mpw.gen;
basmva = mpw.baseMVA;
f0 = mpw.f0;
nbus = numel(bus(:,1));
nlin = numel(branch(:,1));

% Find and convert dynamic generator data
idxgen = arrayfun(@(x)find(bus(:,1)==x,1), gen(:,1));  % indices of each gen (1st->2nd->3rd->...->end) in bus[]
genreqs = [gen(:,1), gen(:,9)/0.85, bus(idxgen,10)];  % Note that gen(:,9) is the vector of generators' Pmax, the divider is the "default nominal cosphi"
         % busID   , Sg_base (MVA), Vn_busID          
genPST = findGen(genreqs, basmva);

% Bus data conversion
busType = bus(:,2);
tmp = busType;
busType(tmp == 1) = 3;  % PQ bus
busType(tmp == 2) = 2;  % PV bus
busType(tmp == 3) = 1;  % slack bus
Pg = zeros(nbus,1);
Qg = zeros(nbus,1);
Qgmax = zeros(nbus,1);
Qgmin = zeros(nbus,1);
pgbus = accumarray(idxgen, gen(:,2));  % cumulative Pgen at bus
[row, ~, pgbus] = find(pgbus);
Pg(row) = pgbus;
qgbus = accumarray(idxgen, gen(:,3));  % cumulative Qgen at bus
[row, ~, qgbus] = find(qgbus);
Qg(row) = qgbus;
qgbushi = accumarray(idxgen, gen(:,4));  % cumulative Plim at bus
[row, ~, qgbushi] = find(qgbushi);
Qgmax(row) = qgbushi;
qgbuslo = accumarray(idxgen, gen(:,5));  % cumulative Qlim at bus
[row, ~, qgbuslo] = find(qgbuslo);
Qgmin(row) = qgbuslo;

bus = [bus(:,1), bus(:,8), bus(:,9), Pg, Qg, bus(:,3), bus(:,4),...
  bus(:,5), bus(:,6), busType, Qgmax, Qgmin, bus(:,10), bus(:,12), bus(:,13)];
bus(:,4:5) = bus(:,4:5)/basmva;  % MW/MVar to p.u. {Pg and Qg}
bus(:,6:7) = bus(:,6:7)/basmva;  % MW/MVar to p.u. {Pload and Qload}
bus(:,8:9) = bus(:,8:9)/basmva;  % MW/MVar to p.u. {Gshunt and Bshunt}
bus(:,11:12) = bus(:,11:12)/basmva;  % MW/MVar to p.u. {Qg,min and Qg,max}

% Line data conversion
lin = [branch(:,1:5), branch(:,9:10), zeros(nlin,3)];
linflow = (abs(branch(:,14)) + abs(branch(:,16)))/2;   % MATPOWER calculates it in MW
linflow = linflow/basmva;   % convert MW to p.u.
trafos = arrayfun(@(x,y)(x~=0 || y~=0), branch(:,9), branch(:,10));
lin = [lin, trafos, linflow];

% Generator data insertion
pst = struct_pst(bus, lin, genPST, basmva, [], f0);
end


function genPST = findGen(sample, baseMVA)
% 
% Syntax: genPST = findGen(tmpl, baseMVA)
%
% Purpose: 
% add dynamic generator data in PST (Power System Toolbox by J. Chow, 
% G. Rogers & L. Vanfretti) format to the loadflow model (e.g. taken from
% MATPOWER, but not only)
%
% Input: 
% tmpl - load flow generator data to find the nearest dynamic model from
%        the database. It consists of 3 columns (all double): 
%       [busID]--[nominal MVA estimate]--[nom. bus kV in lf data]
% baseMVA - Matpower network base power (double)
%
% Output: 
% genPST - nGen x 21 matrix of synchronous generator parameters for each
%          SG in the network model, in PST format (double)          
%
% Remarks:
% For coherency studies, classical generator model is assumed. The implicit
% assumption for this is xd' = xq'.
% For simplicity, we assume that all inserted dynamic generator models are 
% connected to the grid through an IDEAL unit trafo with low voltage equal  
% to the generator nominal voltage and high voltage equal to the nominal 
% voltage of the bus to which the generator is connected in the load flow 
% model specification (e.g. of IEEE300). Transformer short circuit ratio  
% is assumed 0% on the equipment base (generally uk is 6.5 - 14%, for 
% medium-size transformers it's usually 10-11%).
%
% Author: Ilya Tyuryukanov
% Date: 14 July 2015
% Last revision: 23 July 2015
  
% Load data table 'Generators' along with explanation 'GenDataLegend'
load('Generators_All.mat');  % run ExcelMakeGenData to create the corresponding mat file       
genData = table2array(Generators); 

% Introduce some data prefilters based on saliency of generator models
% Any other prefilters can be conceived to exclude some group of generators
symmSS = abs((genData(9,:) - genData(4,:))./genData(9,:)) < 0.125;  % index vector for steady state nonsaliency
symmTr = abs((genData(10,:) - genData(5,:))./genData(10,:)) < 0.5;  % index vector for transient nonsaliency (very rough)
symmSsTr = symmSS & symmTr;   % index vector for steady state AND transient nonsaliency...
genData = genData(:,symmTr);  % to consider all the generators, comment out this line

genData = genData';  % in PST, parameters are in columns instead of rows

% Aggregate generators on the same bus
sample = sortrows(sample,1,'ascend');
buscnt = tabulate(sample(:,1));
ng = numel(buscnt(:,2));
sampl0 = sample;
sample = [];
for i = 1:1:ng
  if buscnt(i,2)==0
    continue
  else
    bus = buscnt(i,1);
    row = find(sampl0(:,1)==bus);
    sample = [sample;[sampl0(row(1),1),sum(sampl0(row,2)),sampl0(row(1),3)]];
  end  
end

ngen = size(sample,1);
Sn_estimate = sample(:,2);
Sn_database = genData(:,1);
C = ipdm(Sn_estimate, Sn_database);  % similarity matrix from each estimate of Sn... 
                                     % to each Sn available in the "database" genData
[~, I] = min(C,[],2);  % compute smallest distance in each row 
genPST = genData(I,:);  % choose the nearest generators from the "database" as network generators

genPST = insertrows(genPST', sample(:,1)', 16)';  % insert bus number in accordance with PST format
genPST = [(1:1:ngen)', sample(:,1), genPST];  % insert generator number (quite arbitrary) and 
                                            % bus number in accordance with PST format
% Prepare the final generator data matrix for output:
genKV = genPST(:,20); %just fyi
busKV = sample(:,3);  %just fyi
genPST(:,20:end) = [];
% 1. Average the values of reactances over d and q axes (we assume to neglect transient saliency):
genPST(:,7) = (genPST(:,7) + genPST(:,12))/2;  % xd' = {average of xd' and xq'}
genPST(:,12) = genPST(:,7);  % xq' = {average of xd' and xq'}
% 2. Transform generator impedances from equipment base to network base
% and set base power to network base power for every generator
genPST(:,4) = genPST(:,4).*baseMVA./genPST(:,3);
genPST(:,5) = genPST(:,5).*baseMVA./genPST(:,3);
genPST(:,6) = genPST(:,6).*baseMVA./genPST(:,3);
genPST(:,7) = genPST(:,7).*baseMVA./genPST(:,3);
genPST(:,8) = genPST(:,8).*baseMVA./genPST(:,3);
genPST(:,11) = genPST(:,11).*baseMVA./genPST(:,3);
genPST(:,12) = genPST(:,12).*baseMVA./genPST(:,3);
genPST(:,13) = genPST(:,13).*baseMVA./genPST(:,3);
genPST(:,16) = genPST(:,16).*genPST(:,3)./baseMVA;
genPST(:,3) = ones(length(genPST(:,3)),1)*baseMVA;    
genPST(:,19) = genPST(:,1);  % an option; can be commented out
end