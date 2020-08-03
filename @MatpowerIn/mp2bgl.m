function mpw = mp2bgl( obj, varargin )

% Initialization
data = loadcase(obj.caseid);
nbus = length(data.bus(:,1));

% Input processing
p = inputParser;
p = Utils.inputParserSetup(p);
p.addParameter('lsc', 'var', @(x) ischar(x) || (Utils.isdbl(x) &&...
  size(x,1)==nbus && size(x,2)==2 && all(x(:)>=0) && all(x(:)<=1)) );
p.addParameter('acdc', 'ac', @(x)ischar(x)&&length(x)==2 );
p.parse(varargin{:});
varinp = p.Results;
trash = p.Unmatched;
assert(isempty(fieldnames(trash)), [mfilename,':WrongKeyValueInput'],...
  ['[%s] Some unexpected key value pairs are detected. Please check your ',...
  'spelling. The allowed keys are: lsc, acdc.'], mfilename);
lsc = varinp.lsc;
acdc = varinp.acdc;

% Data preprocessing
[data.bus(:,1),I] = sort(data.bus(:,1),'ascend');  % sort with respect to busIDs to ensure that... 
data.bus(:,2:end) = data.bus(I,2:end);  
[data.gen(:,1),I] = sort(data.gen(:,1),'ascend');  % busIDs in bus[] and gen[] are in the same order
data.gen(:,2:end) = data.gen(I,2:end);
if isfield(data,'gencost')
  data.gencost = data.gencost(I,:);
end

% Scale P-Q loads at the load buses
if ischar(lsc)
  if strcmp(lsc,'var')
    scalpq = obj.initScalPQ(true);  % retrieve bus load scalings for the current power flow case    
  else
    scalpq = obj.initScalPQ(false);
  end
  scalp = scalpq; 
  scalq = scalpq; 
else
  scalp = lsc(:,1);
  scalq = lsc(:,2);  
end
data.bus(:,3) = data.bus(:,3).*scalp;  % the only thing which can be adjusted...
data.bus(:,4) = data.bus(:,4).*scalq;  % because no "distributed slack" in MATPOWER
mpOpt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1);
if strcmp(acdc,'ac')
  data = runpf(data, mpOpt);
else
  data = rundcpf(data, mpOpt);
end
if ~data.success
  warning([mfilename,':PowerFlowNotConverged'],...
    ['[%s] For the network %s and the config number %d, no power flow ',...
    'solution could be found'], mfilename, obj.caseid, obj.curr-1);
  mpw = struct_mp([], [], [], 0, 0);
  return;  % isempty(mpw.bus)==true <=> no loadflow for ith load scaling
end

if isfield(data,'gencost')
  mpw = struct_mp(data.bus, data.branch, data.gen, data.baseMVA, data.gencost);
else
  mpw = struct_mp(data.bus, data.branch, data.gen, data.baseMVA);
end
mpw.f0 = obj.f0;
end



%data.branch(:,5) = 0;   % remove any reactive shunts
%data.bus(:,[5,6]) = 0;   % remove any reactive shunts
%data.branch(:,[9,10]) = 0;   % remove any (phase shifting) trafos