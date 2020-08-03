classdef MatpowerIn < BaseIn 
  %MatpowerIn: Class for input data required by network partitioning
  %  algorithms. This data representation is used by Matpower, and can be
  %  converted to Power System Toolbox.
  
  % MatpowerIn containing all power flow specifications for each case was
  % chosen instead of an array of scalar MatpowerIn object for each lf as:
  % 1) It stores case resutfolder (which is always one per case). So an 
  %    object is gone -> the assotiated folder will be removed as well 
  %    (if empty)
  % 2) dev_pq + seed combination of properties with the fixed base seed and 
  %    the increments to it equal to the number of the current test case's 
  %    load flow is good for the reproductability of the results.
  
  
  properties (SetAccess = private)    
    slack  % numbers of generator buses that should participate in power balancing (whole number vector)
    seed  % fixed seed for random number generator as in rng(seed)    
  end  
    
  % These properties can be modified for every class object (sensible setters are beneficial)
  properties
    f0  % base network frequency
    curr  % number of the currently considered power flow case (whole number scalar)
    mean_pq  % mean P-Q load scaling of all load buses, one per power flow case (whole number vector)    
    dev_pq  % allowed deviation range of load scaling at each particular bus, one per power flow case (whole number vector)      
  end
  
  properties (Constant, Access = private)
    dflt_f0 = 50;  % default value of test case frequency
  end

  
  methods
    [bus, lin, gen, basmva] = mp2bgl( obj, varargin )
    
    
    function obj = MatpowerIn(varargin)
      % MatpowerIn initializer. It also works with no nargin == 0, but the 
      % maximal nargin is 20 (10 parameter-value pairs) 

      assert(length(varargin) <= 32, [ mfilename(), ':TooManyInputs'],...
        '[%s] A %s object has at most 16 optional parameter-value pairs',...
        mfilename, mfilename);      
      
      % Create input parser and search for the inputs of the superclass
      p = inputParser;
      p = Utils.inputParserSetup(p);
      p.addParameter('caseid', []);
      p.addParameter('n_pf', []);
      p.addParameter('n_coh', []);
      p.parse(varargin{:});
      baseInp = p.Results; 
      mpwInp = p.Unmatched;
      fld = fieldnames(baseInp);
      for i = 1:1:length(fld)
        tmp = baseInp.(fld{i});
        if isempty(tmp)
          baseInp = rmfield(baseInp, fld{i});
        end
      end
            
      % Call the superclass constructor
      obj = obj@BaseIn(baseInp);
      
      % Default values for optional inputs (MatpowerIn-properties only!)    
      p = inputParser;
      p = Utils.inputParserSetup(p);
      p.addParameter('f0', MatpowerIn.dflt_f0);
      p.addParameter('mean_pq', MatpowerIn.dflt_mean_pq(obj.n_pf));
      p.addParameter('dev_pq', MatpowerIn.dflt_dev_pq(obj.n_pf));
      p.addParameter('slack', []);
      p.addParameter('seed', 1, @(x) isscalar(x) && x>0 && Utils.isint(x));
      p.parse(mpwInp);
      mpwInp = p.Results;
      trash = p.Unmatched;
      assert(isempty(fieldnames(trash)), [ mfilename, ':WrongKeyValueInput'],...
        ['[%s] Some unexpected key value pairs are detected. Please check your ',...
          'spelling. The legitimate keys for the constructor are: caseid, n_pf, n_coh, ',...
          'f0, mean_pq, dev_pq, slack, seed. ',...
          'The following properties are generated automatically in the ',...
          'constructor: id, excelfile, rng, curr. The following properties ',...
          'are set through the %s.initScalPQ method based on the object''s ',...
          'caseid, curr, mean_pq, dev_pq, slack, n_coh and rng properties: ',...
          'bus, lin, gen, coh, basmva.'], mfilename, mfilename);
      obj.f0 = mpwInp.f0;
      obj.slack = mpwInp.slack;
      obj.mean_pq = Utils.chkVectorScalar(obj, mpwInp, 'mean_pq', 'n_pf',...
        @MatpowerIn.dflt_mean_pq, mfilename);  % check size of mean_pq
      obj.dev_pq = Utils.chkVectorScalar(obj, mpwInp, 'dev_pq', 'n_pf',...
        @MatpowerIn.dflt_dev_pq, mfilename);  % check size of dev_pq
      
      % Check the relative dimensions of some properties
      assert(all(obj.mean_pq > obj.dev_pq),...
        [ mfilename, ':IncompatibleProperties'],...
        ['[%s] Elements of %s.dev_pq should be less than the corresponding ',... 
        'elements of %s.mean_pq'], mfilename, mfilename, mfilename);     
      
      % Save randomly initialized rng and initialize obj.curr to 1
      obj.seed = rng('shuffle', 'twister');
      obj.seed.Seed = mpwInp.seed;
      obj.curr = 1;
    end
    
    function [scal, slack] = initScalPQ(obj, rng_mode)
      % initScalPQ Generates actual load scalings for each bus for
      %   ith power flow of the system described by caseid, based on 
      %   mean_pq(i) and dev_pq(i), where i = obj.curr
      
      assert(obj.curr<=obj.n_pf, [ mfilename,...
        ':IncompatibleProperties'], ['[%s] %s.curr should be less than ',... 
        'or equal to %s.n_pf. It seems that the %s.initScalPQ method was run ',...
        'too many times for the same object.'], mfilename, mfilename,...
        mfilename, mfilename);      
      
      i = obj.curr;
      caseid = obj.caseid;  % from here we get the actual bus list
      m_pq = obj.mean_pq(i);
      d_pq = obj.dev_pq(i);
      s = rng;  % save rng
      if rng_mode
        rng(obj.seed.Seed + i, obj.seed.Type); % setup rng from the stored object seed for the current i
      else
        rng('default');
      end
      
      netw = loadcase(caseid);  % load network from MATPOWER
      nbus = size(netw.bus, 1);
      ngen = size(netw.gen, 1);
      
      nslack = floor(0.5*ngen);  % default number of power balancing generators
      slack = logical([ones(nslack,1); zeros(ngen-nslack,1)]);
      scal = (m_pq-d_pq) + 2*d_pq*rand(nbus,1);
      if isempty(obj.slack)
        if d_pq == 0
          rng(1, 'twister');  % determitistic slack for ANY object if fixed load scaling
        end
        slack = find(slack(randperm(ngen)));  % if no slack function input, slack generators are randomly (but predictably) chosen
      else
        slack = obj.slack(:);
      end
      rng(s);  % restore rng to its initial setting
    end  
    
    function set.f0(obj, f0)
      assert(isscalar(f0) && (f0==50 || f0==60),...
        [ mfilename, ':InvalidProperty'],...
        ['[%s] %s.f0 should be either 50 or 60 Hz (i.e. the only two ',...
        'meaningful values)'], mfilename, mfilename);
      obj.f0 = f0;      
    end
    
    function set.mean_pq(obj, mean_pq)
      assert(isvector(mean_pq) && Utils.isdbl(mean_pq) && all(mean_pq>0),...
        [ mfilename, ':InvalidProperty'],...
        ['[%s] %s.mean_pq should be a vector of positive doubles (one ',... 
        'number for each caseid bus) or a nonnegative double scalar (one ',...
        'number for ALL caseid buses)'], mfilename, mfilename);
      obj.mean_pq = mean_pq(:);      
    end
    
    function set.dev_pq(obj, dev_pq)
      assert(isvector(dev_pq) && Utils.isdbl(dev_pq) && all(dev_pq>=0),...
        [ mfilename, ':InvalidProperty'],...
        ['[%s] %s.dev_pq should be a vector of nonnegative doubles (one ',...
        'number for each caseid bus) or a nonnegative double scalar (one ',...
        'number for ALL caseid buses)'], mfilename, mfilename);
      obj.dev_pq = dev_pq(:);        
    end
    
    function set.slack(obj, slack)
      assert(isempty(slack) || (isvector(slack) && Utils.isint(slack) && all(slack>0)),...
        [ mfilename, ':InvalidProperty'],...
        '[%s] %s.slack should be a logical vector or empty', mfilename, mfilename);
      obj.slack = slack(:);
    end
    
    function set.curr(obj, curr)
      assert(isscalar(curr) && Utils.isint(curr) && curr>0,...
        [ mfilename, ':InvalidProperty'],...
        '[%s] %s.curr should be a positive whole number scalar', mfilename, mfilename);
      obj.curr = curr;
    end
            
  end
  
  
  methods (Static, Access = private)
    function mean_pq = dflt_mean_pq(n_pf, new_dflt)
      % dflt_mean_pq returns default mean_pq as unity for each power flow
      % case, unless a 2nd input is provided which would define the value 
      % for each power flow case 
      
      switch nargin
        case 1
          mean_pq = ones(n_pf, 1);
        case 2
          mean_pq = new_dflt*ones(n_pf, 1);
        otherwise
          error([ mfilename, ':WrongPrivateInput'],...
            ['[%s] Wrong input to a private method: this shouldn''t ',...
            'ever happen!'], mfilename);
      end
    end
    
    function dev_pq = dflt_dev_pq(n_pf, new_dflt)
      % dflt_dev_pq returns default dev_pq as zero for each power flow 
      % case, unless a 2nd input is provided which would define the value 
      % for each power flow case
      
      switch nargin
        case 1
          dev_pq = zeros(n_pf, 1);
        case 2
          dev_pq = new_dflt*ones(n_pf, 1);
        otherwise
          error([ mfilename, ':WrongPrivateInput'],...
            ['[%s] Wrong input to a private method: this shouldn''t ',...
            'ever happen!'], mfilename);          
      end
    end
  end
end


