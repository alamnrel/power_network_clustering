classdef Utils
  %Utils Static class with utility methods which are used by other classes
  
  methods (Static)
   
    
    function isWhole = isint(x)
      % isint Checks if array x all consists of real whole numbers  
      %   (integers or doubles with zero decimal part)
      
      isWhole = isnumeric(x) && all(~isempty(x(:))) && ...
         all(~isa(x(:), 'char') & ~isinf(x(:)) & ~isnan(x(:)) & ...
         imag(x(:))==0 & mod(x(:), 1)==0);
    end 
    
    function isDouble = isdbl(x)
      % isdbl Checks if array x all consists of real double numbers
      
      isDouble = isnumeric(x) && all(~isempty(x(:))) && ...
        all(~isinf(x(:)) & ~isnan(x(:)) & imag(x(:))==0 & isa(x(:), 'double'));
    end
    
    function str = rnd_str(stLength)
      % rnd_str - Generate a random string
      %
      % When a name is needed but none is provided, never fear
      % rnd_str is here! It returns alphanumeric random strings.
      %
      % Author: StackOverflow & Ilya Tyuryukanov
      % Date of first version: ?? 2015
      % Last revision: 24 December 2015
      
      % Check inputs
      if (nargin < 1)
        stLength = 10; % if no stLength is specified, use 10
      else
        assert(Utils.isint(stLength) && isscalar(stLength) && stLength>0,...
          [ mfilename(), ':InvalidProperty'],...
          '[%s] The input of %s.rnd_str must be a positive integer',...
          mfilename(), mfilename());
      end
      
      % Generate random string
      symbols = ['a':'z' 'A':'Z' '0':'9'];
      nums = randi(numel(symbols), [1 stLength]);
      str = symbols (nums);
    end
    
    function V = normalize_rows(V)
    V = bsxfun(@rdivide, V, sqrt(sum(V.^2,2)));
    end      
    
    function p = inputParserSetup(p)
      %This function takes a just created inputparser object and sets its
      % properties in a proper way
      
      if isprop(p, 'PartialMatching'), p.PartialMatching = false; end
      p.StructExpand = true;
      p.KeepUnmatched = true;
      p.CaseSensitive = false;      
    end    
    
    function prop = chkVectorScalar(obj, data, prop1, prop2, fcn, caller)
      % chkVectorScalar implements one common set of rules for checking  
      %   the length of some vector properties against the values of some 
      %   scalar properties
      
      if numel(data.(prop1)) == obj.(prop2)  % check size of dev_pq
        prop = data.(prop1);
      elseif numel(data.(prop1)) == 1
        prop = fcn(obj.(prop2), data.(prop1));
      else
        error([ caller, ':InvalidPropertySize'],...
        ['[%s] %s.%s should have its length equal to %s.%s. May ',... 
        'be it was provided without %s.%s or with neglect of the %s.%s value'],...
        mfilename, caller, prop1, caller, prop2, caller, prop2, caller, prop2);
      end      
    end   
            
    function [large_idx, cardinality] = large_comp(C, cutoff)
      if nargin < 2
        cutoff = 0;
      end
      comp_idx = unique(C, 'sorted');  % component indices
      cardinality = histcounts(C, 'BinMethod', 'integers');  % # of nodes in components with indices 'comp_idx'
      large = cardinality > cutoff;  % indices (in 'cardinality') of large components
      cardinality = cardinality(large);
      [cardinality, idx] = sort(cardinality, 'descend');
      large_idx = comp_idx(large);
      large_idx = large_idx(idx);  % largest first
    end    
    
  end
  
end
