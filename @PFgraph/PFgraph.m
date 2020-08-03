classdef PFgraph 
  % PFgraph represents a (power flow) graph by adjacency and incidence 
  % matrices, it also contains infos about coherency grouping of generator 
  % nodes, what nodes are generators and what branches cannot be cut  
  % 
  % PFgraph.coh: max(g.coh(2,:)) should always represent the requested number of 
  %   clusters/islands, so a care should be taken when assigning the integer group
  %   labels in g.coh(2,:) for every generator
  
  properties (SetAccess = private)
    adj  % graph adjacency matrix (symmetrical real matrix with nonnegative entries)
    inc  % graph incidence matrix (matrix composed of ones and zeros)    
    scal  % scaling factor for adj to return it back to the user. The real adjacency matrix can be sometimes premultiplied with 1/scal for a better numeric range of values in computations
    Loff  % unnormalized Laplacian
    Lrwk  % "random-walk" normalized Laplacian
    Lsym  % "symmetric" normalized Laplacian   
    merge_map  % mapping vector between the original and the contracted (merged) graph vertices (relevant in case the graph is a result of a node contraction); it is set by default to a 1:1 mapping    
    basmva  % base power (network base) 
  end
  
  properties (SetAccess = public)
    bus  % mapping from bus indices (in adj) to bus labels (vector of initial bus numbers, or cell array of text labels) 
    gen  % for every bus shows if it is generator (logical vector) 
    coh  % for certain nodes shows to which "coherent group" they belong to (two-rows matrix with nodes in the 1st row) 
    ml  % for each branch shows if it is a must-link branch or not (logical vector) 
    vw  % generic vertex weigths (combined with adj, may participate in a normalized spectral clustering problem) 
    powgen  % generated power at each node as a [P, Q]-vector 
    powlod  % load power at each node as a [P, Q]-vector 
    vmag  % power flow solution for voltage of each node 
    vnom  % base voltage of each node 
    vang  %  power flow solution for voltage phase angle 
  end
  
  methods    
    [ok, bran_viols, coh_viols] = chkBranCoh(g)
    trees = coh2sptree(g)  	
    cutind = cutset(g, busind)
	busind = fiedler_bisect(g)
    [T, contig] = kahip_improve_v1(g, I, ref)
    [T_out, minor_cc] = contig(g, T, varargin)
    [idx_out, val] = edges2adj( g, linidx, NDEBUG )
    g_redu0 = reduce_pml( g )
    [eXp, pcut, ixs, adj_sep, num_viols, viols, tot_outl_gen] = ici_info( g, cutind )
    g_cntrct = merge_nodes(g, ixs)
    [cutind, busind, valid, eXp, pcut, n_outl_gen, dif] = metis_info(g, map)
    graph_py = pfgraph2igraph( g, clust ) 
    [Y, v] = spClust(g, nm, nc) 
    pfgraph2hmetis(g, vw)
    [T, UBfactor, n_stat] = ucHMETIS(g, method, fix, imb)  
    [T,           n_stat] = ucGRACLUS(g, varargin)  
    [T] = ucSpEmbGrPart(g, varargin) 
    [T] = coGrReduClust(g, trees, varargin)
    [ixs_out, ixs_ext] = outliers_kwaysp(g, varargin)
    [adj_sep, C, ixs, eXp, pcut] = cc_info(g, cutind)
    [cutind, busind, T] = final_cutset(g, T, merg)
     trees = hclust_constrained_path(g)
    [T] = uckwaycut(g, core_nodes, method, card_min)
    
    function obj = PFgraph(varargin)
      % PFgraph initializer. It also works with no nargin == 0, but the 
      % maximal nargin is 8 (4 parameter-value pairs)         
                      
      assert(length(varargin) <= 38, [ mfilename(), ':TooManyInputs'],...
        '[%s] A %s object has at most 19 parameter-value pairs (''adj'' is mandatory)',...
        mfilename, mfilename);             
  
      % Create input parser and initialize the object
      p = inputParser;
      p = Utils.inputParserSetup(p);
      p.addParameter('adj', []);  
      p.addParameter('inc', []);
      p.addParameter('bus', []);
      p.addParameter('gen', [], @(x) isempty(x) || islogical(x));
      p.addParameter('coh', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,1)==2));
      p.addParameter('ml', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==2));
      p.addParameter('scal', 1, @(x) Utils.isdbl(x) && isscalar(x) && x>0);
      p.addParameter('merge_map', []);
      p.addParameter('powgen', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==2));
      p.addParameter('powlod', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==2));
      p.addParameter('vmag', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==1 && all(x>=0)));
      p.addParameter('vnom', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==1 && all(x>=0)));
      p.addParameter('basmva', [], @(x) isempty(x) || (Utils.isdbl(x) && isscalar(x) && x>0));
      p.addParameter('vw', []);
      p.addParameter('vang', [], @(x) isempty(x) || (Utils.isdbl(x) && size(x,2)==1));
      p.parse(varargin{:});
      varinp = p.Results;
      trash = p.Unmatched;
      assert(isempty(fieldnames(trash)),...
        [ mfilename, ':WrongKeyValueInput'],...
        ['[%s] Some unexpected key value pairs are detected. Please ',...
        'check your spelling. The legitimate keys for the constructor ',...
        'are: adj, inc, bus, gen, coh, ml, scal, merge_map, powgen, ',...
        'powlod, vmag, vang, vnom, vw'], mfilename);
      
      % For now, don't check too much the mutual sizes, what is required
      % and what is not, as well as the datatypes
      obj.scal = varinp.scal;
      if ~isempty(varinp.adj)
        obj.adj = varinp.scal*varinp.adj;
      else
        % error([ mfilename, ':RequiredInputIsMissing'],...
        %   ['[%s] %s.adj is a required property. Please add a meaningful ',...
        %   '''adj'' value to the constructor input'], mfilename, mfilename);
        dflt_adj = tril(ones(5,5), -1);
        dflt_adj = dflt_adj + dflt_adj';
        obj.adj = dflt_adj;  % no error to allow PFgraph()-syntax for object array allocation
      end
      if ~isempty(varinp.inc)
        obj.inc = varinp.inc;  % an unnecessary freedom, but user may put its own inc (it will be checked!)
      else
        obj.inc = GraphUtils.adj2inc( obj.adj );
      end
      obj.gen = varinp.gen;
      obj.coh = varinp.coh;
      obj.ml = varinp.ml;
      obj.merge_map = varinp.merge_map;
      obj.powgen = varinp.powgen;
      obj.powlod = varinp.powlod;
      obj.vmag = varinp.vmag;
      obj.vang = varinp.vang;
      obj.vnom = varinp.vnom;
      obj.basmva = varinp.basmva;
      obj.vw = varinp.vw;
            
      assert(GraphUtils.cmpIncAdj(obj.inc, obj.adj),...
        [ mfilename, ':IncompatibleProperties'],...
        ['[%s] The provided incidence and adjancency matrices (somehow) ',...
         'don''t agree with each other'], mfilename);
             
      siz_adj = size(obj.adj, 1);
      if isempty(varinp.merge_map)
        obj.merge_map = 1:1:siz_adj;  % the map property is set by default to a 1:1 mapping
      else
        obj.merge_map = varinp.merge_map;
      end
      if isempty(varinp.vw)
        obj.vw = sum(obj.adj, 2);
      else
        obj.vw = varinp.vw;  % the vertex weights property is set by default to vertex degrees
      end
      if isempty(varinp.bus)
        obj.bus = 1:1:siz_adj;  % set bus 'names' equal to numbers of corresp. graph vertices
      else
        obj.bus = varinp.bus;
      end
    end
       
    function g = createLaplacian(g, nm)
      %
      % A wrapper around the more general static method GraphUtils.createLaplacian()
      % 
      % Author: Ilya Tyuryukanov
      % Date of first version: 18 August 2015
      % Last revision: 20 January 2015
      
      % Compute Laplacian 
      % As the adj matrix for the object is initialized only once by the 
      % constructor, it doesn't make sense to recompute a Laplacian entry 
      % after it was computed for the first time (add a test for this?)            
            
      switch nm
        case 'rwk'
          if isempty(g.Lrwk)
            L = GraphUtils.createLaplacian(g.adj, 'rwk', g.vw);
            g.Lrwk = L;
          end
        case 'sym'
          if isempty(g.Lsym)
            L = GraphUtils.createLaplacian(g.adj, 'sym', g.vw);
            g.Lsym = L;
          end
        case 'off'
          if isempty(g.Loff)
            L = GraphUtils.createLaplacian(g.adj, 'off', g.vw);
            g.Loff = L;
          end
        otherwise
          error([mfilename,':WrongPublicInput'],...
            ['[%s] Invalid Laplacian normalisation method in %s.createLaplacian(g, nm). ',...
            'Use nm = ''off''/''rwk''/''sym''.'], mfilename, mfilename);
      end
    end   
        
    function obj = set.adj(obj, adj)
      % Setter and input checker for PFgraph.adj
      assert(Utils.isdbl(adj)&&isreal(adj)&&~any(diag(adj)),... 
       [ mfilename, ':InvalidProperty'],...  %&&issymmetric(adj)
       ['[%s] %s.adj should be a double, real, symmetric matrix with zero ',... 
       'diagonal.'], mfilename, mfilename);
      assert(isequal(logical(adj),logical(adj')),...  % relevant if non-symmetric matrices are allowed
        [ mfilename, ':InvalidProperty'],...
       '[%s] The adj-matrix should posess a symmetrical non-zero structure', mfilename, mfilename);      
      % S = GraphUtils.alecconncomp(adj);
      % assert(S==1, [ mfilename, ':InvalidProperty'],...
      %   '[%s] %s.adj should represent a connected graph.', mfilename, mfilename);
      obj.adj = adj;
    end    
    
    function obj = abs_adj(obj)
      % Sets the adjacency matrix (a private property) to abs(obj.adj).
      % This is needed because adj may contain negative entries when it is
      % constructed as a DC power flow admittance matrix for DC power flow  
      % calculations such as MILP islanding. But for some other experiments
      % it may be required to have a compeletely positive adjacency matrix. 
      obj.adj = abs(obj.adj);            
    end
    
    function obj = set.coh(obj, coh)
      % Setter and input checker for PFgraph.coh
      assert(isempty(coh)||(Utils.isint(coh)&&all(coh(1,:)>0)&&size(coh,1)==2),...
        [ mfilename, ':InvalidProperty'],...
        ['[%s] %s.coh should be a matrix of positive integers with two rows ',...
         'or an empty matrix'], mfilename, mfilename);      
      obj.coh = coh;
    end    
    
    function adj = get.adj(obj)
      % Scales adj back to the original
      adj = obj.adj * 1/obj.scal;
    end
  end
  
end

