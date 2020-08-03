function vca_case68()
% 
% SVC case studies on the IEEE 68 bus test system.
% 

% Amend MATLAB path
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath('netw_data',fileparts(fileparts(fileparts(mfilename('fullpath')))),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'zp_clustering'),...
  fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'third_party', 'matlab_bgl'));

% Disable selected warnings
old_warn = warning;   % save warning state
cleanup0 = onCleanup(@() warning(old_warn));
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:eigs:SigmaNearExactEig');
warning('off', 'BaseIn:InvalidFolder');
warning('off', 'max_flow_mex:cutValueNotFlowValue');

debug = true; 
rng('default'); 
study = MatpowerIn('caseid', 'case68', 'n_pf', 1, 'f0', 60,...
     'dev_pq', 0, 'mean_pq', 1);
mpopt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1);
x_width = 240; y_width = 160;

% Iterate through cases
for obj = study
  while obj.curr <= 1
    mpw = mp2bgl( obj, 'lsc', 'var' );
    if strcmp(obj.caseid(1:6), 'case68')
      mpw.gen(:,4) =  1000;  % to avoid neg Qgen in CSVC simulation (this case is not coded for simplicity)
      mpw.gen(:,5) =  -1000;
    end
    SET_PQG = mpw.gen( :, 1 );
    SET_LOD = setdiff( mpw.bus(mpw.bus(:,2)==1,1), SET_PQG );   %(PQ-PQG)
    param.num_clu = 2:numel(SET_PQG)+3;
    param.step = 0.1;
    param.balcost = 1e-5;
    param.dvmetr = 'RMSE';
    
    % Create graph for viz & connectivity
    pst = mp2pst(mpw);
    pst.coh = [];
    g = pst2graph(pst);
    ADJ = double(logical(g.adj));
    
    % Obtain VarCS clusters:
    
    param.tst_dst = 1;
    param.sen_thr = 0.01;
    SET_PQG = SET_PQG(randperm(length(SET_PQG)));
    SET_LOD = SET_LOD(randperm(length(SET_LOD)));
    [J_GL, dQgdVg, dVldQl] = VCA.csvc_sen(mpw, true);
    [J_GL] = VCA.cropS(J_GL, SET_PQG, SET_LOD);
    [S_GL] = VCA.makeS(mpw, {'G'}, 1);
    [S_GL] = VCA.cropS(S_GL, SET_PQG, SET_LOD);
    param.num_clu = 2:16;
    [LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = VCA.vcs_sc(J_GL, [], param, ADJ, true);    % vcs_sc on dVl/dVg
    %}
    
    % Here they are as they should be for 10 zones
    
    LBL1 = [16 18 49 50 51];
    LBL2 = [10 11 30 31 32 33 38 46 47 48 53];
    LBL3 = [12 13 17 34 35 36 39 43 44 45 61];
    LBL4 = [15 42];
    LBL5 = [14 40 41];
    LBL6 = [9 26 28 29];
    LBL7 = [4 5 19 20];
    LBL8 = [6 7 21 22 23 24 67 68];
    LBL9 = [1 8 25 27 37 52 54 55];
    LBL10 = [2 3 56 57 58 59 60 62 63 64 65 66];
    LBLA = {LBL1,LBL2,LBL3,LBL4,LBL5,LBL6,LBL7,LBL8,LBL9,LBL10};
    %}
    
    % 25% reactive load increase:
    mpw0 = mpw;   % the reference values
    idx_lod =  mpw.bus(:,3)~=0 | mpw.bus(:,4)~=0;   % if any P or Q demand
    disturbP = zeros(nnz(idx_lod),1);
    disturbQ = abs(0.25*mpw.bus(idx_lod,4));
    mpw.bus(idx_lod,3) = mpw.bus(idx_lod,3) + disturbP;
    mpw.bus(idx_lod,4) = mpw.bus(idx_lod,4) + disturbQ;
    [mpw, xx] = runpf(mpw, mpopt);
    qc_lvl0 = mpw0.gen(ismember(mpw0.gen(:,1),SET_PQG),3)./mpw0.gen(ismember(mpw0.gen(:,1),SET_PQG),4);
    qc_ini = std( qc_lvl0 );   % before disturbance
    qc_lvl1 = mpw.gen(ismember(mpw.gen(:,1),SET_PQG),3)./mpw.gen(ismember(mpw.gen(:,1),SET_PQG),4);
    qc_dst = std( qc_lvl1 );   % after disturbance
    SET_PQG_0 = SET_PQG;
    
    % Our choices:
    PLT_OUR_2 = [18 32 35 42 41 28 20 21 54 60]; SET_PQG = setdiff(SET_PQG_0,PLT_OUR_2); PLT_OUR_2 = setdiff(PLT_OUR_2,SET_PQG_0);
    [mpw_our, dv_our_2, qlvl_our_2, obj_our_2] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_2, param, false);
    qc_our_2 = std(qlvl_our_2);
    %}
    
    % Pilot bus refinement from an initial choice
    PLTB = zeros(1,numel(LBLA));
    % PLT0 = [50 32 35 42 41 28 20 21 55 60];  % 94
    PLT0 = [50 32 35 42 40 28 20 21 54 60];  % 94
    for k = 1:numel(LBLA)
      LBL = LBLA{k};
      dv_opt = Inf;
      x_best = LBL(1);
      for x = LBL
        PLT_OUR_0 = PLT0; PLT_OUR_0(k) = x; SET_PQG = setdiff(SET_PQG_0,PLT_OUR_0); PLT_OUR_0 = setdiff(PLT_OUR_0,SET_PQG_0);
        [mpw_our, dv_our_0, qlvl_our_0, obj_our_0] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_0, param, true);
        qc_our_0 = std(qlvl_our_0);
        if dv_our_0(1)<dv_opt
          dv_opt = dv_our_0(1);
          x_best = x;
          PLTB(k) = x;
        end
      end
    end
    %}
    
    fprintf('%d power flow cases were processed\n', obj.curr);
    obj.curr = obj.curr + 1;
  end
  
end
