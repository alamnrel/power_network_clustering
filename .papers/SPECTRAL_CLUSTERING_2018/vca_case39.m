function vca_case39()
% 
% SVC case studies on the IEEE 39 bus test system.
% 
% Reference: I. Tyuryukanov, M. Popov, M.A.M.M. van der Meijden, and V. 
% Terzija. "Discovering Clusters in Power Networks From Orthogonal  
% Structure of Spectral Embedding". In: IEEE Trans. Power Syst. 33.6 
% (2018), pp. 6441-6451.
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
study = MatpowerIn('caseid', 'case39_0', 'n_pf', 1, 'f0', 60,...
     'dev_pq', 0, 'mean_pq', 1);
mpopt = mpoption('out.all', 0, 'verbose', 0, 'out.suppress_detail', 1);
x_width = 240; y_width = 160;
   
% Iterate through cases
for obj = study
  while obj.curr <= 1 
    
    mpw = mp2bgl( obj, 'lsc', 'var' );
    if strcmp(obj.caseid(1:6), 'case39')
      isw = find(mpw.bus(:,2)==3);
      mpw.bus(isw, 2) = 2;   
      mpw.bus( 39, 2) = 3;   % set the slack at 39!
      gsw = find(mpw.gen(:,1)==39);
      mpw.gen(gsw,4) = 1.1*mpw.gen(gsw,4);   % this is just necessary for 25% load increase
      qps = mpw.gen(:,5)>0;
      mpw.gen(qps,5) = 0;   % high bottom limits of Qgen are too strict! 
    end
    SET_PQG = mpw.gen( :, 1 ); 
    SET_PQG = setdiff( SET_PQG, 39 );
    SET_LOD = setdiff( mpw.bus(mpw.bus(:,2)==1,1), SET_PQG );   % (PQ - PQG)
    param.num_clu = 2:10;
    param.step = 0.1;
    param.balcost = 1e-5;
    param.dvmetr = 'RMSE';    
    
    % Create graph for viz
    pst = mp2pst(mpw);
    pst.coh = [];
    g = pst2graph(pst);
    ADJ = double(logical(g.adj));    
    
    % Conejo-1993 pilots study
    %{
    [S_GL, S_QG, S_LL] = VCA.csvc_sen(mpw, debug);
    [S_GL] = VCA.cropS(S_GL, SET_PQG, SET_LOD);
    [S_LL] = VCA.cropS(S_LL, SET_LOD, SET_LOD);
    i_lod = ismember(mpw.bus(:,1), S_LL.busrow);
    q_lod = abs(mpw.bus(i_lod,4));
    q_lod(q_lod==0) = min(nonzeros(q_lod))/100;
    C_LL = diag(q_lod);
    Q_x = C_LL;
    [plt_set, plt_qlt] = VCA.pilots_conejo1993(S_GL, S_LL, C_LL, Q_x, [5,12,16,20], 2, true);    
    %}
    
    % Run generator-power-level balancing OPF (or the base-case is a good reference?)
    
    idx_pqg = ismember(mpw.gen(:,1), SET_PQG);
    mpw = VCA.set_gencosts_svc(mpw, idx_pqg);          
    mpw = runopf(mpw, mpopt);         
    %}
    
    % Test various zoning methods.....
    %{
    param.tst_dst = 1;
    param.sen_thr = 0.01;
    [J_GL, dQgdVg, dVldQl] = VCA.csvc_sen(mpw, true);
    [J_GL] = VCA.cropS(J_GL, SET_PQG, SET_LOD);
    [LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = VCA.vcs_sc(J_GL, [], param, ADJ, debug);    % vcs_sc on dVl/dVg
    [S_GL] = VCA.makeS(mpw, {'G'}, 1);
    [S_GL] = VCA.cropS(S_GL, SET_PQG, SET_LOD);
    [LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = VCA.vcs_sc(S_GL, [], param, ADJ, debug);    % vcs_sc on dVl/dQg
    [LBL, avg_clu_dst] = VCA.vcs(S_GL, 'average');               % vcs on dVl/dQg
    [S_LL] = VCA.makeS(mpw, {'L'}, 1);
    [S_LL] = VCA.cropS(S_LL, SET_LOD, SET_LOD); 
    [dVldQl] = VCA.cropS(dVldQl, SET_LOD, SET_LOD);
    assert( norm(S_LL.S-dVldQl.S,'fro')<1e-2 );                  % compare dVl/dQl with 2 methods
    [LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = VCA.vcs_sc(S_GL, S_LL, param, ADJ, debug);    % vcs_sc on both dVl/dQg AND dVl/dQl     
    for k = 2:1:10
      param.num_clu = k;
      [LBL1, ~, ~, ~, ~, eXp1] = VCA.vcs_sc(J_GL, [], param, ADJ);
      [LBL2, ~, ~, ~, ~, eXp2] = VCA.vcs_sc(S_GL, [], param, ADJ);
      LBL1 = cellfun(@(x)x(:)',LBL1,'un',false);
      LBL2 = cellfun(@(x)x(:)',LBL2,'un',false);
      for j = 1:1:numel(LBL1)
        fprintf('{%d: ', j);
        fprintf('%d ', LBL1{j});
        fprintf('%-1s, ', '}');
      end
      fprintf('\n');
      for j = 1:1:numel(LBL2)
        fprintf('{%d: ', j);
        fprintf('%d ', LBL2{j});
        fprintf('%-1s, ', '}');
      end
      fprintf('\n');            
    end
    fprintf('\n');
    [plt, plt_qt] = VCA.pilots(S_LL, LBL, 'strongest_bus', 3);
    [plt, plt_qt] = VCA.pilots(S_LL, LBL, 'max_zone_sens', 3);
     plt = cellfun(@(x) x(1), plt, 'un', true);
    %}
    
    % Case study #1: Automatic clusters and their testing by Q up +25%...
    [~, dQgdVg, dVldQl, S_GL] = VCA.csvc_sen(mpw, true);   % agreed that dVi/dVj is the best!
    [S_GL] = VCA.cropS(S_GL, SET_PQG, SET_LOD);         
    [LBL, LBL_GEN, LBL_LOD, RST, CQ1, eXp] = VCA.vcs_sc(S_GL, [], param, ADJ, false);
    [LBL, T] = VCA.vca_rest(g.adj, g.bus, LBL);        
    
    % 25% reactive load increase:
    mpw0 = mpw;   % the reference values
    % idx_lod = mpw.bus(:,1)==12;
    % disturbP = 0.5*mpw.bus(idx_lod,3);     
    % disturbQ = 0.5*mpw.bus(idx_lod,4);   
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
       
    % The final choices
    debug = false;
    PLT_OUR_0 = [ 5, 16, 20, 28, 3 ];
    [mpw_our, dv_our_0, qlvl_our_0, obj_our_0] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_0, param, true); 
    qc_our_0 = std(qlvl_our_0);    
    PLT_OUR_2 = [ 5, 16, 20, 28, 3, 12 ];
    [mpw_our, dv_our_2, qlvl_our_2, obj_our_2] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_2, param, debug);  
    qc_our_2 = std(qlvl_our_2);        
    % Now test with pilots from Conejo (set 3 + 1)
    PLT_CNJ_5 = [ 5, 12, 16, 20, 28 ];
    [mpw_cnj_5, dv_cnj_5, qlvl_cnj_5, obj_cnj_5] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_CNJ_5, param, debug);    
    qc_cnj_5 = std(qlvl_cnj_5);        
    [mpw_min_3, dv_min_3, qlvl_min_3, obj_min_3] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, [5,16,28], param, debug);    
    qc_min_3 = std(qlvl_min_3);      
    
    SET_PQG = [SET_PQG; 39];
    % And pilots from Sun-Guo:
    PLT_GUO_0 = [1, 3, 6, 24, 28];    
    [mpw_guo_0, dv_guo_0, qlvl_guo_0] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_0, param, true);
    qc_guo_0 = std(qlvl_guo_0);     
    PLT_GUO_2 = [1, 3, 5, 12, 24, 28];    
    [mpw_guo_2, dv_guo_2, qlvl_guo_2] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_2, param, debug);
    qc_guo_2 = std(qlvl_guo_2);      
    PLT_GUO_1 = [1, 3, 6, 24, 28, 20];    
    [mpw_guo_1, dv_guo_1, qlvl_guo_1] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_1, param, debug);
    qc_guo_1 = std(qlvl_guo_1);     
    PLT_GUO_3 = [1, 3, 5, 12, 24, 28, 20];    
    [mpw_guo_3, dv_guo_3, qlvl_guo_3] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_3, param, debug);
    qc_guo_3 = std(qlvl_guo_3);      
    SET_PQG(end) = [];
    
    % Case study #3: a topological change!
    mpw = mpw0;
    lin = mpw.branch(:,1:2);
    lin_out = (lin(:,1)==6 | lin(:,2)==6)&(lin(:,1)==11 | lin(:,2)==11);
    mpw.branch(lin_out,:) = [];
    %lin_out = (lin(:,1)==4 | lin(:,2)==4)&(lin(:,1)==14 | lin(:,2)==14);
    %mpw.branch(lin_out,:) = [];
    [mpw, sx] = runpf(mpw, mpopt);   % re-running pf is mandatory    
    pst = mp2pst(mpw);   % update graph for viz
    pst.coh = [];
    g = pst2graph(pst);     
    [~, dQgdVg, dVldQl, S_GL] = VCA.csvc_sen(mpw, true);  % agreed that dVi/dVj is the best!
    [S_GL] = VCA.cropS(S_GL, SET_PQG, SET_LOD);         
    % param.num_clu = 2:8;
    [LBL, LBL_GEN, LBL_LOD, RST, CQ2, eXp] = VCA.vcs_sc(S_GL, [], param, ADJ, debug);
    [LBL, T] = VCA.vca_rest(logical(g.adj), g.bus, LBL);    
    if debug && exist('CQ1', 'var') && exist('CQ2', 'var')
      figure;
      plot(1:numel(CQ1), CQ1, 'rs-', 'MarkerSize', 8, 'LineWidth', 1.5);
      hold on;
      plot(1:numel(CQ2), CQ2, 'bo-', 'MarkerSize', 8, 'LineWidth', 1.5);
      xlabel('Nr. of eigenvectors, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
      ylabel('$C^*$, [-]', 'Interpreter', 'Latex', 'FontSize', 13);
      set(gcf, 'pos', [100 100 1.5*x_width 1.5*y_width]);
      set(gca, 'XTick', 1:numel(CQ1));
      xlim([1,numel(CQ1)]);
      set(gca, 'XTickLabel', cellfun(@num2str, num2cell(param.num_clu), 'unif', false));
      set(gca, 'xgrid', 'on');
    end
    
    % 25% reactive load increase:
    mpwX = mpw; % the reference values
    mpw.bus(idx_lod,3) = mpw.bus(idx_lod,3) + disturbP; 
    mpw.bus(idx_lod,4) = mpw.bus(idx_lod,4) + disturbQ;
    [mpw, xx] = runpf(mpw, mpopt);
    qc_lvl = mpwX.gen(ismember(mpwX.gen(:,1),SET_PQG),3)./mpwX.gen(ismember(mpwX.gen(:,1),SET_PQG),4);
    qc_ini = std( qc_lvl );  % before disturbance
    qc_lvl = mpw.gen(ismember(mpw.gen(:,1),SET_PQG),3)./mpw.gen(ismember(mpw.gen(:,1),SET_PQG),4);
    qc_dst = std( qc_lvl );  % after disturbance      
    [mpw_our_1, dv_our_7, qlvl_our_7, obj_our_7] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_0, param, true);
    qc_our_7 = std(qlvl_our_7);
    [mpw_our_2, dv_our_8, qlvl_our_8, obj_our_8] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_OUR_2, param, true);
    qc_our_8 = std(qlvl_our_8);
    SET_PQG = [SET_PQG; 39];
    [mpw_guo_4, dv_guo_4, qlvl_guo_4] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_0, param, true);
    qc_guo_4 = std(qlvl_guo_4);        
    [mpw_guo_5, dv_guo_5, qlvl_guo_5] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_2, param, true);
    qc_guo_5 = std(qlvl_guo_5);          
    [mpw_guo_6, dv_guo_6, qlvl_guo_6] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_1, param, true);
    qc_guo_6 = std(qlvl_guo_6);       
    [mpw_guo_7, dv_guo_7, qlvl_guo_7] = VCA.test_vca_plt(mpw0, mpw, SET_PQG, PLT_GUO_3, param, true);
    qc_guo_7 = std(qlvl_guo_7);          
     
    fprintf('%d power flow cases were processed\n', obj.curr);
    obj.curr = obj.curr + 1;    
  end
end


end