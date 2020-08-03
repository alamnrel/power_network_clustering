function main()
%
% This simple method is mostly suitable for offline studies and may benefit
% a lot from repeated runs with different initializations.
% 
% A more robust alternative is in the folder DPSP2018.
%
% Reference: I. Tyuryukanov, J. Quir?s-Tort?s, M. Naglic, M. Popov, M. van 
% der Meijden, and V. Terzija. ?Controlled Islanding of Power Networks 
% based on Graph Reduction and Spectral Clustering?. In: Proc. 10th 
% Mediterranean Conference on Power Generation, Transmission, Distribution 
% and Energy Conversion (MedPower 2016). Belgrade, Serbia, 2016, pp. 1?6.
%

% Terminate all Excel processes
[~, cmdout] = system('taskkill /f /im excel.exe');
disp(['EXCEL process status (in start.m): ''', cmdout(1:end-2), '''']);

% Disable selected warnings
old_warn = warning;  % save warning state
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:eigs:SigmaNearExactEig');
warning('off', 'MATLAB:eigs:IllConditionedA');
cleanup0 = onCleanup(@() warning(old_warn));
% Amend MATLAB path
old_matlabpath = path;
addpath(fullfile('..','..'), fullfile('..','..','third_party'));
cleanup1 = onCleanup(@() path(old_matlabpath));

% INITIALIZE SEVERAL CONSECUTIVE TESTCASES
cases = {'case300', 'case1354pegase', 'case2869pegase', 'case39', 'case2383wp', 'case89pegase'}; 
n_coh = {[3,2], 2, [3;4], [4;2;3], [3;4], [4;2;3]};  % number of coherent generator groups for each case and each power flow (an integer vector for each case)
n_pf = cellfun(@numel, n_coh, 'UniformOutput', true);  % number of sampled power flows for each case
f0 = [50, 60, 50, 60, 50, 50];
pvstruct1 = struct('sc_alg', 'sym');
pvstruct2 = struct(); 
study(1, numel(cases)) = MatpowerIn();
for i = 1:1:numel(cases)
  study(1, i) = MatpowerIn('caseid', cases{i}, 'n_pf', n_pf(i), 'n_coh',...
    n_coh{i}, 'f0', f0(i));
end

% ITERATE THROUGH THE CASES
for obj = study
  %xl = ExcelUtils(obj.excelfile);
  %xl.write_header(obj.caseid);
  while obj.curr <= obj.n_pf
    clearvars -global;  % this cleans up PST        
    mpw = mp2bgl(obj);
    pst = mp2pst(mpw);         
    pst = PST.coh_PST(pst, obj.n_coh(obj.curr), 1, 0.98);
    g = pst2graph(pst);
    if strcmp(obj.caseid, 'case39') && obj.n_coh(obj.curr) == 4 
      g.coh = [30 31 32 33 34 35 36 37 38 39;...
                4  1  1  2  2  3  3  4  4  1];
    end       
    if strcmp(obj.caseid, 'case1354pegase') && obj.n_coh(obj.curr) == 2
      [~,idx,~] = intersect(g.bus(g.coh(1,:)), [1341, 8950, 7049, 4850, 2886, 9108,...
        5237, 8976, 1478, 4952, 2481, 3611, 7327, 4624, 583, 1083, 8676, 5781,...
        6734, 3492, 6516, 6947, 8807, 7697, 5971, 972, 2786, 7776, 7998,...
        3205, 4125, 8564, 5940, 7495, 5564, 7036, 7495, 616, 3346, ...
        2446, 2627, 4024, 9067, 2085, 7056, 5814, 795, 8311,...
        2719, 4880, 5144, 3513, 1043, 2985, 8158, 3240, 6331, 9180, 7159,...
        8795, 4056, 3661, 4419, 5486, 776, 6168, 6888, 5994, 2050, 8522,...
        4128, 1100, 8473, 4331, 3353, 2197, 7115, 851, 5110, 8961, 221,...
        453, 5533, 4231, 6831, 2600, 6807, 8997, 7842, 3218, 2797, 2426,...
        6153, 8903, 1295, 5482, 2550, 6820, 778, 7209, 2934, 2468, 338,...
        3133, 682, 2291, 3656, 8225, 3114, 4823, 5004, 6351, 639, 8043,...
        3436, 1721, 150, 8818, 1808, 6376]);
      g.coh(2,:) = 1;
      g.coh(2,idx) = 2;
    end
    if strcmp(obj.caseid, 'case89pegase') && obj.n_coh(obj.curr) == 2
      g.coh = [7 21 24 37 42 45 56 60 65 71 81 89;...
               1  1  2  1  1  1  2  2  1  2  2  2];
    end
    if strcmp(obj.caseid, 'case89pegase') && obj.n_coh(obj.curr) == 3
      g.coh = [7 21 24 37 42 45 56 60 65 71 81 89;...
               1  1  3  1  3  1  2  2  3  2  2  2];
    end    
    if strcmp(obj.caseid, 'case89pegase') && obj.n_coh(obj.curr) == 4
      g.coh = [7 21 24 37 42 45 56 60 65 71 81 89;...
               1  1  3  1  4  4  2  2  4  2  3  2];
    end     
    [ok, bran_viols, coh_viols] = g.chkBranCoh();
    if ~ok
      g.coh(2,coh_viols(:,1)) = g.coh(2,coh_viols(:,2));  % or 2 <- 1
      g.ml(bran_viols, :) = [];
    end
    trees = coh2sptree(g);    
    f = @() g.coGrReduClust(trees, pvstruct1);
    t_exe1 = timeit(f, 1);
    T1 = g.coGrReduClust(trees, pvstruct1);  % now repeat for results
    cut1 = final_cutset(g, T1);
    [eXp1, pcut1, ~, ~, n_viols1, ~] = g.ici_info( cut1 );
    %xl.write_results(2*obj.curr-1, 1, 'coGrReduClust', pvstruct1.sc_alg, pcut1,...
    %   obj.n_coh(obj.curr), eXp1(1,:), max(eXp1(1,:)), t_exe1, obj.caseid,...
    %   n_viols1, eXp1(2,:), max(eXp1(2,:)), min(eXp1(3,:)));
         
    obj.curr = obj.curr + 1;
  end
  %delete(xl);  % destroy the object, i.e. shutdown Excel
end

end
