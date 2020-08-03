%
% CGP.m does all the work. hMetis is not available under Windows.
% Therefore, the third input to CGP() should be [1,0,1] or [0,0,1] under
% Windows or Mac.
%
% Reference: I. Tyuryukanov, A.C. Karagiannis, M. Popov, M.A.M.M. van der 
% Meijden, and V. Terzija. "Generator grouping cutset determination based 
% on tree construction and constrained spectral clustering". In: The 
% Journal of Engineering 2018.15 (2018), pp. 1309-1314.
%

clc;
clear all;

% Disable warnings
old_warn = warning;   % save warning state
warning('off');   % is optional..
cleanup0 = onCleanup(@() warning(old_warn));

% Amend MATLAB path
old_matlabpath = path;
cleanup1 = onCleanup(@() path(old_matlabpath));
addpath('..', fullfile('..','..'), fullfile('..','..','third_party'),...
  fullfile('..','..','third_party','balancedKCuts_V1_1'),...
  fullfile('..','..','third_party','balancedKCuts_V1_1','Ncut_9'),...
  fullfile('..','..','third_party','balancedKCuts_V1_1','GraphDemos'),...
  fullfile('..','..','third_party','gaimc'));  
clc;
javaaddpath(pwd);   % for java.util.PriorityQueue

% INITIALIZE TESTCASES PARAMETERS
n_ovrl = 2050;  
n_pf = n_ovrl;
n_coh = ones(1,n_pf);
f0 = 60;

% Do the study for the MATPOWER cases
cases = {'case118', 'case1354pegase', 'case2869pegase'};  
for i = 1:1:numel(cases)  
  fprintf('CURRENT CASE NAME: %s\n', cases{i});  
  obj = MatpowerIn('caseid', cases{i}, 'n_pf', n_pf, 'n_coh', n_coh, 'f0',...
    f0, 'dev_pq', 0.5*ones(1,n_pf), 'mean_pq', ones(1,n_pf), 'seed', 2);      
  mpw = mp2bgl(obj);
  pst = mp2pst(mpw);
  pst.gen(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model
  if isempty(pst.bus), error('AN ERROR.'); end
  case_stat = CGP(pst, 2:1:18, [1,0,1]);
  for k = 1:1:size(case_stat,1)
    fprintf('%4d ', case_stat(k,:));
    fprintf('\n');
  end
end

% Do the study for the 48 machine NPCC test system
fprintf('CURRENT CASE NAME: %s\n', 'case140_npcc'); 
load('datanp48.mat');
% Prepare the data of PST to be consistent with PFgraph
[bus(:,1),I] = sort(bus(:,1),'ascend');  % sort with respect to busIDs to ensure that... 
 bus(:,2:end) = bus(I,2:end);  
[mac_con(:,2),I] = sort(mac_con(:,2),'ascend');  % busIDs in bus[] and gen[] are in the same order
ncol_g = size(mac_con,2);
 mac_con(:,[1,3:ncol_g]) = mac_con(I,[1,3:ncol_g]);
tol = 2e-13; iter_max = 30; 
vmin = 0.5;  vmax = 1.5; acc = 1.0; 
[DUMMY, bus, line_flw] = evalc('PST.LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,''n'',2);'); 
line(:,12) = eps;   % dummy line flow
bus(:,13) = 1;      % dummy nominal voltage
mac_con(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model
pst = struct_pst(bus, line, mac_con, 100, [], f0);   % 100 MVA base in PST toolbox coh_np48.m 
case_stat = CGP(pst, 2:1:18, [1,0,1]);
for k = 1:1:size(case_stat,1)
  fprintf('%4d ', case_stat(k,:));
  fprintf('\n');
end

% Do the timing study
cases = {'case39', 'case300', 'case118', 'case1354pegase', 'case2869pegase'};
clk = zeros(2, numel(cases));
for i = 1:1:numel(cases)  
  fprintf('CURRENT CASE NAME: %s\n', cases{i});  
  obj = MatpowerIn('caseid', cases{i}, 'n_pf', n_pf, 'n_coh', n_coh, 'f0',...
    f0, 'dev_pq', 0.5*ones(1,n_pf), 'mean_pq', ones(1,n_pf), 'seed', 2);      
  mpw = mp2bgl(obj);
  pst = mp2pst(mpw);
  pst.gen(:,[4:6,8:15,18,20:21]) = 0;  %ensure classical model
  if isempty(pst.bus), error('AN ERROR.'); end
  [case_stat, clk(:,i)] = CGP(pst, 5, [0,0,1]);  
end
plot([118, 300, 1354, 2383, 2869], clk(1,:), '.--', 'Color', [0.7,0.7,0.7], 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
plot([118, 300, 1354, 2383, 2869], clk(2,:), '.-', 'Color', [0.7,0.7,0.7], 'LineWidth', 1.5, 'MarkerSize', 10);
plot([118, 300, 1354, 2383, 2869], clk(1,:)+clk(2,:), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 10);
set(gca, 'XTick', [118, 300, 1354, 2383, 2869]);
xlim([118 2869]);
x_width = 241; y_width = 161;
set(gcf, 'pos', [100 100 1.75*x_width 1.5*y_width]);
xlabel('Number of nodes, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('Execution time, [s]', 'Interpreter', 'Latex', 'FontSize', 12);
legend({'Distance Graph', 'Tree heuristic', 'Total'}, 'Box', 'on', 'Location', 'Best');
set(gca, 'xgrid', 'on');
set(gca, 'ygrid', 'on');
% set(gca, 'yscale', 'log'); 
% ylim([0.01 20])
% set(gca,'YTick',[0.01,0.1,1,20])

