old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath('..', fullfile('..','..'), fullfile('..','..','third_party'));

load('timing_all_nodes.mat');
num_nod = timing(:,end);
timing(:,end) = [];
num_net = size(timing,1);
num_clus = size(timing,2);
timing([3,5],:) = [];
num_nod([3,5],:) = [];

hold on;
plot(num_nod, timing(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(num_nod, timing(:,7), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 10);

set(gca, 'XTick', num_nod);
xlim([num_nod(1) num_nod(end)]);
ylim([0 2.75]);
%set(gca, 'yscale', 'log'); 
x_width = 241; y_width = 161;
set(gcf, 'pos', [100 100 1.75*x_width 1.5*y_width]);
xlabel('Number of nodes, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('Execution time, [s]', 'Interpreter', 'Latex', 'FontSize', 12);
legend({'k=3', 'k=8'}, 'Box', 'on', 'Location', 'Best');
set(gca, 'xgrid', 'on'); 
set(gca, 'ygrid', 'on'); 