% 
% Plots outlier statistics from the specified mat-file with results
% (see the results folder).
% 
%
% bad_clu_ou - nr. of undesirable too small clusters (these cannot be fully 
%   precluded by outlier filtering)
% pp_clu_ou  - nr. of outlier clusters (these should be fully precluded)
%

nr1 = 3;
nr2 = 4;

old_matlabpath = path;
cleanup = onCleanup(@() path(old_matlabpath));
addpath('..', fullfile('..','..'), fullfile('..','..','third_party'));

load('ou_case2383wp_classic.mat');
num_tst = size(pp_clu_ou,2);
n_tst = 150;

% Number of small clusters (undesired partitioning outcomes)
x_vec0 = [1:1:num_tst; 1:1:num_tst; 1:1:num_tst];
x_vec0 =  x_vec0(:);
y_vec3 = [zeros(1,num_tst); bad_clu_ou(nr1,:); zeros(1,num_tst)];
y_vec3 = y_vec3(:);
y_vec4 = [zeros(1,num_tst); bad_clu_ou(nr2,:); zeros(1,num_tst)];
y_vec4 = y_vec4(:);
figure;
plot(x_vec0, y_vec4, 'Color', rgb('lightgray'), 'MarkerSize', 2, 'LineWidth', 2);
hold on;
plot(x_vec0, y_vec3, 'Color', rgb('blue'), 'MarkerSize', 2, 'LineWidth', 2);
xlabel('Nr. of power flow case, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('Nr. of small clusters, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
set(gca, 'XTick', [1,n_tst:n_tst:20*n_tst])
xlim([1 num_tst]);
x_width = 241; y_width = 161;
set(gcf, 'pos', [100 100 1.75*x_width 1.25*y_width]);
for i = 1:1:round(num_tst/100)-1
  plot([i*n_tst, i*n_tst], ylim, '-.', 'LineWidth', 0.5, 'Color', [0.2 0.2 0.2]);
end
legend({'HSC', 'HSC-pp'}, 'Box', 'on', 'Location', 'Best');


% Average outlier size:
%{
y_vec1 = [zeros(1,num_tst); av_siz_ou(1,:); zeros(1,num_tst)];
y_vec1 = y_vec1(:);
y_vec3 = [zeros(1,num_tst); av_siz_ou(3,:); zeros(1,num_tst)];
y_vec3 = y_vec3(:);
y_vec4 = [zeros(1,num_tst); av_siz_ou(4,:); zeros(1,num_tst)];
y_vec4 = y_vec4(:);
figure;
plot(x_vec0, y_vec4, 'Color', rgb('orange'), 'MarkerSize', 2, 'LineWidth', 2);
hold on;
plot(x_vec0, y_vec1, 'Color', rgb('cyan'), 'MarkerSize', 2, 'LineWidth', 2);
plot(x_vec0, y_vec3, 'Color', rgb('blue'), 'MarkerSize', 2, 'LineWidth', 2);
xlabel('case number, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('$n_{small}$, [-]', 'Interpreter', 'Latex', 'FontSize', 12);
% legend({'HSC', 'k-means', 'HSC-pp'}, 'Box', 'off', 'Location', 'Best');
for i = 1:1:round(num_tst/100)
  plot([i*100, i*100], ylim, 'k--', 'LineWidth', 1);
end
%}