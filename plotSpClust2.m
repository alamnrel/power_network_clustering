function [] = plotSpClust2(varargin)
%
% Plot different pictures related to spectral embedding. The type of
% representation is selected based on the input types.


clrs = {rgb('SlateGray'), rgb('ForestGreen'), 'b', rgb('Purple'), rgb('Orange'), 'm', 'c', 'y',  ...  %1-7
  rgb('DarkKhaki'), rgb('Brown'), rgb('Salmon'), 'r', 'k', rgb('Teal'), rgb('Maroon'), rgb('Purple'),... %8-13
  rgb('MidnightBlue'), rgb('Tomato'), rgb('Wheat'), rgb('DarkMagenta'),... % 14-17
  rgb('Azure'), rgb('DeepPink'), rgb('Tan'), 'g', rgb('Plum'),... % 18-22
  rgb('Indigo'), rgb('DarkSlateGrey')};  % 23-26
mrks = {'.', 'x', '*',  'o', 's',  'x',  'd', '^', 'v', '>', '<', 'p', 'h'};

switch length(varargin)
  case 2
    mode = 'eigvec';
    Y = varargin{1};
    if size(Y,2)>size(Y,1)
      Y = Y.';
    end
    busmap = varargin{2};
    
  case 3
    mode = 'genlod';
    Y = varargin{1};
    busmap = varargin{2};
    coh = varargin{3};
    
  otherwise
    disp('plotSpClust.m: The supplied set of arguments is not valid, nothing will be plotted.');
    return;
end

switch mode
  case 'eigvec'
    drawEigVec(Y, busmap);
  case 'genlod'
    drawGenLod(Y, busmap, coh, clrs, mrks);
  otherwise
    error('plotSpClust.m: By any means, we shouldn''t be here!');
end
end

function drawEigVec(Y, busmap)
assert(isa(Y,'double')&&ismatrix(Y)&&isreal(Y),...
  'plotSpClust.m: 1st input should be real double matrix');
assert(isa(busmap,'double')&&ismatrix(busmap)&&isreal(busmap)&&...
  size(busmap,1)<=2,...
  'plotSpClust.m: 2nd input should be real double matrix with 1 or 2 rows');

n = size(Y,2);
m = size(Y,1);
busmap = busmap(1,:);
for i = 1:n
  subplot(n,1,i);
  plot(1:1:m, Y(:,i), 'k.-');
  for j = 1:m
    text(j, Y(j,i), num2str(busmap(1,j)), 'FontSize', 8, 'Color', 'k', 'FontWeight',...
      'normal', 'VerticalAlignment','bottom', 'HorizontalAlignment', 'right');
  end
  ylabel([num2str(i),' Eigenvector'], 'FontSize', 9.5);
  xlabel('Node Number', 'FontSize', 9.5);
end
end

function drawGenLod(Y, busmap, buslbl, clrs, mrks)
assert(isa(Y,'double')&&ismatrix(Y)&&isreal(Y),...
  'plotSpClust.m: 1st input should be real double matrix');
assert(isa(busmap,'double')&&ismatrix(busmap)&&isreal(busmap)&&...
  size(busmap,1)==1,...
  'plotSpClust.m: 2nd input should be real double matrix with 1 row');
assert(isa(buslbl,'double')&&ismatrix(buslbl)&&isreal(buslbl),...
  (size(buslbl,1)==2 || size(buslbl,1)==3),...
  'plotSpClust.m: 3rd input should be real double matrix with 2 or 3 rows');
Y = checkDims(Y);
if isempty(Y), return; end
n_buses = numel(busmap);
area = unique(buslbl(2,:));
n_area = numel(area);
n_clrs = length(clrs);
if n_area>n_clrs
  warning('The number of colors available is less than the number of colors requested, define less areas.');
  return;
end
if size(buslbl,1)>2
  area_mrk = unique(buslbl(3,:));
  n_area_mrk = numel(area_mrk);
  n_mrks = length(mrks);
  if n_area_mrk>n_mrks
    warning('The number of markers available is less than the number of markers requested, only the marker x will be used instead');
    mrks = repmat({'x'},1,n_area_mrk);
  end
end

if size(Y,2) == 1
  [Y, idx] = sort(Y, 'descend');
  buslbl = buslbl(:,idx);
  busmap = busmap(idx);
  numDi = 1;
elseif all(Y(:,1)/Y(1,1) < 1.001) && all(Y(:,1)/Y(1,1) > 0.998)
  numDi = 2;
else
  numDi = 3;
end

figure;
hold on;
std_mark_siz = 75;
std_line_wdth = 3;
for i = 1:n_area
  idx_area = find(buslbl(2,:)==area(i));
  if size(buslbl,1)==2
    switch numDi
      case 1
        plot(idx_area, Y(idx_area,1), '.', 'Color', clrs{i}, 'MarkerSize', 12);
      case 2
        scatter(Y(idx_area,2), Y(idx_area,3), std_mark_siz, clrs{i}, 'x', 'LineWidth', std_line_wdth);
      case 3
        scatter3(Y(idx_area,1), Y(idx_area,2), Y(idx_area,3), std_mark_siz, clrs{i}, 'x', 'LineWidth', std_line_wdth);
    end
  else
    switch numDi
      case 1
        plot(idx_area, Y(idx_area,1), mrks{i}, 'Color', clrs{i});
      case 2
        scatter(Y(idx_area,2), Y(idx_area,3), std_mark_siz, clrs{i}, mrks{i}, 'LineWidth', std_line_wdth);
      case 3
        scatter3(Y(idx_area,1), Y(idx_area,2), Y(idx_area,3), std_mark_siz, clrs{i}, mrks{i}, 'LineWidth', std_line_wdth);
    end
  end
end

% Print bus labels
for i = 1:n_area
  k = buslbl(1,(buslbl(2,:)==area(i)));
  k = ismember(busmap,k);
  k = find(k);
  for j = 1:numel(k)
    switch numDi
      case 1
        xx = find(buslbl(2,:)==area(i));
        text(xx(j), Y(k(j),1), num2str(busmap(1,k(j))), 'FontSize', 11, 'Color', clrs{i},...
          'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');        
      case 2
        text(Y(k(j),2), Y(k(j),3), num2str(busmap(1,k(j))), 'FontSize', 11, 'Color', clrs{i},...
          'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
      case 3
        text(Y(k(j),1), Y(k(j),2), Y(k(j),3), num2str(busmap(1,k(j))), 'FontSize', 15, 'Color', clrs{i},...
          'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    end
  end
end
setView(Y);
xl = xlim; yl = ylim; zl = zlim;
epsilon = eps('double');
[~, ix] = min(abs(xl)); [~, iy] = min(abs(yl)); [~, iz] = min(abs(zl));
sx = sign(xl);
if sx(1) == sx(2), xl(ix) = -epsilon*sx(ix); end
mx = xl(ix);
sy = sign(yl);
if sy(1) == sy(2), yl(iy) = -epsilon*sy(iy); end
my = yl(iy);
% yl = abs(yl).*sy;
sz = sign(zl);
if sz(1) == sz(2), zl(iz) = -epsilon*sz(iz); end
mz = zl(iz);

switch numDi
  case 1
    plot([0 n_buses+1], [0 0], 'k-', 'LineWidth', 2);
    xlim(xl);
  case 2
    plot([-100 100], [0 0], 'k-', 'LineWidth', 2);
    plot([0 0], [-100 100], 'k-', 'LineWidth', 2);
    xlim(xl); ylim(yl);
    text(xl(xl~=mx),0, 'X', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    text(0,yl(yl~=my), 'Y', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
  case 3
    plot3([-100 100], [0 0], [0 0], 'k-', 'LineWidth', 2);
    plot3([0 0], [-100 100], [0 0], 'k-', 'LineWidth', 2);
    plot3([0 0], [0 0], [-100 100], 'k-', 'LineWidth', 2);
    xlim(xl); ylim(yl); zlim(zl);
    text(xl(xl~=mx),0,0, 'X', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    text(0,yl(yl~=my),0, 'Y', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    text(0,0,zl(zl~=mz), 'Z', 'FontSize', 15, 'FontWeight', 'bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end
end


function setView(Y)
twoDi = all(Y(:,1)/Y(1,1) < 1.001) && all(Y(:,1)/Y(1,1) > 0.998);
n = size(Y, 2);
if n==1
  return;
end
if (n == 2) || twoDi
  view(0,90);  % view xOy
else
  view(-37.5, 30);  % default 3D view
  axis equal;
end
grid on;
end

function [Y, n, m] = checkDims(Y)
n = size(Y, 2);
m = size(Y, 1);
if (n<1) || (n>3)
  disp('plotSpClust.m: Plotting of spectral embedding for n>3 or n<1 is not supported.');
  Y = [];
elseif (n==2)
  Y = [ones(m, 1), Y];  % pad the 3rd dimension with 1s (for conformity!)
end
end