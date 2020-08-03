function pst = coh_PST(pst, ns, meth, tol)
% 
% Syntax: pst = coh_PST(pst, ns, meth, tol)
% 
% Purpose:
% Add the coherency information field to a PST-formatted structure.
%
% Input:
% pst: a pst-formatted struct
%   It should contain  bus, line and gen data for the input network in the 
%   PST format (with the solution for the load flow).
% ns:  the desired number of groups.
% meth: slow coherency calculation method (positive whole double)
%   1 for slow coherency, 2 for tight coherency
% tol: (optional) 'tolerance' for coherency in methods 2 and 3.
%
% Output:
% pst: a pst-formatted struct
%   Updated input object with the new field 'coh' (double matrix)
%
% Author: Ilya Tyuryukanov
% Date: 24 July 2015
% Previous revision: 08 January 2016, 06 June 2016
% Last revision: 14 Jan 2020 (decoupling from the original PST files)

assert(ns >= 2, [ mfilename, ':WrongInput'],...
  ['[%s] The ns input should be greater or equal than 2 ',...
  '(even if it is to be overridden by the method, depending on the ',...
  'value of its second input).'], mfilename);

if nargin < 3
  meth = 1;
else
  assert(isscalar(meth) && Utils.isint(meth) && meth>0 && meth<4,...
    [ mfilename, ':WrongPublicInput'],...
    '[%s] %s.%s(pst, ns, meth, tol) only accepts a whole number between ',...
    '1 and 3 as its third input', mfilename, 'PST', mfilename);
end

% Default values for tol are from coh_np16.m in coherency toolbox folder
if nargin < 4
  switch meth
    case 1
      tol = 0;
    case 2
      tol = 0.95;
    case 3
      tol = 0.05;
    otherwise
      error('coh_PST.m: Unsupported coherency method, choose from meth = 1..3');
  end
end
fo = pst.f0;

% Save the current state of PST globals
pststate = PST.savepst();

% Set up global variables of the PST
clearvars -global;
PST.pst_var
bus = pst.bus;
line = pst.lin(:,1:10); % data in lin(:,11:12) doesn't belong to the original PST format, 'modified' PST format is described in mp2pst.m
mac_con = pst.gen;
basmva = pst.basmva;
basrad = 2*pi*fo;  % system frequency in rad/s
syn_ref = 1;  % machine 1 (arbitrary) is reference

% Apply functions from Coherency Toolbox
switch meth
  case 1
    [V_s, lambda] = PST.pstVslow(pst, ns);
    [area,nmach_a] = PST.L_group(V_s);
  case 2
    [V_s, lambda] = PST.pstVslow(pst, ns);
    coh_grp = PST.coh_loose(V_s, tol); %
    grps = tabulate(coh_grp(2,:));
    ngrp = size(grps,1);
    nmax = max(grps(:,2));
    area = zeros(ngrp,nmax);
    nmach_a = grps(:,2);
    for i = 1:1:ngrp
      ind = coh_grp(2,:)==grps(i,1);
      area(i,1:nmach_a(i)) = coh_grp(1,ind);
    end
  otherwise
    error('coh_PST.m: Unsupported coherency method, choose from meth = 1..3');
end
nCA = size(area,1);  % number of Coherent Areas
nmax = size(area,2);
genID = mac_con(:,1);
busID = mac_con(:,2);
coh = zeros(nmax, nCA);
for i = 1:nCA
  cohIDX = ismember(genID, area(i,1:nmach_a(i)));
  busCoh = busID(cohIDX);
  coh(1:nmach_a(i),i) = busCoh;
end

% (robust implementation) Delete all possible and "impossible" non-unique buses
% "impossible": if generators at one bus go to different coherent groups -> columnwise unique() isn't enough
tmp = coh(:);
[~, i_uniq, ~] = unique(tmp);
i_rep = setdiff(1:1:numel(tmp), i_uniq);
tmp(i_rep) = 0;
coh = reshape(tmp, size(coh));

% Delete the possible groups with no generators
empty_groups = all(coh==0,1);
coh(:, empty_groups) = [];

% Write the output
pst.coh = coh;

% Restore the original PST state
PST.loadpst(pststate);
end

