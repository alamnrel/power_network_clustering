function pststate = savepst()
% Saves the current state of pst variables.

PST.pst_var;
varlist = whos('global');
varnmbr = numel(varlist);
for i = 1:1:varnmbr
  var = varlist(i);
  pststate.(var.name) = eval(var.name);
end

end