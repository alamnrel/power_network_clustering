function loadpst(pststate)
% Restores the previosly saved state of pst variables

fields = fieldnames(pststate);
fldnmbr = numel(fields);
PST.pst_var;
for i = 1:1:fldnmbr
 assignin('base',fields{i},pststate.(fields{i}))
end
end
