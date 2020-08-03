function graphviz(path, gv_command)
%
% This function executes a graphviz command gv_command on all .dot files in
% the directory under path and produces .svg graphical files with the same
% names as the respective .dot files.
% In case no gv_command has been provided, gv_command will be set to 
% 'sfdp -Gsize=50! -Goverlap=prism '.
%
% Author: Ilya Tyuryukanov
% Date: 15 September 2016
% Last revision: 15 September 2016

if ispc
  assert( ischar(path) && strcmp(path(2), ':'),...
      [mfilename,':WrongPublicInput'],...
      ['[%s] The 1st input should be a full path to the foldr of in ',...
      'binary index vectors of partitions in its rows.'], mfilename);
else
  assert( ischar(path) && strcmp(path(1), '/'),...
      [mfilename,':WrongPublicInput'],...
      ['[%s] The 1st input should be a full path to the foldr of in ',...
      'binary index vectors of partitions in its rows.'], mfilename);
end

old_pwd = pwd;
cleanup = onCleanup(@() cd(old_pwd));
cd(path);

if nargin==1
  gv_command = 'neato -Gsize=50! -Goverlap=prism '; %sfdp
end

files = dir('*.dot');
system(['cd ', path(1:2), '1>NUL 2>NUL']); % go to the target local disk
system(['cd ', path(3:end)]); %, '1>NUL 2>NUL'
for i=1:length(files),
    dot_file = files(i).name;
    dot_file(end-3:end) = [];
    full_command = sprintf([gv_command, '-Tsvg %s.dot > %s.svg '],...
      dot_file, dot_file); % 1>NUL 2>NUL
    system(full_command);
end


end


