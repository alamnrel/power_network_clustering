function compileDir_simple(Cdir)
if nargin<1
    Cdir=pwd;
end

files = dir(fullfile(Cdir,'*.c*'));

oldDir=pwd;
cd(Cdir);
for j=1:length(files)
    try
%         cm = sprintf('mex %s',files(j).name);
        if ispc
          cm = sprintf('mex -g -largeArrayDims %s',files(j).name);
        else
          cm = sprintf('mex -v GCC=''/usr/bin/gcc-4.9'' -largeArrayDims %s',files(j).name);
        end
        disp(cm);
        eval(cm);
    catch
        disp(lasterr);
        disp('IGNORE if the file is a C++ file which is not a mex file (ie without a mexFunction inside)');
    end
end

cd(oldDir);