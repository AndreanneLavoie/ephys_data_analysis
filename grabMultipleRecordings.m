function [filenamescell] = grabMultipleRecordings(filetype, intanpath)
%INPUT
%if th

if nargin < 1
    filetype = '*amp-A*.dat';
end

directorycell = {};
filenamescell={};

while 1
    directoryname = uigetdir(intanpath, 'Pick IntanData Directory');%change Intan Directory
    if ~directoryname
        break;
    else
        directorycell{end + 1} = directoryname;
        cd(directoryname);
        tempfilenames=dir(filetype);
        filenamescell{end+1} = tempfilenames;
    end
end 
filelength=[];
for i=1:length(filenamescell)
    filelength=[filelength length(filenamescell{i})];
end

if length(unique(filelength))>1
   print('the number of files does not match');
   return;
end


