function [ops,savePath] = mergeIntanOFPCSystem3(savePath)
% the folder containing intan data files should be named after XXXX_ephys; works for one file per channel intan aquisition (OFPC)
% There cannot be any space in the path, otherwise error message "The system cannot find the file specified" will show
% Select files to merge


recording_dir = uigetdir('Pick IntanData Directory');
cd(recording_dir);
channel_dat=dir('*amp-A*.dat');%find all files with *amp-A*.dat format
channel_dat=struct2cell(channel_dat);
channel_dat=string(channel_dat(1,:));%make file name a string array to match following codes

% extract exp_id
fullPath = recording_dir(1:(strfind(recording_dir,'_ephys')+5));  
lastchar=strfind(fullPath,'_ephys')-1;
firstchar = strfind(fullPath,'\');
firstchar=firstchar(firstchar<lastchar);
firstchar = firstchar(end)+1;
exp_id = fullPath(firstchar:lastchar);  

% fPath    defining the save path for the contatenated file (ops.fbinary)
ampX = channel_dat{1};
ampX = ampX(5);
fname = strcat(exp_id,'_amp',ampX,'.dat');  % contatenated file is named after this rule, for example 18-138_LPrmGtACR_ampA.dat
tempfname = strcat(exp_id,'_amp',ampX,'_temp.dat');
ops.fbinary = fullfile(savePath,exp_id,fname);  % the kilosort input file; in fullfile(savePath,exp_id); this is the full path including file name
savePath = fullfile(savePath,exp_id);           % where the contatenated file is saved; this is the directory path

% Make sure fnamestomerge is sorted in order
chNum=[];
for i=1:size(channel_dat,2)
    chName = channel_dat{i};
    chNum(end+1) = str2num(chName(7:9));
end
[~,~,ic]=unique(chNum);
channel_dat = channel_dat(ic);

numCh = length(channel_dat);

% Check to see if the same file already exists
if exist(ops.fbinary,'file')
    % if yes, ask to overwrite
    prompt = 'Merged file already exists. Overwrite? y/n: ';
    str = input(prompt,'s');
    if strcmp(str,'y')
        system('SET COPYCMD=/Y')
        mergeFlag = true;
    else
        mergeFlag = false; % exit code
    end
else
    mergeFlag = true;
end

%% Select probe configuration ****guide me to select
[amp, okstatus]=electrodeSelectorSingle;

if okstatus
    if ~strcmp(amp, 'None')
        switch amp
            case 'A1x16'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x16.mat'; % need to be changed
                ops.NchanTOT = 16;
                ops.Nchan = 16;
            case 'A1x16-Poly2'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x16Poly2.mat'; % need to be changed REMOVED POLY2 AT THE END OF PATH BECAUSE FOR SOME REASON MY OWN ADDED CHAN MAP DOES NOT WANT TO BE READ IN THE GUI!
                ops.NchanTOT = 16;
                ops.Nchan = 16;
            case 'A1x32-Poly2'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Poly2.mat'; % need to be changed
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A1x32-Edge'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Edge.mat'; % need to be changed
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A2x32'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA2x32.mat'; % need to be changed
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A2x32-Poly5'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA2x32-Poly5.mat'; % need to be changed
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A1x64-Poly2'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x64Poly2.mat'; % need to be changed
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A1x32'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Classic.mat'; % need to be changed
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A1x32-Poly3'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Poly3.mat'; % need to be changed
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A2x16-Poly2'
                ops.chanMap = 'E:\extarcellular\Andreanne\chanMaps\chanMapA2x16-Poly2.mat'; % need to be changed
                ops.NchanTOT = 32;
                ops.Nchan = 32;
        end
        ops.Nfilt = ops.NchanTOT*4; %ops.Nfilt = ops.Nchan*2;
    else
        error('No electrode configuration selected')
    end
else
    error('Process cancelled by user')
end

if numCh~=ops.NchanTOT
    error('Probe configuration does not match the number of files selected')
end

%%
if mergeFlag
    % merge files and save
    mkdir(fullfile(savePath,exp_id));
    if ispc  %if it is pc
        cmdline = 'copy /B ';
        if numCh>1
            for i=1:numCh-1
                chPath = fullfile(recording_dir,channel_dat{i}); % Path to the individual files
                cmdline = [cmdline chPath ' + '];
            end
            chPath = fullfile(recording_dir,channel_dat{end}); % Path to the individual files
            cmdline = [cmdline chPath ' ' ops.fbinary];
        end
    elseif isunix || ismac
        if numCh>1
            cmdline = 'cat ';
            for i=1:numCh-1
                chPath = fullfile(recording_dir,channel_dat{i}); % Path to the individual files
                cmdline = [cmdline chPath ' '];
            end
            chPath = fullfile(recording_dir,channel_dat{end}); % Path to the individual files
            cmdline = [cmdline chPath ' > ' ops.fbinary];
        else
            cmdline = ['mv ' channel_dat{1} ' ' ops.fbinary];
        end
    else
        error('Cannot determine the OS')
    end
    
    status = system(cmdline); % call Win command prompt to concatenate files; it is faster than reading file into the memory than writing back to harddrive 
    
    if status == 0   
        % transpose the data % becuase the data need to be read in blocks
        % including all channels.
        fileinfo = dir(ops.fbinary);
        nSamp = fileinfo.bytes/2; % int16 has two bytes (16bits)
        mmf = memmapfile(ops.fbinary, 'Format', {'int16', [nSamp/numCh numCh], 'x'},'Writable',true); %memory remaping on harddrive; much more efficient
        mmf.Data.x = mmf.Data.x';
    end
    % data is now channels x time; [numCh sampleperCh]
    fprintf('File merging complete \n');

end
