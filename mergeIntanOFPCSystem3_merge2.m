function [ops,savePath] = mergeIntanOFPCSystem3_merge2(savePath, intanpath)
% the folder containing intan data files should be named after XXXX_ephys; works for one file per channel intan aquisition (OFPC)
% There cannot be any space in the path, otherwise error message "The system cannot find the file specified" will show
% Select files to merge in RECORDING ORDER
if nargin < 2
    intanpath = 'E:\extarcellular\Andreanne\IntanData';
end

tmp =strfind(savePath,'\');
chanMap_parent = savePath(1:(tmp(end)-1));
batch_size = 8; % can take care of four recording session for batch size 16; reduce to 8 to take care of 8 batches

%grab multiple folders
[filenamescell] = grabMultipleRecordings('*amp-A*.dat', intanpath);
num_folders = length(filenamescell);
num_files = length(filenamescell{1});
fpathtomerge = cell(num_folders, 1);
fnamestomerge = cell(num_files, 1);

for s=1:num_folders
    fpathtomerge{s} = filenamescell{s}(1).folder;
end

for t=1:num_files %all channel numbers must be identical between recordings
    fnamestomerge{t} = filenamescell{1}(t).name;
end
fnamestomerge=cellstr(fnamestomerge);


%Data folders in IntanData should be named as "cage number - mouse number -
%recording location - stimulation type _ ephys _ XXXXX(date & time)"
recordings = [];
for m = 1:num_folders 
    fullPath = fpathtomerge{m};
    firstchar = strfind(fullPath,'\');% find "\"
    firstchar = firstchar(end)+1;% the letter after last "\" 
    mergechar = strfind(fullPath,'-');% find "-"
    mergechar = mergechar(end)+1;% the letter after last "-"
    lastchar=strfind(fullPath,'_ephys')-1;% the letter before "_ephys"
    if m ==1
        recordings = [recordings fullPath(mergechar:lastchar)];
    else
        recordings = [recordings 'AND' fullPath(mergechar:lastchar)];
    end
end
% the information between the first and the merge char should be the same
% for all folders, so just use the last one
exp_id = [fullPath(firstchar:mergechar-1) recordings];

% fPath    defining the save path for the contatenated file (ops.fbinary)
ampX = fnamestomerge{1};
ampX = ampX(5);
fname = strcat(exp_id,'_amp',ampX,'.dat');  % contatenated file is named after this rule, for example 18-138_LPrmGtACR_ampA.dat
tempfname = strcat(exp_id,'_amp',ampX,'_temp.dat');
ops.fbinary = fullfile(savePath,exp_id,fname);  % the kilosort input file; in fullfile(savePath,exp_id); this is the full path including file name
savePath = fullfile(savePath,exp_id);           % where the contatenated file is saved; this is the directory path

% Make sure fnamestomerge is sorted in order
chNum=[];
for i=1:num_files
    chName = fnamestomerge{i};
    chNum(end+1) = str2num(chName(7:9));
end
[~,~,ic]=unique(chNum);
fnamestomerge = fnamestomerge(ic);

numCh = length(fnamestomerge);

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
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x16.mat']; 
                ops.NchanTOT = 16;
                ops.Nchan = 16;
            case 'A1x16-Poly2'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x16Poly2.mat']; % need to be changed REMOVED POLY2 AT THE END OF PATH BECAUSE FOR SOME REASON MY OWN ADDED CHAN MAP DOES NOT WANT TO BE READ IN THE GUI!
                ops.NchanTOT = 16;
                ops.Nchan = 16;
            case 'A1x32-Poly2'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x32Poly2.mat']; 
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A1x32-Edge'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x32Edge.mat'];
                ops.NchanTOT = 32;
                ops.Nchan = 32;
            case 'A2x32'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA2x32.mat']; 
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A2x32-Poly5'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA2x32-Poly5.mat']; 
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A1x64-Poly2'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x64Poly2.mat']; 
                ops.NchanTOT = 64;
                ops.Nchan = 64;
            case 'A1x32'
                ops.chanMap = [chanMap_parent '\chanMaps\chanMapA1x32Classic.mat']; 
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
ite = floor(numCh/batch_size);
remainder = rem(numCh, batch_size);
tempathcell={};

%this really long and complicated loop was created because command prompt
%required to merge these files cannot handle a very long command. So the
%mergin is broken up into batches

if mergeFlag
    mkdir(savePath); %create the output path file
    if ispc  %if it is pc     
        if numCh>1  
            for p = 1:ite %for each batch
                cmdline = 'copy /B ';
                tempathcell{p} = fullfile(savePath,['temp' num2str(p) '.dat']); % temp destination
                for i=1:batch_size %all, but not the last channel
                    for j=1:num_folders
                        chPath = fullfile(fpathtomerge{j},fnamestomerge{(p-1)*batch_size+i}); % Path to the individual files
                        if i == batch_size && j == num_folders
                            cmdline = [cmdline chPath ' ' tempathcell{p}];
                        else
                            cmdline = [cmdline chPath ' + '];
                        end
                    end 
                end
                status = system(cmdline); % call Win command prompt to concatenate files; it is faster than reading file into the memory than writing back to harddrive 
            end 
            if remainder
                tempathcell{p+1} = fullfile(savePath,['temp' num2str(p+1) '.dat']); 
                cmdline = 'copy /B ';
                for i=1:remainder %catch if number of channels is not factor of batchsize
                     for j=1:num_folders
                        chPath = fullfile(fpathtomerge{j},fnamestomerge{p*batch_size+i}); % Path to the individual files
                        if i == remainder && j == num_folders
                            cmdline = [cmdline chPath ' ' tempathcell{p+1}];
                        else
                            cmdline = [cmdline chPath ' + '];
                        end
                    end 
                end
                status = system(cmdline);
            end
            % call Win command prompt to concatenate files; it is faster than reading file into the memory than writing back to harddrive
            %last file must close the merging process (ie add the last chan
            %for each recording, but not the very last file to merge
            cmdline = 'copy /B ';
            for i=1:length(tempathcell)
                if i == length(tempathcell)
                    cmdline = [cmdline tempathcell{i} ' ' ops.fbinary];
                else
                    cmdline = [cmdline tempathcell{i} ' + '];
                end               
            end
            status = system(cmdline); % call Win command prompt to concatenate files; it is faster than reading file into the memory than writing back to harddrive 
        end
        
    elseif isunix || ismac
        if numCh>1
            cmdline = 'cat ';
            for i=1:numCh-1
                chPath = fullfile(fpathtomerge,fnamestomerge{i}); % Path to the individual files
                cmdline = [cmdline chPath ' '];
            end
            chPath = fullfile(fpathtomerge,fnamestomerge{end}); % Path to the individual files
            cmdline = [cmdline chPath ' > ' ops.fbinary];
        else
            cmdline = ['mv ' fnamestomerge{1} ' ' ops.fbinary];
        end
        status = system(cmdline); % call Win command prompt to concatenate files; it is faster than reading file into the memory than writing back to harddrive 
    else
        error('Cannot determine the OS')
    end
    
    
    if status == 0   
        % transpose the data % becuase the data need to be read in blocks
        % including all channels.
        fileinfo = dir(ops.fbinary);
        nSamp = fileinfo.bytes/2; % int16 has two bytes (16bits)
        mmf = memmapfile(ops.fbinary, 'Format', {'int16', [nSamp/numCh numCh], 'x'},'Writable',true); %memory remaping on harddrive; much more efficient
        mmf.Data.x = mmf.Data.x';
        
        fprintf('Deleting temp files \n');
        %deleted temp files
        for v=1:length(tempathcell)
            delete(tempathcell{v})
        end
        
        % data is now channels x time; [numCh sampleperCh]
        fprintf('File merging complete \n');
    else
        fprintf('File merging incomplete. \n');
    end
    
    %% merge analog inputs : includes trigger time, stim time, led, cam, 
%     analogfilespathcell = grabMultipleRecordings('board-ADC*.dat', intanpath);
    filetype='board-ADC*.dat';
    cd(filenamescell{1}(1).folder); % go into a selected folder
    analogfilespath_temp=dir(filetype);  % search of board-ADC files
    %concatilate board-adc files
    for j=1:length(analogfilespath_temp) % number of analog channel
        cmdline = 'copy /B ';
        destinationpath = fullfile(savePath, analogfilespath_temp(j).name);
        num_folders=length(filenamescell);
        for i=1:num_folders % folder/recording session
            cd(filenamescell{i}(1).folder);
            analogfilespath=dir(filetype);
            chPath = fullfile(analogfilespath(j).folder, analogfilespath(j).name);
            if i == num_folders
                cmdline = [cmdline chPath ' ' destinationpath];
            else
                cmdline = [cmdline chPath ' + '];
            end               
        end
        status = system(cmdline);
    end
    
    %% merge time.dat
    destinationpath = fullfile(savePath, 'time.dat');
    fid_dest = fopen(destinationpath, 'w'); 
    num_folders=length(filenamescell);
    lastnum=[];
    firstnum=[];
    for i=1:num_folders % folder/recording session
        timeFullFileName = fullfile(filenamescell{i}(1).folder, 'time.dat'); 
        fileinfo = dir(timeFullFileName);
        num_samples = fileinfo.bytes/4; % int32 = 4 bytes 
        fid = fopen(timeFullFileName, 'r'); 
        timestamps = fread(fid, num_samples, 'int32'); 
        fclose(fid);        
        if i == 1
            lastnum=timestamps(end);
        else
           firstnum=timestamps(1);
           timestamps=timestamps-(firstnum-lastnum)+1;
           lastnum=timestamps(end);
        end        
        fwrite(fid_dest,timestamps,'int32');
    end
    fclose(fid_dest);
 
    %% copy the info.rhd to kilosort folder 
    status = system(['copy /B ' fullfile(filenamescell{1}(1).folder, 'info.rhd') ' ' fullfile(savePath, 'info.rhd')]);
end
