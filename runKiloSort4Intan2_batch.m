function runKiloSort4Intan2_batch(user,selection_mode)
%%user
%the initial of the code user: AL, HH
%this will determined the Intan path, Kilosort path, and the electrode
%configurations (if two people are using the same number of channels but
%different electrode configuration)

%%selection_mode: 
%'manual': a selection window will pop up for manually select folders
%'auto': will automatically find un-analyzed folders start with "vgt-" and "c57-"
%       to add or remove the mouse strain, edit inputs for function "startsWith" under "Paths" section 
%       the string needs to be under "" for startsWith to recognize as patterns

%%the function will generate a pop-up window to display error messages for
%%each folders lacking amp dat files

%%if multiple intan folders exist a specific mouse at a specific recording
%%loaction, the function will run on the larges recording

%%to update proconfiguration type edit in "recording-wise analysis"
%%section, amp value

if nargin<2
    selection_mode='auto'; %auto selection by default
end

switch user
    case 'AL'
        kilosortpath = 'E:\extarcellular\Andreanne\KilosortData'; % an experiment-specific subdirectory will be created under this, where kilosort input data files and sorting result files go
        intanpath='E:\extarcellular\Andreanne\IntanData';     % an experiment-specific subdirectory will be created under this
        electrode_config ='A1x32-Poly3';%'A2x16-Poly2'; %the 32 channel electrode Andreanne is using
    case 'HH'
        kilosortpath = 'E:\extarcellular\Hayley\KilosortData'; % an experiment-specific subdirectory will be created under this, where kilosort input data files and sorting result files go
        intanpath='E:\extarcellular\Hayley\IntanData';     % an experiment-specific subdirectory will be created under this
        electrode_config ='A1x32-Poly3'; %the 32 channel electrode Hayley is using
end
  

%% Paths
%parent paths
addpath(genpath('C:\Users\admin-liubaohu\Documents\KiloSort-master')) % path to kilosort folder, in case it's not on path
addpath(genpath('C:\Users\admin-liubaohu\Documents\npy-matlab-master')) % path to npy-matlab scripts downloaded from GitHub
cd(intanpath);

%paths for individual recordings
switch selection_mode
    case 'manual' %manual selection of Intan folders, exit selection by choosing "cancel"
        directorycell = {};
        while 1 
            directoryname = uigetdir(intanpath, 'Pick IntanData Directory');
            if ~directoryname
                break;
            else
                directorycell{end + 1} = directoryname;
            end
        end
    case 'auto'
        %current NOT-DTN recordings all starts with vgt or c57
        %select the largest recordings if multiple intan folders for the same mice & location exist
        all_intan_folders = dir(intanpath); %information for all folders in intan folder
        list_intan=string({all_intan_folders([all_intan_folders.isdir]).name})';%the names for all folders in intan folder
        list_intan=list_intan(startsWith(list_intan,["vgt-","c57-"])); %find all recordings starts with vgt or c57
        all_kilosort_folders=dir(kilosortpath);%information for all folders in kilosrot folder
        list_kilosort=string({all_kilosort_folders([all_kilosort_folders.isdir]).name})';
        list_kilosort=list_kilosort(startsWith(list_kilosort,["vgt-","c57-"])); %find all recordings starts with vgt or c57
        %check if a recording has been analyzed by kilosort
        intan_to_analyze={};
        %intan_to_analyze:the first colomn saves the intan name and the second column saves the calculated kilosort name
        for i = 1:length(list_intan)
            intan_name = list_intan{i};
            lastchar=strfind(intan_name,'_ephys')-1;
            kilosort_name = intan_name(1:lastchar); %caculate the corresponding killosort name 
            %find all recordings that don't have a corresponding kilosort folder
            if ~any(strcmp(list_kilosort,kilosort_name))
                if isempty(intan_to_analyze)
                    intan_to_analyze = {intan_name kilosort_name};
                elseif ~any(strcmp(intan_to_analyze(:,2),kilosort_name))
                    %the first colomn saves the intan name and the second column saves the calculated kilosort name
                    intan_to_analyze(end+1,:) = {intan_name kilosort_name};
                else %compare the size of the new intan folder to the one on the list, and keep the largest folder
                    %pull out the intan folder name for comparison
                    intan_name_compare = intan_to_analyze{strcmp(intan_to_analyze(:,2),kilosort_name),1};
                    temp_old=dir(intan_name_compare);
                    temp_new=dir(intan_name);
                    %sum up size for the files (not directory)
                    size_old=sum(cell2mat({temp_old(~[temp_old.isdir]).bytes}));
                    size_new=sum(cell2mat({temp_new(~[temp_new.isdir]).bytes}));
                    if size_old<size_new %replace the old folder by the new folder
                       intan_to_analyze(strcmp(intan_to_analyze(:,2),kilosort_name),:)={intan_name, kilosort_name};
                    end                    
                end
            end
        end
        if isempty(intan_to_analyze)
            fprintf('no new intan folders to run\n')
            return
        else
            directorycell=strcat(intanpath, '\', intan_to_analyze(:,1));
        end
end
%% recording-wise analysis
for i = 1:length(directorycell)
    %select chanel configuration based on channel number
    cd(directorycell{i})
    if size(dir('amp-A*.dat'),1)==0
        f=msgbox(['no amp dat exist in ' directorycell{i}]);%generate a pop-up window for error messages
        continue
    elseif size(dir('amp-A*.dat'),1)==64
        amp = 'A2x32-Poly5';
    elseif size(dir('amp-A*.dat'),1)==16
        amp = 'A1x16';
    elseif size(dir('amp-A*.dat'),1)==32
        amp = electrode_config;
    end
    
%% Select probe configuration and merge Intan one file per channel dat files into one
    [ops,savePath] = mergeIntanOFPCSystem3_batch(directorycell{i}, kilosortpath , amp); 
                                                  %%% the input: fpathtomerge is the name of the folder (one for each recording session) which contains the raw data (one file per channel)
                                                  %%%savePath is the parent folder where all kilosort data folders (experiment specific) are
                                                  %%% the output savePath is the path fof the experiment specific folder (automatically generated)
    %% Load the configuration file
    run('E:\extarcellular\Andreanne\code\kilosort code\StandardConfig_SKM.m')

    %% kilosort
    if ops.GPU     
        gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
    end

    [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
    rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
    rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

    % AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
    % rez = merge_posthoc2(rez);

    % save matlab results file
    save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

    % save python results file for Phy
    rezToPhy(rez, ops.root);

    % remove temporary file
    delete(ops.fproc);
    
    %% create_spikes
    cd (savePath)
    try 
        spikes=create_spikes(directorycell{i});
    catch ME
        f=msgbox({'line1';'line2'},directorycell{i}, ME.message)
        continue
    end
    fprintf([num2str(length(directorycell)-i) ' recordings left\n'])
end