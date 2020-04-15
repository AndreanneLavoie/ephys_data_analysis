function runKiloSort4Intan2

% For other users: change savePath, make new StandardConfig

%% Paths

addpath(genpath('C:\Users\admin-liubaohu\Documents\KiloSort-master')) % path to kilosort folder, in case it's not on path
addpath(genpath('C:\Users\admin-liubaohu\Documents\npy-matlab-master')) % path to npy-matlab scripts downloaded from GitHub
savePath = 'E:\extarcellular\Andreanne\KilosortData'; % an experiment-specific subdirectory will be created under this, where kilosort input data files and sorting result files go
intanpath='E:\extarcellular\Andreanne\IntanData';     % an experiment-specific subdirectory will be created under this
cd(intanpath);
%% Select probe configuration and merge Intan one file per channel dat files into one
[ops,savePath] = mergeIntanOFPCSystem3(savePath); %%% the input savePath is the parent folder where all kilosort data folders (experiment specific) are
                                                  %%% the output savePath is the path fof the experiment specific folder (automatically generated)
%% Load the configuration file
run('E:\extarcellular\Andreanne\code\kilosort code\StandardConfig_SKM.m')

%%
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