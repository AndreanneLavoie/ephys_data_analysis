tuning = 'sf';

%% get compiled data
datapath='E:\extarcellular\Andreanne\KilosortData';%parent folder of all the data

data_table = readtable('E:\extarcellular\Andreanne\data_analysis\phototagging_recordings_summary.xlsx');

switch tuning
    case 'sf'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 1);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 1, :);
    case 'tf'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 3);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 3, :);
    case 'ori'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 2);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 2, :);
end

%get spikes.mat path
data2compile={};
for i = 1: length(recordinglist)
        data2compile{end+1} = [datapath '\' recordinglist{i} '\' 'spikes.mat']; 
end


%% loop



for i=1:length(data2compile)
    load(data2compile{i});%load spike construct to workspace
    
        psth_depth_any_probe_liu(spikes, tuning , 'chanMapA2x32-Poly5.mat');
        savefig(['E:\extarcellular\Andreanne\data_analysis\photo-tag\depth recording analyses\' tuning '\' recordinglist{i}]);

end
    