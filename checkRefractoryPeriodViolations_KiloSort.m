function checkRefractoryPeriodViolations_KiloSort(varargin)

% pass a vector of cluster IDs to check. if empty, it reports violations for all good
% units
%% Violation threshold in ms
thrVio = 1;
fprintf('Inter spike interval violation defined as < %d ms \n',thrVio);
%% Choose the directory containing npy files from KiloSort
expDir = uigetdir('E:\Documents\DATA','Select directory with KiloSort results');

%% Read in parameters
rez = load([expDir '\rez.mat']);
params = rez.rez.ops;

%% Read spike_times.npy
spike_times = double(readNPY([expDir '\spike_times.npy'])); % samples
% Convert to ms using samping rate
spike_times = spike_times/params.fs*1e3;

%% Read spike_clusters.npy
spike_clusters = double(readNPY([expDir '\spike_clusters.npy']));

%% Calculate violations
if ~nargin
    formatSpec = '%f%C';
    T = readtable([expDir '\cluster_groups.csv'],'Format',formatSpec);
    goodID = T.cluster_id(T.group=='good',:);
    for i=1:length(goodID)
        st = spike_times(spike_clusters == goodID(i));
        inspk = diff(st);
        numspk = length(st);
        numVio = sum(inspk<thrVio);
        fprintf('cluster %d, %d violations in %d spikes \n',goodID(i),numVio,numspk)
    end
else
    spkIDs = varargin{1};
    for i=1:length(spkIDs)
        if any(ismember(spkIDs(i),spike_clusters))
            st = spike_times(spike_clusters == spkIDs(i));
            inspk = diff(st);
            numspk = length(st);
            numVio = sum(inspk<thrVio);
            fprintf('cluster %d, %d violations in %d spikes \n',spkIDs(i),numVio,numspk)
        else
            fprintf('cluster %d does not exist \n',spkIDs(i));
        end
    end
end
        
    