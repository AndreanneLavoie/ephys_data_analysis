function [refractory_period_violation, clusters]= checkRefractoryPeriodViolations_KiloSort_AL(spikes, assigns, tasks)
%INPUT
%   -spikes construct
%   -assigns (cluster IDs of interest in a vector); use [] for 'all' case
%   -task: 'all' check rp for all clusterIDs, 'merge' to check the 
% pass a vector of cluster IDs to check. if empty, it reports violations for all good
% units

%% Parameters/data
% Violation threshold in ms
thrVio = 2;%1.5;
thrRefViol = 0.2; %percent
% Convert spiketimes to ms
spike_times = (spikes.abs_spiketimes)*1e3;

% Read spike_clusters.npy
spike_clusters = double(readNPY([spikes.kilosort_path 'spike_clusters.npy']));
clusters=spike_clusters';

%convert tsv to csv (as we cannot readtable does not read tsv files)
% if ~isfile([spikes.kilosort_path 'cluster_groups.csv'])
%     copyfile([spikes.kilosort_path 'cluster_group.tsv'], [spikes.kilosort_path 'cluster_groups.csv']);
% end
clusterID = unique(readNPY([spikes.kilosort_path 'spike_clusters.npy'])'); %identity of cluster 
%[clusterID, labels] = readClusterGroupsCSV([spikes.kilosort_path 'cluster_groups.csv']);
%initialize variables

refractory_period_violation = nan(length(clusterID), 5); 
refractory_period_violation(:, 1) = clusterID;

%% Calculate violations
switch tasks
    case 'all'
        
        for i=1:length(clusterID)
            st = spike_times(spike_clusters == clusterID(i)); %extract all the spiketimes that are associated with given unit
            inspk = diff(st); %calculates the interspike intervals 
            numspk = length(st); %determines how many intervals are smaller than refractory period (1ms)
            numVio = sum(inspk<thrVio); %add up the number of violations
            percent_vio =  numVio/numspk*100; %percent of violations
            consider_bad = percent_vio > thrRefViol;
            refractory_period_violation(i,2:5) = [numVio, numspk, percent_vio, consider_bad];

            %fprintf('cluster %d, %d violations in %d spikes \n',goodID(i),numVio,numspk)
        end

    case 'select'
        
        refractory_period_violation = nan(length(assigns),5); 
        
        for i=1:length(assigns)
            if any(ismember(assigns(i),spike_clusters))
                st = spike_times(spike_clusters == assigns(i));
                inspk = diff(st);
                numspk = length(st);
                numVio = sum(inspk<thrVio);
                percent_vio =  numVio/numspk*100; %percent of violations
                consider_bad = percent_vio > thrRefViol;
                refractory_period_violation(i,:) = [assigns(i), numVio, numspk, percent_vio, consider_bad];
                %fprintf('cluster %d, %d violations in %d spikes \n',spkIDs(i),numVio,numspk)
            else
                fprintf('cluster %d does not exist \n',assigns(i));
            end
        end
    case 'merge'
        
%         spikeflag=zeros(size(spike_times));
%         for i=1:length(assigns)
%             spikeflag = spikeflag | (spike_clusters == assigns(i));
%         end
        spike_index = [];
        for i=1:length(assigns) %find the index values of spiketimes for each clusterID and combine them together in a single list
            spike_index = [spike_index; find(spike_clusters == assigns(i))];
        end
        sorted_spike_index = sort(spike_index); %sort the list in order
        st = spike_times(sorted_spike_index); %spikeflag used to be used as index (logic values)
        inspk = diff(st);
        numspk = length(st);
        numVio = sum(inspk<thrVio);
        percent_vio =  numVio/numspk*100; %percent of violations
        consider_bad = percent_vio > thrRefViol;
        refractory_period_violation = [-1, numVio, numspk, percent_vio, consider_bad];    
end
        
    