function [spike_times,spike_clusters, spike_groups] = extract_spikes(sp_fileName, sc_fileName, sg_fileName ,filedir)
%This function will extract .npy files from kilosort/phy
%inlcuding spike_times, spike_cluster 
%it will also extract the good clusters from cluster_groups.tsv

%NOTE: make sure to have npy to matlab setup from 

sp_FullfileName = fullfile(filedir, sp_fileName);
sc_FullfileName = fullfile(filedir, sc_fileName);
sg_FullfileName = fullfile(filedir, sg_fileName);

spike_times = readNPY(sp_FullfileName);
spike_clusters = readNPY(sc_FullfileName);
spike_groups = readNPY(sg_FullfileName);

end

