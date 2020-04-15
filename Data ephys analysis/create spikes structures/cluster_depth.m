function cluster_depths = cluster_depth(spikes)
%extract cluster channel depth based on average of individual spikes
%amplitudes for each cluster

cluster_depths = zeros(length(spikes.labels),1);

for i=1:length(spikes.labels)

    filt_spikes = filtspikes(spikes, 0, 'assigns', spikes.labels(i,1));
    cluster_depths(i) = mean(filt_spikes.channel_depth);    
end

end

