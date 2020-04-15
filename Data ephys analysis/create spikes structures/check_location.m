%%check for correct NOT/DTN recording location & population led activity
% i.e. look at population psth (and tuning - for ori only)
%will automatically save figures in correct folder with proper names

%must load spikes construct manually 

%filter out noise
assigns = spikes.labels(spikes.labels(:,2) ~= 0);
good_spikes = filtspikes(spikes, 0, 'assigns', assigns);

%psth
figure;
h = psth_bl(good_spikes)
savefig([spikes.kilosort_path 'pop_psth(noise)']) 

%tuning for 'ori'
if ~isempty(strfind(spikes.kilosort_path, 'ori'))
    figure;
    l = compute_tuning(spikes, assigns, 'ori', 1, 1)
    savefig([spikes.kilosort_path 'pop_tuning']) 
end

 
