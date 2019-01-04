function spikes = create_spikes()
%this code will create a structure called spikes which will hold all of the
%data required for analysis 

%data from different sources should be collected together in a folder located in filedir; 
%includes: 
%from intan 'time.dat', 'board-ADC-00.dat' (00-04),
%(record of analog inputs for different stim; led, eye track, vs timing, behaviour etc.)
%'info.rhd'; 
%from psychopy: order of stimuli '..._vs_date.npy',
%from kilosort/phy: 'spike_times.npy', 'spike_clusters.npy' and
%'cluster_groups.tsv

%must have .npy to matlab installed from https://github.com/kwikteam/npy-matlab

%location data stored: CHANGE
read_Intan_RHD2000_file %the chosen folder path wil have the same exp name as vs >> use this instead

stim_order_file = 'test2_vs_ori_20181214_105526.npy';

spikes.times = readNPY(fullfile(path, 'spike_times.npy')); %absolute time point for each spike   
spikes.stimorder = readNPY(fullfile(path,'stim_order_file'));

spikes.vstiming = vs_time2('time.dat', 'board-ADC-00', path);
spikes.trigger = vs_time2('time.dat', 'board-ADC-01', path); %visual stim trigger times

%Both spiketimes, trials and stimcond are genereated together; 
%spiketimes: relative time of spike to trigger;
%stimcond: the type of stimulus associated with each spike, given the closest trigger vsstim
%trigger identity
[spikes.spiketimes, spikes.trials, spikes.stimcond] = make_rel_spiketimes(spikes.times, spikes.trigger, spikes.stimorder);

spikes.label = tdfread(fullfile(path,'cluster_group.tsv')); %identity of cluster (ie good, MUA (multi-unit activity) or noise)
spikes.assign = readNPY(fullfile(path, 'spike_clusters.npy'));

%spikes.fileID = if multiple files are generated;

end


