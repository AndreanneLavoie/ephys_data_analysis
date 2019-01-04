function spikes2 = create_spikes()
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
read_Intan_RHD2000_file

spikes.param.Fs = frequency_parameters.amplifier_sample_rate;
%spikes.param.detect_method = 
%spikes.param.thresh =
%spikes.param.window_size = 
%spikes.param.cross_time =
%spikes.param.refractory_period =
%spikes.param.max_jitter =
%spikes.param.agg_cuttoff =
%spikes.param.kmeans_clustersize =
%spikes.param.display =
%spikes.param.initial_split_figure_panel = 

%spikes.info.

spikes.spiketimes = readNPY('spike_times.npy')
spike.waveforms = 

spikes.vstiming = vs_time2('time.dat', 'board-ADC-00', path);
spikes.vstrigger = vs_time2('time.dat', 'board-ADC-01', path);

end


