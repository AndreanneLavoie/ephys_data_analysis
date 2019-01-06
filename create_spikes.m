function spikes = create_spikes()
%this code will create a structure called spikes which will hold all of the
%data required for data analysis 

%data from different sources should be collected together in a folder located in variable filedir;
%filedir will be automatically generated when you select this folder
%(through GUI) when running this code;

%For the code to work, this SECLECTED DATA FOLDER MUST INCLUDE: 
%from intan:
%'time.dat', 'board-ADC-00.dat' (00-04),
%(record of analog inputs for different stim; led, eye track, vs timing, behaviour etc.)
% and 'info.rhd'; 
%from psychopy: order of stimuli '..._vs_date.csv',
%from kilosort/phy: 'spike_times.npy', 'spike_clusters.npy' and
%'cluster_groups.tsv

%BUG: csv file does not want to automatically be selected; works in comand line; 
%TEMP SOLUTION; manually change name of vs sitm order file; currently
%marked as %CHANGE (line 31)

%must have .npy to matlab installed from https://github.com/kwikteam/npy-matlab

%%
%location data stored: 
read_Intan_RHD2000_file %the chosen folder path will be extracted in this line when select info.rhd file

%name of vs stim order file 
%get_stim_order = dir(fullfile(path, '*.csv')); %automatically select the only .csv file in path folder DOES NOT WORK YET
%stim_order_file = get_stim_order.name;
stim_order_file = 'savetest3_vs_loc_20190105_180311.csv' %CHANGE

spikes.times = readNPY(fullfile(path, 'spike_times.npy')); %absolute time point for each spike   
spikes.stimorder = csvread(fullfile(path, stim_order_file));

spikes.vstiming = vs_time('time.dat', 'board-ADC-00', path);
spikes.trigger = vs_time('time.dat', 'board-ADC-01', path); %visual stim trigger times

%Both spiketimes, trials and stimcond are genereated together; 
%spiketimes: relative time of spike to trigger;
%stimcond: the type of stimulus associated with each spike, given the closest trigger vsstim
%trigger identity
[spikes.spiketimes, spikes.trials, spikes.stimcond] = make_rel_spiketimes(spikes.times, spikes.trigger, spikes.stimorder);

spikes.label = tdfread(fullfile(path,'cluster_group.tsv')); %identity of cluster (ie good, MUA (multi-unit activity) or noise)
spikes.assign = readNPY(fullfile(path, 'spike_clusters.npy'));

%spikes.fileID = if multiple files are generated;

end


