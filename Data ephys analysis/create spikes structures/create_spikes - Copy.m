function spikes = create_spikes(intan_data_path)
%this code will create a structure called spikes which will hold all of the
%data required for data analysis after Phy (manual spike sorting)

%data from different sources should be collected together in a folder three seperate folders, 
%organised by location/software used to generate the data; 

    % This includes: 
% -> From IntanData -- folder containing raw ephys data, and vs/led timing ('info.rhd', 'time.dat', 'board-ADC-00.dat' (00-04)(record of analog inputs for different stim; led, eye track, vs timing, behaviour etc.)
% -> From Labview -- order of stimuli ( '20trk1DATE_expName.eye' ) or for ori tuning .csv file form stim vis computer,
% -> From Kilosort -- 'spike_times.npy', 'spike_clusters.npy' ETC.
% -> From Phy --  'cluster_group.tsv'

%must have .npy to matlab installed from https://github.com/kwikteam/npy-matlab

%%%

%location of Intan data stored: 
%[intandata_path, frequency_parameters] = read_Intan_RHD2000_file; %the chosen folder path wil have the same exp name as vs
intandata_path = intan_data_path;

%% extract experiment name and extrapolate path for different types of data using the Intan info.rhd file selected;

%find index of key markers in the intandata filename
gen_path_ind = strfind(intandata_path, 'IntanData'); %28 (index of 'I' in IntaData) 
exp_ind = strfind(intandata_path, '_ephys'); %index of '_' in _ephys 

%using indexes and intandata filename to assign path of the other data
%folders (visual stimulus) and set new path to kilosort file where analysed
%data will be saved
general_path = intandata_path(1:gen_path_ind - 1); %'E:\extarcellular\Andreanne\'
exp_name = intandata_path(gen_path_ind + 10 : exp_ind - 1); %everything between 'IntanData\' and '_ephys' 
kilosortdata_path = [general_path 'KilosortData\' exp_name '\'];
date = intandata_path((exp_ind + 6):(exp_ind + 12)); %e.g. '_190223'
vspara_path = [general_path 'VSpara\'];

tuning_param_ind = strfind(exp_name,'-');
tuning_param_ind = tuning_param_ind(end)+1;
tuning = exp_name(tuning_param_ind:end);

%grab proper vsparam; if comes from eye tracking computer .eye format; if
%comes from vis stim computer uses .csv format (only for orientation tuning)
if ~strcmp(tuning(1:2), 'or')
    vspara_filename = ['20trk' date(6:7) date(4:5) date(2:3) '_' exp_name '.eye'];
%load the vs stimulation parameters
    [VarParam spikes.labview_param] = getPsychStimParameters_BL(vspara_path, vspara_filename);
    spikes.vs_params = spikes.labview_param.stimpara';
    spikes.vs_params = spikes.vs_params(1:spikes.labview_param.stimnum, :);
else
    vspara_filename = erase(exp_name,'-');
    spikes.vs_params = csvread([vspara_path vspara_filename]);
end

spikes.raw_data_path = intandata_path;
spikes.kilosort_path = kilosortdata_path;
 
%few required static parameters:
%info_path=[intandata_path 'info.rhd'];
%frequency_parameters =read_Intan_info(info_path);
spikes.frequency_parameters.amplifier_sample_rate= 30000; %frequency_parameters.amplifier_sample_rate;
threshold=1.8;%used to detect rising phase in vs_time2 function
digfilter_winT=0.0006; %in sec


%read kilosort and intan data files into spikes structures
spikes.abs_spiketimes = double(readNPY([kilosortdata_path 'spike_times.npy']))/double(spikes.frequency_parameters.amplifier_sample_rate); %absolute time point for each spike   
spikes.vstiming = extract_signal_time('time.dat', 'board-ADC-00.dat', intandata_path, threshold,digfilter_winT, spikes.frequency_parameters.amplifier_sample_rate);
spikes.trigtiming= extract_signal_time('time.dat', 'board-ADC-01.dat', intandata_path, threshold,digfilter_winT, spikes.frequency_parameters.amplifier_sample_rate); %visual stim trigger times

%% Both spiketimes, trials and stimcond are genereated together; 
%spiketimes: relative time of spike to trigger;
%stimcond: the type of stimulus associated with each spike, given the closest trigger vsstim
%stimorder: trigger identity
spikes = make_spike_struct(spikes);

%for amplitude/depth extracting function:
templates = readNPY([kilosortdata_path 'templates.npy']);
whitening_mat = readNPY([kilosortdata_path 'whitening_mat.npy']);
channel_coords = double(readNPY([kilosortdata_path 'channel_positions.npy']));
ycoords = channel_coords(:,2);
spike_templates = readNPY([kilosortdata_path 'spike_templates.npy']);
tempScalingAmps = readNPY([kilosortdata_path 'amplitudes.npy']);

%extract amplitude and depth for each spikes
[amplitudes, channel_depth, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(templates, whitening_mat, ycoords, spike_templates, tempScalingAmps);

%change orientation
spikes.waveform = waveforms;
spikes.amplitudes = amplitudes';
spikes.channel_depth = round(channel_depth)';

%cluster / unit identity of each spikes
spikes.assigns = readNPY([kilosortdata_path 'spike_clusters.npy'])';

%% import LED info
if sum(spikes.vs_params(:,17)) > 0
    spikes.led = spikes.trials; 
    spikes.led(spikes.trials > 0) = spikes.vs_params(spikes.trials(spikes.trials > 0), 17); %on = 1; off = 0;
    
    %spikes.led_ID(spikes.trials > 0) = led_ID(spikes.trials(spikes.trials > 0));
    spikes.led_abstime = extract_signal_time('time.dat', 'board-ADC-02.dat', intandata_path, threshold, 0); %visual stim trigger times
    %[length_led, spikes.led_spiketimes, led_hist] = make_first_led_spiketimes(spikes, [], 1, [025], 1, 0);
   
    %%%REPLACE INFO BELOW 
    %spikes.spiketimes_led = rem(spikes.spiketimes - spikes.vs_params(1, 11)/1000, spikes.vs_params(1, 13)/1000); %spiketimes relative to led time (using delay and period as bin)
    %spikes.led_ID = floor((spikes.spiketimes - spikes.vs_params(1, 11)/1000) / (spikes.vs_params(1, 13)/1000));
    %spikesflag=spikes.spiketimes<spikes.vs_params(1, 11)/1000 | spikes.spiketimes >= (spikes.vs_params(1, 11)/1000+ spikes.vs_params(1, 13)/1000*spikes.vs_params(1, 16));
    %spikes.spiketimes_led(spikesflag)=NaN;
    %spikes.led_ID(spikesflag)=NaN;
end

%% Extract spikes.labels field that marks cluster as either: good, mua, noise 
% also create or modify cluter_group.tsv file so that ref_viol clusters are marked as 'mua'

%NOTE: After opening sorting data in Phy -> automatically create
%cluster_groups.csv file with labels (which includes unlabelled clusters)

%in constrast, Phy2 creates cluster_group.tsv which DOES NOT store
%unlabelled clusters therefore, we must get all clusterID from spikes.assigns

%get all clusterIDs and label them all as 3 (unsorted)
spikes.labels = unique(spikes.assigns)';
spikes.labels(:,2) = 3;

%Label units with refractory violoation as mua
spikes.ref_period_violation = checkRefractoryPeriodViolations_KiloSort_AL(spikes, [], 'all');
refviolflag=spikes.ref_period_violation(:,5); %this will be used later below

%read in tsv file and the mannually labeled units
%replace previous labels if labeled change
%convert tsv to csv (as readtable does not read tsv files)
if isfile([spikes.kilosort_path 'cluster_group.tsv'])
    if isfile([spikes.kilosort_path 'cluster_groups.csv'])
       delete([spikes.kilosort_path 'cluster_groups.csv']); %must deleted existing .csv file (because came from old spikes construct
    end
    %create new .csv file that matches info from most recently saved .tsv
    %file (updated in Phy2)
    copyfile([spikes.kilosort_path 'cluster_group.tsv'], [spikes.kilosort_path 'cluster_groups.csv']);
    
end

formatSpec = '%s%s';
csvpath=[spikes.kilosort_path 'cluster_groups.csv'];
T = readtable(csvpath,'Format',formatSpec);
noiseflag=strcmp(T{:,2},'noise');
goodflag = strcmp(T{:,2},'good');
muaflag = strcmp(T{:,2},'mua');
spikes.labels(find(goodflag),2) = 2; %label good before refviolflag to remove any dirty units previously labeled as good
spikes.labels(find(refviolflag),2) = 1;
spikes.labels(find(noiseflag),2) = 0;
spikes.labels(find(muaflag),2) = 1;

%generate new tsv file (for phy2) with only the sorted ones (noise,mua,good)
label_temp = spikes.labels(spikes.labels(:,2) ~= 3, :);

%convert label_temp into table format which can then be saved as .csv file
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
cluster_id = num2cell(label_temp(:,1)); %label_temp in cell format (cluster_id is the name of the 1st colum in .tsv file)
group = cell(length(label_temp)); %label_temp in cell format (group is the name of the 2nd column in .tsv)
for i=1:length(label_temp)
   if label_temp(i,2) == 2
        group(i) = {'good'};
   elseif label_temp(i,2) == 1
        group(i) = {'mua'};
   elseif label_temp(i,2) == 0
        group(i) = {'noise'};
   end
end
%create table and save as .csv
label_temp_table = table(cluster_id, group);
writetable(label_temp_table, csvpath, 'Delimiter', 'tab');

%convert .csv to .tsv to where changes can be seen in phy2
if isfile([spikes.kilosort_path 'cluster_groups.csv'])
    delete([spikes.kilosort_path 'cluster_group.tsv']);
end

copyfile([spikes.kilosort_path 'cluster_groups.csv'], [spikes.kilosort_path 'cluster_group.tsv']);



% if isfile([kilosortdata_path 'cluster_group.tsv'])
% 
%     
%     [clusterID, labels] = readClusterGroupsCSV([kilosortdata_path 'cluster_group.tsv']); %identity of cluster 
%     spikes.labels = [clusterID' labels'];
%     %cluster_depths = cluster_depth(spikes); %generate cluster depth too
%     %spikes.labels(:,3) = round(cluster_depths);
%     spikes.cluster_annotation = cell(6, 2);
%     annotation_list = {'noise' 'mua' 'good' 'unsorted' 'fast' 'reg'};
% 
%     for i = 1:size(spikes.cluster_annotation, 1)
%         spikes.cluster_annotation{i,1} = i - 1;
%         spikes.cluster_annotation{i,2} = annotation_list{i};
%     end
% 
% 
% end

%Summary of visual stimilation paramaters (t_bef, t_dur, t_aft, t_trig, tot_time)
spikes.vs_label = cell(7,1);
spikes.vs_values = zeros(7,1);
spikes.vs_label(:, 1) = {'t_trig' 't_bef' 't_dur' 't_aft' 'vs_start' 'vs_end' 'tot_time'};
spikes.vs_values(1:4, 1) = [spikes.vs_params(1, 9)/1000 spikes.vs_params(1, 6)/1000 spikes.vs_params(1, 7)/1000 spikes.vs_params(1, 8)/1000];
spikes.vs_values(5, 1) = [(spikes.vs_params(1, 9) + spikes.vs_params(1, 6))/1000];
spikes.vs_values(7, 1) = [sum(spikes.vs_params(1, 6:9))/1000];
spikes.vs_values(6, 1) = [spikes.vs_values(5, 1) + spikes.vs_values(3, 1)];

%Check wether each cluster violates the refractory period and assign a
%number to quantify the 'contamination' of the cluster;




%save spike construct in kilosort path
save([kilosortdata_path 'spikes.mat'], 'spikes');
end



