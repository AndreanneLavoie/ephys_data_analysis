function spikes = create_spikes_merge(kilosort_data_path, probe)
%this code will create a structure called spikes which will hold all of the
%data required for data analysis after Phy (manual spike sorting)

%data from different sources should be collected together in a folder three seperate folders, 
%organised by location/software used to generate the data; 

    % This includes: 
% -> From IntanData -- folder containing raw ephys data, and vs/led timing ('info.rhd', 'time.dat', 'board-ADC-00.dat' (00-04)(record of analog inputs for different stim; led, eye track, vs timing, behaviour etc.)
% -> From Labview -- order of stimuli ( '20trk1DATE_expName.eye' ),
% -> From Kilosort -- 'spike_times.npy', 'spike_clusters.npy' ETC.
% -> From Phy --  'cluster_groups.csv'

%must have .npy to matlab installed from https://github.com/kwikteam/npy-matlab

%%%

%location of Intan data stored: 
%[intandata_path, frequency_parameters] = read_Intan_RHD2000_file; %the chosen folder path wil have the same exp name as vs
percT_thres = 0.15;% 15% of the led / vs duration difference for ramping down
preset_trigdur = 1/60*3;%60 is the frame rate of the monitor which delivers vs, trigger lasts 3 frames

%threshold and digital filter for led, vs and trigger signals:
threshold=1.8;%used to detect rising phase in vs_time2 function
digfilter_winT_fast=0.0006; %in sec
digfilter_winT_slow=0.001; %in sec

%% extract experiment name and extrapolate path for different types of data using the Intan info.rhd file selected;
spikes.kilosort_path = [kilosort_data_path '\'];

%find index of key markers in the intandata filename
gen_path_ind = strfind(spikes.kilosort_path, 'KilosortData'); %28 (index of 'K' in KilosortData) 
general_path = spikes.kilosort_path(1:gen_path_ind - 1); %'E:\extarcellular\Andreanne\'
kilo_name_len = length('KilosortData');
exp_name = spikes.kilosort_path(gen_path_ind + kilo_name_len + 1 : end -1); %everything after 'KilosortData\' except the last \ (hence he -1)

%each recording session name (ie. sf1 or ori3 etc.)
unmerged_exp_name_ind = strfind(exp_name,'AND');
num_merged_files = length(unmerged_exp_name_ind) + 1;
exp_session_name_ind = strfind(exp_name,'-');
exp_session_name_ind = exp_session_name_ind(end) + 1;
exp_name_base = exp_name(1:exp_session_name_ind -2);
tuning = {};
recording_sess_exp_name = {};
for i = 1:num_merged_files
    if i == 1
        tuning{i} = exp_name(exp_session_name_ind :unmerged_exp_name_ind(i)-1);
    elseif i > length(unmerged_exp_name_ind)
        tuning{i} = exp_name(unmerged_exp_name_ind(end)+3:end);
    else
        tuning{i} = exp_name(unmerged_exp_name_ind(i)+3:unmerged_exp_name_ind(i+1));
    end
    recording_sess_exp_name{i} = [exp_name_base '-' tuning{i}];
end

%grab on of the original raw intan data path to get date in order to find
%extract vis stim parameter file name (which has data in title)
cd(fullfile(general_path, 'IntanData'));
raw_path = dir([exp_name(1:unmerged_exp_name_ind-1) '*.*']);
exp_ephys_ind = strfind(raw_path(1).name, '_ephys'); %index of '_' in _ephys 
date = raw_path(1).name((exp_ephys_ind + 6):(exp_ephys_ind + 12)); %e.g. '_190223'
vspara_path = [general_path 'VSpara\'];

%concatilating vsparams from different recordings
spikes.recoding_sess_vsparam ={};

for j = 1:length(recording_sess_exp_name)
    %grab proper vsparam; if comes from eye tracking computer .eye format; if
    %comes from vis stim computer uses .csv format (only for orientation tuning)
    if ~strcmp(tuning{j}(1:2), 'or') %for sf and tf tuning
        vspara_filename = ['20trk' date(6:7) date(4:5) date(2:3) '_' recording_sess_exp_name{j}];
    
        %load the vs stimulation parameters depending on wether a recording has 
        %multiple vspara files(e.g. led condition change within trial)
        if isfile([vspara_path vspara_filename '.eye']) % a single vspara file for recording
            [~,spikes.labview_param] = getPsychStimParameters_BL(vspara_path, [vspara_filename '.eye']);
            temp = spikes.labview_param.stimpara';
            spikes.recoding_sess_vsparam{j} = temp(1:spikes.labview_param.stimnum, :);
        else
            cd(vspara_path)
            all_files = dir; %information for all folders in current directory
            file_list=string({all_files(~[all_files.isdir]).name})';%the names for all files in the directory
            if any(startsWith(file_list,string(vspara_filename)))%list folder names for all the recodings starts with 'vgt' or 'c57'
                vspara_files_name=file_list(startsWith(file_list,string(vspara_filename)));
                vspara_files={};
                for i=1: length(vspara_files_name)
                    temp=dir(vspara_files_name{1});
                    vspara_files = [vspara_files struct2cell(temp)];%after converting to cell, each column is a file
                    %row: (1)name (2)folder (3)date (4)bytes (5)isdir (6)datenum
                end
                [~,vs_files_ind] = sort(datetime(vspara_files(3,:)),'ascend');%sort by time from earlier to later
                %vs_files_ind saves the index of vspara_files from earliest (left) to latest(right)
                spikes.vs_params=[];
                for k = 1:length(vs_files_ind)
                    %save labview parameter structure for all stimulus conditions under different subfields in labview_param structure
                    [~,spikes.labview_param.(['condition' num2str(vs_files_ind(k))])] = getPsychStimParameters_BL(vspara_path, char(vspara_files(1,vs_files_ind(k))));
                    %concatenate vs_params files in time order from the earliest to the latest 
                    % vertically
                    vs_params_temp=spikes.labview_param.(['condition' num2str(vs_files_ind(k))]).stimpara';
                    vs_params_temp=vs_params_temp(1:spikes.labview_param.(['condition' num2str(vs_files_ind(k))]).stimnum,:);
                    vs_params_temp(:, end + 1) = k;
                    spikes.recoding_sess_vsparam{j} =[spikes.vs_params; vs_params_temp];
                end
            else
                f=msgbox(['no vs file exist for ' recording_sess_exp_name{j}])
                return
            end
        %%FOR RECORDINGS THAT LED CONDITION IS NOT SAVED IN vs_param COLUMN 17 yet
        %different led power used within recording; low = 1.9 V; med = 3.7; V; high = 5 V;
        %[~,labview_param_low] = getPsychStimParameters_BL(vspara_path, [vspara_filename '-low.eye']);
%         [~,labview_param_med] = getPsychStimParameters_BL(vspara_path, [vspara_filename '-med.eye']);
%         [~,labview_param_med_decay] = getPsychStimParameters_BL(vspara_path, [vspara_filename '-med-d.eye']);
%         [~,labview_param_high] = getPsychStimParameters_BL(vspara_path, [vspara_filename '-high.eye']);
%         [~,labview_param_high_decay] = getPsychStimParameters_BL(vspara_path, [vspara_filename '-high-d.eye']);
%         %vs_params_low = labview_param_low.stimpara';
%         vs_params_med = labview_param_med.stimpara';
%         labview_param_med_decay = labview_param_med_decay.stimpara';
%         vs_params_high = labview_param_high.stimpara';
%         labview_param_high_decay = labview_param_high_decay.stimpara';
%         %vs_params_low = vs_params_low(1:labview_param_low.stimnum,:);
%         vs_params_med = vs_params_med(1:labview_param_med.stimnum,:);
%         vs_params_med_decay = vs_params_med_decay(1:labview_param_med_decay.stimnum,:);
%         vs_params_high = vs_params_high(1:labview_param_high.stimnum,:);
%         vs_params_med_decay = vs_params_med_decay(1:labview_param_med_decay.stimnum,:);
%         %add extra column to mark which power was used; low = 1; med = 2;high = 3 
%         vs_params_low(:, end+1) = 1; 
%         vs_params_med(:, end+1) = 2; 
%         vs_params_high(:, end+1) = 3; 
%         %concatilate all three vs_params
%         spikes.vs_params = [vs_params_low; vs_params_med; vs_params_high];
        end
    elseif strcmp(tuning{j}(1:2), 'or')
        
    %for ori / direction selectivity parameters
    vspara_filename = erase(recording_sess_exp_name{j},'-');
    
    if isfile([vspara_path vspara_filename '.csv'])   
        spikes.vs_params = csvread([vspara_path vspara_filename '.csv']);
    else
        cd(vspara_path)
        all_files = dir; %information for all folders in current directory
        file_list=string({all_files(~[all_files.isdir]).name})';%the names for all files in the directory
        if any(startsWith(file_list,string(vspara_filename)))
            vspara_files_name=file_list(startsWith(file_list,string(vspara_filename)));
            vspara_files={};
            for i=1: length(vspara_files_name)
                temp=dir(vspara_files_name{i});
                vspara_files=[vspara_files; struct2cell(temp)'];%after converting to cell, each column is a file
                %row: (1)name (2)folder (3)date (4)bytes (5)isdir (6)datenum
            end        
            [~,vs_files_ind] = sort(datetime(vspara_files(:,3)),'ascend');%sort by time from earliest to latest
            
            spikes.vs_params=[];
            if length(vs_files_ind)>1
                for k = 1:length(vs_files_ind)
                    %save labview parameter structure for all stimulus conditions under different subfields in labview_param structure
                    vs_params_temp = csvread([vspara_path vspara_files{vs_files_ind(k),1}]);
                    %concatenate vs_params files in time order from the earliest to the latest 
                    % vertically
                    vs_params_temp(:, end + 1) = k;
                    spikes.vs_params=[spikes.vs_params; vs_params_temp];
                end
            else
                spikes.recoding_sess_vsparam{j}  = csvread([vspara_path vspara_files{1,1}]);
            end
        else 
            f=msgbox(['no vs file exist for ' recording_sess_exp_name{j}])
            return 
        end
        
    end
    end
end

 %final concatilation
 spikes.vs_params =[];
for i=1:length(spikes.recoding_sess_vsparam)
    dimention_rec_sess_col = size(spikes.recoding_sess_vsparam{i});
    spikes.vs_params = [spikes.vs_params; spikes.recoding_sess_vsparam{i}(:, 1:17) i*ones(dimention_rec_sess_col(1),1)];
end
 
% read in info.rhd 
spikes.intan_info=read_Intan_RHD2000_file_ana('info.rhd',[spikes.kilosort_path '\']);
amplifier_sample_rate =spikes.intan_info.frequency_parameters.amplifier_sample_rate;

threshold=1.8;%used to detect rising phase in vs_time2 function
digfilter_winT=0.0006; %in sec


%read kilosort and intan data files into spikes structures
spikes.abs_spiketimes = double(readNPY([spikes.kilosort_path 'spike_times.npy']))'/double(amplifier_sample_rate); %absolute time point for each spike   
spikes.vstiming = extract_signal_time('time.dat', 'board-ADC-00.dat', spikes.kilosort_path, threshold,digfilter_winT, amplifier_sample_rate);
spikes.trigtiming= extract_signal_time('time.dat', 'board-ADC-01.dat', spikes.kilosort_path, threshold,digfilter_winT, amplifier_sample_rate); %visual stim trigger times

%% Both spiketimes, trials and stimcond are genereated together; 
%spiketimes: relative time of spike to trigger;
%stimcond: the type of stimulus associated with each spike, given the closest trigger vsstim
%stimorder: trigger identity
spikes = make_spike_struct(spikes);
spikes.triggers = spikes.trials; %triggers = 1 to n trials throughout ALL merged recording sessions 
spikes.vs_types = zeros(size(spikes.triggers));
spikes.files = zeros(size(spikes.triggers));

recording_start = 1;
for n = 1:length(recording_sess_exp_name)
    recording_end = sum(spikes.vs_params(:,18)==n) + (recording_start - 1);
    index = find(spikes.triggers >= recording_start & spikes.triggers <= recording_end);
    
    % file number 
    if n == 1
        spikes.files(find(spikes.triggers <= recording_end & spikes.triggers ~= 0)) = n;
    elseif n == length(recording_sess_exp_name)
        spikes.files(find(spikes.triggers >= recording_start)) = n;
    else n ~= 1 & n ~= length(recording_sess_exp_name)
        spikes.files(index) = n;
    end

    spikes.trials(index) = spikes.trials(index) - (recording_start - 1) ; 
    
    switch tuning{n}
        case 'sf'
            spikes.vs_types(index) = 1;
        case 'ori'
            spikes.vs_types(index) = 2;
        case 'tf'
            spikes.vs_types(index) = 3;
    end
    
    recording_start = recording_end + 1;
end

%cluster / unit identity of each spikes
spikes.assigns = double(readNPY([spikes.kilosort_path 'spike_clusters.npy'])');

%% import LED info
if  sum(spikes.vs_params(:,17)) > 0
    spikes.led = nan(size(spikes.trials)); 
    spikes.led(spikes.trials > 0) = spikes.vs_params(spikes.trials(spikes.trials > 0), 17); %on = 1; off = 0;
    if abs(spikes.vs_params(1,16)-1)<10e-8
        %spikes.led_ID(spikes.trials > 0) = led_ID(spikes.trials(spikes.trials > 0));
        spikes.led_abstime = extract_signal_time('time.dat', 'board-ADC-02.dat', spikes.kilosort_path, threshold, digfilter_winT_slow, amplifier_sample_rate); %visual stim trigger times
    else
        spikes.led_abstime = extract_signal_time('time.dat', 'board-ADC-02.dat', spikes.kilosort_path, threshold, digfilter_winT_fast, amplifier_sample_rate); %visual stim trigger times
    end
        %[length_led, spikes.led_spiketimes, led_hist] = make_first_led_spiketimes(spikes, [], 1, [025], 1, 0);
   if any(abs(spikes.led_abstime(:,3)-spikes.vs_params(1,12)/1000)>spikes.vs_params(1,12)*percT_thres/1000)
      fprintf ('ERROR in led time calculation, led period longer than set\n'); %ERROR CHECK for led time calculation 
   end
   spikes = make_led_spike_struct(spikes, []);%assign spikes to led events and calculate the relative time from the assigned spikes
end
%% Extract spikes.cluster and spikes.labels field that marks cluster as either: good, mua, noise 
% save cluser info(id, amplitude, channel, depthm firing_rate (spk/s),group id, number of spikes and shank) in a structure
% table format cause error when doing filtspikes
% all numbers (id, amplitude, channel, depthm firing_rate, n_spikes, shank) are in double format, the group is a char array 
if isfile([spikes.kilosort_path 'cluster_info.tsv'])
    spikes.clusterinfo = tdfread([spikes.kilosort_path 'cluster_info.tsv'],'tab');
    spikes.clusterinfo.firing_rate=str2double(regexp(string(spikes.clusterinfo.firing_rate), '[\d.]+','match'));
    
    %create spikes.channel structure (for each spiketime assign which
    %channel it has max amplitude
    load(['E:\extarcellular\Andreanne\chanMaps\' probe])
    spikes.channel = spikes.assigns;
    clusters = spikes.clusterinfo.id;
    for i = 1:length(clusters)
        spikes.channel(abs(spikes.assigns - clusters(i)) < 1e-6) = Intan_chan_num((spikes.clusterinfo.channel(i)+1)); %both phy channel info and Intan channels are 0 based, but need to be 1 based for indexing
    end
end

%assign corresponding depth for each spiketime:
%%TO ADD

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
    formatSpec = '%s%s';
    csvpath=[spikes.kilosort_path 'cluster_groups.csv'];
    T = readtable(csvpath,'Format',formatSpec);
elseif isfile([spikes.kilosort_path 'cluster_groups.csv'])
    formatSpec = '%s%s';
    csvpath=[spikes.kilosort_path 'cluster_groups.csv'];
    T = readtable(csvpath,'Format',formatSpec);
    %take out unsorted clusters in csv (Phy1 data)
    unsortedflag=strcmp(T.group,'unsorted');
    T(unsortedflag,:)=[];
else
    T=table;%initiate empty table when none exist
    csvpath=[spikes.kilosort_path 'cluster_groups.csv'];
end

numofclusVI=sum(refviolflag);
if isempty(T)
    if numofclusVI>0 && numofclusVI>0
        clusIDVI=spikes.labels(find(refviolflag),1);
        cluster_id=strsplit(num2str(clusIDVI'))';
        group=cell(size(cluster_id,1),1);
        group(1:end,1)={'mua'};
        T=table(cluster_id, group);
    end
else
    if numofclusVI>0
        clusIDVI=spikes.labels(find(refviolflag),1);
        clusIDtsv=str2num(char(T.cluster_id));
        grouptsv=T.group;
        for i=1:numofclusVI
            Subtr=clusIDtsv - clusIDVI(i);
            if ~any(abs(Subtr)<1e-9)% check if any element in sub equals zero
                if  any(Subtr>0)>0
                    Ind=min(find(Subtr>0));
                    clusIDtsv=[clusIDtsv(1:(Ind-1)); clusIDVI(i); clusIDtsv(Ind:end)];
                    grouptsv=[grouptsv(1:Ind-1); {'mua'}; grouptsv(Ind:end)];
                else
                    clusIDtsv=[clusIDtsv; clusIDVI(i)];
                    grouptsv=[grouptsv; {'mua'}];
                end
            end
        end
        T.cluster_id(1:length(clusIDtsv))=strsplit(num2str(clusIDtsv'));
        T.group(1:length(clusIDtsv))=grouptsv;
    end
end

%updata spikes.labels after mannual labeling
Ttemp=T;
mis_num=setdiff(spikes.labels(:,1),str2double(string(Ttemp.cluster_id)));
for i=1:length(mis_num)
    mis_ind=find(abs(spikes.labels(:,1)-mis_num(i))<1e-8);
    Ttemp=[Ttemp(1:mis_ind-1,:); {'unsorted' 'unsorted'}; Ttemp(mis_ind:end,:)];
end

noiseflag=strcmp(Ttemp{:,2},'noise'); 
goodflag = strcmp(Ttemp{:,2},'good');
muaflag = strcmp(Ttemp{:,2},'mua');
spikes.labels(find(goodflag),2) = 2; %label good before refviolflag to remove any dirty units previously labeled as good
spikes.labels(find(noiseflag),2) = 0;
spikes.labels(find(muaflag),2) = 1;


%create table and save as .csv
writetable(T, csvpath, 'Delimiter', 'tab');

%convert .csv to .tsv to where changes can be seen in phy2
if isfile([spikes.kilosort_path 'cluster_groups.tsv'])
    delete([spikes.kilosort_path 'cluster_group.tsv']);
end

copyfile([spikes.kilosort_path 'cluster_groups.csv'], [spikes.kilosort_path 'cluster_group.tsv']);
delete([spikes.kilosort_path 'cluster_groups.csv']);

%save spike construct in kilosort path
save([spikes.kilosort_path 'spikes.mat'], 'spikes');
end



