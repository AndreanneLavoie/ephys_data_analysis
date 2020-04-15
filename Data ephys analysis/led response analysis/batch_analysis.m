%%batch analysis template for auto-analyzing multiple recordings
%add and(or) edit recordings to be analyzed in the recordinglist
%add and(or) edit analyses to be done on each recording in the "analysis for each recording" section
%and and(or) edit analyses to be done on all recordings compiled together
%in the "compiled analysis for all the recordings" section

%% recordings list to be analyzed
datapath='';%parent folder of all the data
electrode_config= ; % the electrode configuration 
%update recording list for batch analysis
recordinglist={'vgt-0009-m0-t1-sf'...
               'vgt-0009-m0-t2-sf'...
               'vgt-0009-m0-t3-sf'...
               'vgt-0009-m0-t1-ori'...
               'vgt-0009-m0-t2-ori'...
               'vgt-0009-m1-t4-sf5'...
               'vgt-0009-m0-t1-sf2'...
               'vgt-0009-m0-t2-sf2'...
               'vgt-0009-m0-t1-sf3'...
               'vgt-0009-m1-t3-sf'...
               'vgt-0015-m0-t4-sf1'...
               'vgt-0015-m0-t5-sf1'...
               'vgt-0015-m0-t6-sf1'...
               'vgt-0015-m0-t4-ori'...
               'vgt-0015-m0-t5-ori'...
               'vgt-0015-m0-t6-ori'...
               'vgt-0015-m1-t3-sf'...
               'vgt-0015-m1-t4-sf'...
               'vgt-0015-m2-t1-sf'...
               'vgt-0015-m2-t7-ori'...
               'vgt-0016-m2-t3-ori1'...
               'vgt-0016-m2-t4-sf'...
               'vgt-0016-m2-t5-sf'...
               }';   
%get spikes.mat path
data2compile=strcat(datapath, '\', recordinglist, '\spikes.mat');

%% analysis for each recording
for i = 1:length(data2compile)
    
   %%load spike construct to workspace 
   load(data2compile{i});
    
   %%get tuning for each recording
   fullpath = data2compile{i};
   hyphen=strfind(fullpath,'-');%find hyphen position
   tuning = fullpath(hyphen(4)+1:hyphen(4)+2);%extract tuning information from text after the last underscore
   
   if tuning == 'or'
       tuning = 'ori';
   end
   
   %%update spikes construct
   spikes=create_spikes(spikes.raw_data_path);
   
 
   %%cluster analysis
   if  all(spikes.vs_params(:,17)) && abs(spikes.vs_params(1,13)-spikes.vs_params(1,12))> spikes.vs_params(1,12)*0.2
       tagging = 1;
   else
       tagging = 0;
   end
   
   spikes=cluster_analysis(spikes,[],tuning,electrode_config, tagging);
       
   if  tagging==1
       %%is_led_response analysis
       spikes.led_response_cluster=is_led_response(spikes, [], tuning);
   end
   
   %%save the spikes construct to its folder
   save([spikes.kilosort_path 'spikes.mat'], 'spikes');
end

%% compiled analysis for all the recordings (regardless of tuning type in the recording list)

% scatter plots of jitter over peakdis for all clusters with led response
binwidth=0.002;
led_jitter_peakdis=[];%(1)cluster id (2)jitter (3)peak distance
peak_dis_all=[];
for i = 1:length(data2compile)
    load(data2compile{i});
    if  sum(spikes.vs_params(:,17)) > 0
        if isfield(spikes, 'led_response_cluster')
            for j = 1:length(spikes.led_response_cluster)
                assign=spikes.led_response_cluster(j);
                jitter=spikes.latency(spikes.latency(:,1)==assign,2);
                peak_dis=spikes.peak_dis(spikes.peak_dis(:,1)==assign,2);
                led_jitter_peakdis = [led_jitter_peakdis; assign jitter peak_dis];%jitter and peak_dis for tagged clusters
            end
        end
    end
    if isfield(spikes, 'peak_dis')
        peak_dis_all=[peak_dis_all; spikes.peak_dis];%peak_dis for all clusters
    end
    
end

figure
subplot(1,2,1)%scatter plots of jitter vs. peak_dis for tagged cluster
scatter(led_jitter_peakdis(:,3),led_jitter_peakdis(:,2))
ylabel('jitter(ms)')
xlabel('spontaneous peak distance(ms)')
title('led responsive clusters')
subplot(1,2,2)%histogram for peak_dis distribution of all clusters
binsize=round((max(peak_dis_all(:,2))-min(peak_dis_all(:,2)))/binwidth);
hist(peak_dis_all(:,2),binsize);%in ms
xlabel('peak distance (ms)')
title('all clusters')