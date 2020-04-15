function [summary_led, summary_non_led] = tuning_distribution_AL(tuning, jitter_cutoff)
%% analysis parameter
% jitter_cutoff: the cut off for jiter in (ms)
%led=unique(spikes.led);% for tagging recordings in compute_tuning
if nargin<3
   jitter_cutoff = 2;
end

if nargin < 3
   include_mediocre_flag = 0;
end

led = 1;
DSI_cutoff = 0.1;
led_driven_alpha = 0.01;
 %include 

% if nargin < 3
%    tuning_firing_rate_min = 0.005; 
% end
%% get compiled data
datapath='E:\extarcellular\Andreanne\KilosortData';%parent folder of all the data

data_table = readtable('E:\extarcellular\Andreanne\data_analysis\phototagging_recordings_summary.xlsx');

switch tuning
    case 'sf'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 1);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 1, :);
    case 'tf'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 3);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 3, :);
    case 'ori'
        recordinglist = data_table.recordings(data_table.include == 1 & data_table.tuning == 2);
        data_table = data_table(data_table.include == 1 & data_table.tuning == 2, :);
end

%get spikes.mat path
data2compile={};
for i = 1: length(recordinglist)
        data2compile{end+1} = [datapath '\' recordinglist{i} '\' 'spikes.mat']; 
end

% confirm tuning type
tuning_data={tuning};
for i = 1:length(data2compile) 
    fullpath = data2compile{i};
    hyphen=strfind(fullpath,'-');%find hyphen position
    tuning_temp = fullpath(hyphen(4)+1:hyphen(4)+2);%extract tuning information from text after the last underscore
    if tuning_temp == 'or'
        tuning_temp='ori';
    end
    tuning_data{end+1}=tuning_temp;
end

if length(unique(tuning_data))>1 %check if all recordings are of the same type of tuning stimulus
   error('ERROR: tunings in the recordinglsit do not agree')
end
%% compiling tuning data (max tuning distribution and average tuning plot)
%max tuning: the orientation with the maximum firing rate for each cluster
%average tuning is the evoked firing rate

preferred_stim_compile_low=[]; %stimuli where max_tuning happens for cluster with jitter lower than the cut-off
preferred_stim_compile_high=[];%stimuli where max_tuning happens cluster with jitter higher than the cut-off
tuning_curve_norm_compile_low=[]; %normalized tuning curve for cluster with jitter lower than the cut-off
tuning_curve_norm_compile_high=[];%normalized tuning curve for cluster with jitter higher than the cut-off
DSI_compile_low = [];
DSI_compile_high = [];%ratio to examine the amplitude/sig 2 noise ratio of direction selectivity
cluster_ID_low = [];
cluster_ID_high = [];
recording_complile_low = table();
recording_complile_high = table();
num_cluster_low=0;
num_cluster_high=0;
is_vis_dri_compile_low = [];
is_vis_dri_compile_high = [];

for i=1:length(data2compile)
    load(data2compile{i});%load spike construct to workspace
    
    %filter spikes stuct
    spikes = is_visual_driven_spikes(spikes, tuning); %%MAYBE ADD A VIS DRIVEN FILTER TOO;
    spikes = is_led_driven(spikes, led_driven_alpha); %add .is_led_driven field to struct; alpha =0.01
    rem_depth = str2num(data_table.rem_chan{i});
    rem_channels = depth2chan(spikes, rem_depth, include_mediocre_flag);
    channels = unique(spikes.channel);
    good_channels = channels(~ismember(channels, rem_channels));
    spikes = filtspikes(spikes, 0, 'channel', good_channels);%filter for clusters that are inside notdtn
    good_assings = unique(spikes.assigns);
    
    %seperate clusters by jitter and led driven activity
    cluster_low=spikes.latency((spikes.latency(:,2)<jitter_cutoff & spikes.is_led_driven(:,2)) & spikes.is_visual_driven(:,2),1);%get clusterID for cluster with jitter lower than the cut-off
    cluster_high=spikes.latency(~spikes.is_led_driven(:,2) & spikes.is_visual_driven(:,2), 1);%get clusterID for cluster with non led-driven activity
    
    %remove cluster that were excluded from analysis (spike filt does not filter this field)
    cluster_low = cluster_low(ismember(cluster_low, good_assings)); 
    cluster_high = cluster_high(ismember(cluster_high, good_assings));
    
    norm_tuning_curve_low=[];%each line is the normalized tuning curve for one unit
    norm_tuning_curve_high=[];
    DSI_low=nan(length(cluster_low),2); %ratio for ds 
    DSI_high=nan(length(cluster_high),2);
    preferred_stim_low=nan(length(cluster_low),2);%prefered tuning 
    preferred_stim_high=nan(length(cluster_high),2);
    is_vis_driven_low = [];
    is_vis_driven_high = [];
  
    for m = 1:length(cluster_low)
        [stim_cond_list, tuning_curve_fr_total, tuning_curve_fr_spon] = compute_tuning(spikes, cluster_low(m), tuning, 0 , led);%0 for do not plot
        tuning_curve_fr_evoked = tuning_curve_fr_total-tuning_curve_fr_spon; %remove 'spontanous spikes' 
        if any(tuning_curve_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
           tuning_curve_fr_evoked = tuning_curve_fr_evoked-min(tuning_curve_fr_evoked);
           fprintf('tuning_fr_evoked for jitter_low contains negative value, force minimum to be zero\n');
        end

        switch tuning
            case {'sf' , 'tf'} %normalization and calculate preffered tuning value (ie max)
                [max_fr, max_fr_ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr;
                norm_tuning_curve_low = [norm_tuning_curve_low;tuning_fr_evoked_norm'];%each line is the normalized tuning for one unit
                preferred_stim_low(m,1:2)= [cluster_low(m) stim_cond_list(max_fr_ind)];%last column should be full of nan (needed for ori but not sf or tf)
                [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, cluster_low(m), tuning);
                is_vis_driven_low = [is_vis_driven_low; visual_activity_flag];
            case 'ori'
                %save normalized tuning curve fr for each unit
                [max_fr max_fr_ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr; %normalizing fr for all orientations by dividing by preferred direction
                norm_tuning_curve_low = [norm_tuning_curve_low; tuning_fr_evoked_norm']; %concatylating cluter orientation fr with all other units
                
                %calculate visual driven activity
                [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, cluster_low(m), tuning);
                is_vis_driven_low = [is_vis_driven_low; visual_activity_flag];
                
                %calculate DSI
                pref_stim_angle = stim_cond_list(max_fr_ind);
                null_angle = wrapTo360(pref_stim_angle + 180);
                    if null_angle > 330 %stim condition list from 0 to 330; to prevent error when 360 occurs:
                        null_ind = stim_cond_list == 0 ;
                    else
                        null_ind = stim_cond_list == null_angle ;
                    end
                null_fr = tuning_curve_fr_evoked(null_ind);
                DSI = (max_fr-null_fr)/(max_fr+null_fr);
                DSI_low(m,:) = [cluster_low(m) DSI];
                
                % above DSI threshhold calculate preferd direction of unit
                if DSI > DSI_cutoff

                    %calculate prefered direction by converting vectors into complex numbers
                    tuning_curve_fr_real = cosd(stim_cond_list)'.*tuning_curve_fr_evoked';
                    tuning_curve_fr_imag = sind(stim_cond_list)'.*tuning_curve_fr_evoked';
                    tuning_curve_fr_complex = complex(tuning_curve_fr_real, tuning_curve_fr_imag);
                    preferred_stim_complex = sum(tuning_curve_fr_complex); %adding up all each stimulus vector (in complex number format) to find preferred direction vector
                    pref_stim_angle = rad2deg(wrapTo2Pi(angle(preferred_stim_complex)));%extract the angle in deg from the complex number
                    
                    %save
                    preferred_stim_low(m,:)= [cluster_low(m) pref_stim_angle]; %save pref angle 
                    
                else
                    preferred_stim_low(m,:)= [cluster_low(m) nan];
%                     fprintf(['excluded unit: ' num2str(cluster_low(m)) '\n']);
                end
            end
    end 

    
    for n = 1:length(cluster_high)
        [stim_cond_list, tuning_curve_fr_total, tuning_curve_fr_spon] = compute_tuning(spikes, cluster_high(n), tuning, 0 , led);
        tuning_curve_fr_evoked = tuning_curve_fr_total-tuning_curve_fr_spon;
        if any(tuning_curve_fr_evoked<0)
           %if there is a negative value, force the smallest value to be 0
           tuning_curve_fr_evoked = tuning_curve_fr_evoked-min(tuning_curve_fr_evoked);
           fprintf('tuning_fr_evoked for non led-driven units contains negative value, force minimum to be zero\n');
        end
        switch tuning
            case {'sf' , 'tf'} %normalization and calculate preffered tuning value (ie max)
                [max_fr max_fr_ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr;
                norm_tuning_curve_high = [norm_tuning_curve_high;tuning_fr_evoked_norm'];%each line is the normalized tuning for one unit
                preferred_stim_high(n,1:2)= [cluster_high(n) stim_cond_list(max_fr_ind)];%last column should be full of nan (needed for ori but not sf or tf)
                [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, cluster_high(n), tuning);
                is_vis_driven_high = [is_vis_driven_high; visual_activity_flag];
            
            case 'ori'
                %save normalized tuning curve fr for each unit
                [max_fr max_fr_ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr; %normalizing fr for all orientations by dividing by preferred direction
                norm_tuning_curve_high = [norm_tuning_curve_high; tuning_fr_evoked_norm']; %concatylating cluter orientation fr with all other units
                
                %calculate visual driven activity (paired t-test)
                [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, cluster_high(n), tuning);
                is_vis_driven_high = [is_vis_driven_high; visual_activity_flag];
                
                %calculate DSI
                pref_stim_angle = stim_cond_list(max_fr_ind);
                null_angle = wrapTo360(pref_stim_angle + 180);
                    if null_angle > 330 %stim condition list from 0 to 330; to prevent error when 360 occurs:
                        null_ind = stim_cond_list == 0 ;
                    else
                        null_ind = stim_cond_list == null_angle ;
                    end
                null_fr = tuning_curve_fr_evoked(null_ind);
                DSI = (max_fr-null_fr)/(max_fr+null_fr);
                DSI_high(n,:) = [cluster_high(n) DSI];
                
                % above DSI threshhold calculate preferd direction of unit
                if DSI > DSI_cutoff

                    %calculate prefered direction by converting vectors into complex numbers
                    tuning_curve_fr_real = cosd(stim_cond_list)'.*tuning_curve_fr_evoked';
                    tuning_curve_fr_imag = sind(stim_cond_list)'.*tuning_curve_fr_evoked';
                    tuning_curve_fr_complex = complex(tuning_curve_fr_real, tuning_curve_fr_imag);
                    preferred_stim_complex = sum(tuning_curve_fr_complex); %adding up all each stimulus vector (in complex number format) to find preferred direction vector
                    pref_stim_angle = rad2deg(wrapTo2Pi(angle(preferred_stim_complex)));%in deg
                    
                    %save
                    preferred_stim_high(n,:)= [cluster_high(n) pref_stim_angle]; %save pref angle 
                    
                else
                    preferred_stim_high(n,:)= [cluster_high(n) nan];
%                     fprintf(['excluded unit: ' num2str(cluster_low(m)) '\n']);
                end
        end
    end
    
    %compile data from multiple recordings
    cluster_ID_low = [cluster_ID_low; cluster_low];
    cluster_ID_high = [cluster_ID_high; cluster_high];
    preferred_stim_compile_low=[preferred_stim_compile_low; preferred_stim_low]; 
    preferred_stim_compile_high=[preferred_stim_compile_high; preferred_stim_high];
    num_cluster_low=num_cluster_low +length(cluster_low);
    num_cluster_high=num_cluster_high +length(cluster_high);
    tuning_curve_norm_compile_low=[tuning_curve_norm_compile_low; norm_tuning_curve_low];
    tuning_curve_norm_compile_high=[tuning_curve_norm_compile_high; norm_tuning_curve_high];
    DSI_compile_low = [DSI_compile_low; DSI_low];
    DSI_compile_high=[DSI_compile_high; DSI_high];
    is_vis_dri_compile_low = [is_vis_dri_compile_low; is_vis_driven_low];
    is_vis_dri_compile_high = [is_vis_dri_compile_high; is_vis_driven_high];

    
    % make summary table
    if ~ isempty(cluster_low)
        temp_recording_low = table(cellstr(repmat(recordinglist{i}, length(cluster_low), 1)));
        recording_complile_low = [recording_complile_low; temp_recording_low];
    end
    
    if ~ isempty(cluster_high)
        temp_recording_high = table(cellstr(repmat(recordinglist{i}, length(cluster_high), 1)));
        recording_complile_high = [recording_complile_high; temp_recording_high];
    end
    
end

%% calculate values for figures

%store values in table
%rename var to make table headers cleaner
    recording = recording_complile_low.Var1;
    cluster = cluster_ID_low;
    pref_stim = preferred_stim_compile_low(:,2);  
    tuning_curve=tuning_curve_norm_compile_low;
    DSI = DSI_compile_low(:,2);
    is_vis_dri = is_vis_dri_compile_low;
    summary_led = table(recording, cluster, tuning_curve, pref_stim, DSI, is_vis_dri); %create table
    
    recording = recording_complile_high.Var1;
    cluster = cluster_ID_high; %first column is cluster ID
    pref_stim= preferred_stim_compile_high(:,2);
    tuning_curve= tuning_curve_norm_compile_high;
    DSI = DSI_compile_high(:,2);
    is_vis_dri = is_vis_dri_compile_high;
    summary_non_led = table(recording, cluster, tuning_curve, pref_stim, DSI, is_vis_dri);

    

%remove inf and convert to nan
tuning_curve_norm_compile_low(find(isinf(tuning_curve_norm_compile_low))) = nan;
tuning_curve_norm_compile_high(find(isinf(tuning_curve_norm_compile_high))) = nan;

%calculating averaged tuning
ave_tuning_norm_compile_low=nanmean(tuning_curve_norm_compile_low); %remove nan
ave_tuning_norm_compile_high=nanmean(tuning_curve_norm_compile_high);  %remove nan
ave_tuning_norm_compile_low=ave_tuning_norm_compile_low/max(ave_tuning_norm_compile_low);%normalized against the max mean
ave_tuning_norm_compile_high=ave_tuning_norm_compile_high/max(ave_tuning_norm_compile_high);

%% PLOT FIGURES
switch tuning
    case 'sf'
        %plot pperecentage of preffered stim
        h = figure;
        set(h,'position', [100, 100, 1200, 900]);%set figure position and size
        C = categorical({'0.04', '0.08', '0.16', '0.32', '0.45'});
        subplot(1,3,1)%plot distribution
        %by default clusters with jitter<2ms is blue, the clusters with jitter>2ms is red
%         bar(C, [percent_pref_stim_low(:,2) percent_pref_stim_high(:,2)])
%         %bar(C,[percent_pref_stim_low(:,2) percent_pref_stim_high(:,2)], 1)
%         xlabel('spatial frequency', 'FontSize', 25)
%         ylabel('% units', 'FontSize', 25)
%         title('percenage of perferred spatial frequency for each population')
        
        %plot relative cummulative graph
        hold on
        ecdf(preferred_stim_compile_low(:,2))
        ecdf(preferred_stim_compile_high(:,2))
        
        subplot(1,2,2)%plot average normalized tuning curve
        plot(stim_cond_list', ave_tuning_norm_compile_low,'b','linewidth',2)
        hold on
        plot(stim_cond_list', ave_tuning_norm_compile_high,'r','linewidth',2)
        xlabel('spatial frequency', 'FontSize', 25)
        ylabel('norm firing rate', 'FontSize', 25)
        title('averaged normalize tuning curve')
        
        %set legend position
        hL = legend(['jitter < ' num2str(jitter_cutoff) 'ms'],['non led-driven units']);
        legendPosition = [0.8 0.95 0.2 0.05]; 
        set(hL, 'Position', legendPosition);
        suptitle('spatial frequency');
    
        %STATISTICS
        [P, H] = ranksum(preferred_stim_compile_low(:,2), preferred_stim_compile_high(:,2));
        
        %SAVE
        save(['E:\extarcellular\Andreanne\data_analysis\photo-tag\population tuning\summary_led_' tuning datestr(now, 'dd-mmm-yyyy') user '.mat'], 'summary_led', 'summary_non_led', ...
            'preferred_stim_compile_low', 'preferred_stim_compile_high', 'stim_cond_list', 'ave_tuning_norm_compile_low', 'ave_tuning_norm_compile_high', 'jitter_cutoff');
        
    % case 'tf'
    
    case 'ori'
        h = figure;
        set(h,'position', [100, 100, 1200, 900]);%set figure position and size
       
        % histogram percentage preferred stim
        subplot(1,3,1)
        n_bins = [(stim_cond_list - 15); stim_cond_list(end) + 15];
        hold on
        histogram(preferred_stim_compile_high(:,2), n_bins, 'FaceColor', [1 0 0], 'Normalization','probability', 'FaceAlpha',0.8)
        yticklabels(yticks*100)
        histogram(preferred_stim_compile_low(:,2), n_bins , 'FaceColor', [0 0 1], 'Normalization','probability', 'FaceAlpha',0.6)
        yticklabels(yticks*100)


        ylabel('% units');
        xlabel('direction selectivity (degrees)');

        %polarplot of average tuning curve
        subplot(1,3,2);
                
        theta = stim_cond_list'/180*pi;
        theta = [theta theta(1)];
        
        ave_tuning_norm_compile_low = [ave_tuning_norm_compile_low ave_tuning_norm_compile_low(1)];
        ave_tuning_norm_compile_high = [ave_tuning_norm_compile_high ave_tuning_norm_compile_high(1)];
        polar(theta, ave_tuning_norm_compile_low,'b')  
        hold on
        polar(theta, ave_tuning_norm_compile_high,'r')    
        title('Averaged direction tuning curve', 'FontSize', 15)
        text(-0.4, 0.55, 'normalized fr');
  
        %set legend position
        hL = legend(['jitter < ' num2str(jitter_cutoff) 'ms'], ['non led-driven units']);
        legendPosition = [0.4 0.8 0.2 0.1]; 
        set(hL, 'Position', legendPosition)  
        suptitle('Direction Selectivity of Inhibitory NOT/DTN neurons');
        
        %bar graff of DSI; 
        preferred_over_null_compile_high_mean=nanmean(DSI_compile_high(:,2));
        preferred_over_null_compile_high_std=nanstd(DSI_compile_high(:,2));
        preferred_over_null_compile_high_stem=preferred_over_null_compile_high_std/sqrt(length(DSI_compile_high(:,2)));
        
        preferred_over_null_compile_low_mean=nanmean(DSI_compile_low(:,2));
        preferred_over_null_compile_low_std=nanstd(DSI_compile_low(:,2));
        preferred_over_null_compile_low_stem=preferred_over_null_compile_low_std/sqrt(length(DSI_compile_low(:,2)));

        [H, P] =ttest2(DSI_compile_high(:,2), DSI_compile_low(:,2));
       
        %plot DSI
        subplot(1,3,3)
        hold on
        bar((1), [preferred_over_null_compile_high_mean], 'r'); 
        errorbar((1), [preferred_over_null_compile_high_mean], [preferred_over_null_compile_high_stem], 'black');
        bar((0), [preferred_over_null_compile_low_mean], 'b'); 
        errorbar((0), [preferred_over_null_compile_low_mean], [preferred_over_null_compile_low_stem], 'black');
        
        title('Direction Selectivity Index')
        ylabel('(prefferd - null)/(preferred + null)')
        set(gca,'xtick',[]);
        xlabel('NOT/DTN populations');
        
        if H
            text(-0.4, max(preferred_over_null_compile_high_mean + preferred_over_null_compile_high_stem*1.1), '----------- * -----------')
        else
            text(-0.4, max(preferred_over_null_compile_high_mean + preferred_over_null_compile_high_stem*1.1), '----------- n.s.-----------')
        end

        %% save
       save(['E:\extarcellular\Andreanne\data_analysis\photo-tag\population tuning\summary_led_' tuning datestr(now, 'dd-mmm-yyyy') user '.mat'], 'data_table', 'summary_led', 'summary_non_led', ...
            'stim_cond_list', 'ave_tuning_norm_compile_low', 'ave_tuning_norm_compile_high', 'jitter_cutoff');
       
       savefig(['E:\extarcellular\Andreanne\data_analysis\photo-tag\population tuning\' tuning '-jitter-' num2str(jitter_cutoff) '-' datestr(now, 'dd-mmm-yyyy')]);
end
end

