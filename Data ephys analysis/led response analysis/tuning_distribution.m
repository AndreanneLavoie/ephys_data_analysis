function [summary_led, summary_non_led] = tuning_distribution(user, tuning, jitter_cutoff)
%% analysis parameter
% jitter_cutoff: the cut off for jiter in (ms)
%led=unique(spikes.led);% for tagging recordings in compute_tuning
if nargin<3
   jitter_cutoff = 2;
end
led = 1;

% if nargin < 3
%    tuning_firing_rate_min = 0.005; 
% end
%% get compiled data
switch user
    case 'AL'
        datapath='E:\extarcellular\Andreanne\KilosortData';%parent folder of all the data
        %recording library, mannually add new recordings to list.
        switch tuning
            case 'sf'
                recordinglist={'vgt-0009-m0-t1-sf'...
                               'vgt-0009-m0-t2-sf'...
                               'vgt-0009-m0-t3-sf'...
                               'vgt-0009-m1-t4-sf5'...
                               'vgt-0020-m1-t4-sf'...
                               'vgt-0024-m2-t2-sf'...
                               'vgt-0024-m0-t3-sf'...
                               'vgt-0024-m0-t4-sf'...
                               'vgt-0024-m0-t5-sf'...
                               'vgt-0026-m0-t4-sf'...
                               };          
            case 'tf'
                recordinglist={};
            case 'ori'
                recordinglist={'vgt-0009-m0-t1-ori'...
                               'vgt-0009-m0-t2-ori'...
                               'vgt-0009-m0-t3-ori'...
                               'vgt-0020-m0-t10-ori'...
                               'vgt-0020-m1-t4-ori'...
                               'vgt-0024-m0-t3-ori'...
                               'vgt-0024-m0-t4-ori'...
                               'vgt-0024-m0-t5-ori'...
                               }; 
        end
    case 'HH'
        datapath='E:\extarcellular\Hayley\KilosortData';
        switch tuning
            case 'sf'
                recordinglist={
                               };          
            case 'tf'
                recordinglist={};
            case 'ori'
                recordinglist={'c57-0047-m1-t8-ori'...
                               'c57-0046-m1-t4-ori'...
                               'c57-0047-m1-t1-ori1'...
                               'c57-0047-m1-t8-ori'...
                               'c57-0047-m1-t9-ori'...
                               'c57-0047-m2-t4-ori'...
                               'c57-0047-m2-t5-ori'...
                               'c57-0045-m1-t3-ori'...
                               'c57-0045-m1-t4-ori'...
                               'c57-0045-m1-t8-ori'...
                               'c57-0045-m2-t4-ori'...
                               'c57-0045-m2-t8-ori'...
                               }; 
        end
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
preferred_over_null_compile_low = [];
preferred_over_null_compile_high = [];%ratio to examine the amplitude/sig 2 noise ratio of direction selectivity
cluster_ID_low = [];
cluster_ID_high = [];
recording_complile_low = table();
recording_complile_high = table();
num_cluster_low=0;
num_cluster_high=0;

for i=1:length(data2compile)
    load(data2compile{i});%load spike construct to workspace
    spikes = is_led_driven(spikes, 0.01) %add .is_led_driven field to struct; alpha =0.01
    %[spikes, spikes.is_notdtn, chan_depth_map_index, notdtn_range] = psth_depth_any_probe(spikes, tuning); %adds spikes.is_notdtn (and more fields and plots psth_depth)
    %spikes = filtspikes(spikes, 0, 'is_notdtn', 1);%filter for clusters that are inside notdtn
    cluster_low=spikes.latency((spikes.latency(:,2)<jitter_cutoff & spikes.is_led_driven(:,2)),1);%get clusterID for cluster with jitter lower than the cut-off
    cluster_high=spikes.latency(~spikes.is_led_driven(:,2),1);%get clusterID for cluster with non led-driven activity
    norm_tuning_curve_low=[];%each line is the normalized tuning curve for one unit
    norm_tuning_curve_high=[];
    preferred_over_null_low=nan(length(cluster_low),2); %ratio for ds 
    preferred_over_null_high=nan(length(cluster_high),2);
    preferred_stim_low=nan(length(cluster_low),2);%prefered tuning 
    preferred_stim_high=nan(length(cluster_high),2);
%     summary_clusters_analyzed{i,1} = data2compile{i}; %reports which data files and units were used for analysis
%     summary_clusters_analyzed{i,2} = num2cell(cluster_low);
%     summary_clusters_analyzed{i,3} = num2cell(cluster_high);
  
    for m = 1:length(cluster_low)
        [stim_cond_list, tuning_curve_fr_total, tuning_curve_fr_spon] = compute_tuning(spikes, cluster_low(m), tuning, 0 , led);%0 for do not plot
        tuning_curve_fr_evoked = tuning_curve_fr_total-tuning_curve_fr_spon; %remove 'spontanous spikes' 
        if any(tuning_curve_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
           tuning_curve_fr_evoked = tuning_curve_fr_evoked-min(tuning_curve_fr_evoked);
           fprintf('tuning_fr_evoked for jitter_low contains negative value, force minimum to be zero\n');
        end

        switch tuning
            case {'sf' , 'tf'} %normalization and calculate preffered tuning value (ie max)
                [max_fr, ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr;
                norm_tuning_curve_low = [norm_tuning_curve_low;tuning_fr_evoked_norm'];%each line is the normalized tuning for one unit
                preferred_stim_low(m,1:2)= [cluster_low(m) stim_cond_list(ind)];%last column should be full of nan (needed for ori but not sf or tf)
            case 'ori'
                %find prefered direction by converting vectors into complex numbers
                tuning_curve_fr_real = cosd(stim_cond_list)'.*tuning_curve_fr_evoked';
                tuning_curve_fr_imag = sind(stim_cond_list)'.*tuning_curve_fr_evoked';
                tuning_curve_fr_complex = complex(tuning_curve_fr_real, tuning_curve_fr_imag);
                preferred_stim_complex = sum(tuning_curve_fr_complex); %adding up all each stimulus vector (in complex number format) to find preferred direction vector
                pref_stim_angle = wrapTo2Pi(angle(preferred_stim_complex));%in rad
                
                [pref_angle_stim_cond_difference, pref_ind] = min(abs(stim_cond_list - rad2deg(pref_stim_angle)));
                pref_stim_angle = stim_cond_list(pref_ind); %deg
                max_fr = tuning_curve_fr_evoked(pref_ind);
                %[max_fr ind]= max(tuning_curve_fr_evoked);
                %pref_stim_angle = deg2rad(stim_cond_list(ind)); %%
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr; %normalizing fr for all orientations by dividing by preferred direction
                norm_tuning_curve_low = [norm_tuning_curve_low; tuning_fr_evoked_norm']; %concatylating cluter orientation fr with all other units
                
                %calculate ratio better prefered/null directions
                null_angle = wrapTo360(pref_stim_angle + 180);
                [null_angle_stim_cond_difference, null_ind] = min(abs(stim_cond_list - null_angle));
                null_fr = tuning_curve_fr_evoked(null_ind);
                ratio = (max_fr-null_fr)/(max_fr+null_fr);
                preferred_over_null_low(m,:) = [cluster_low(m) ratio];
                
                %%%remove non direction tuned cells (where max fr is super low) for quantification of preferd direction
                if ratio > 0.1
                    preferred_stim_low(m,:)= [cluster_low(m) pref_stim_angle]; %save pref angle 
                else
                    fprintf(['excluded unit: ' num2str(cluster_low(m)) '\n']);
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
                [max_fr ind]= max(tuning_curve_fr_evoked);
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr;
                norm_tuning_curve_high = [norm_tuning_curve_high;tuning_fr_evoked_norm'];%each line is the normalized tuning for one unit
                preferred_stim_high(n,1:2)= [cluster_high(n) stim_cond_list(ind)];%last column should be full of nan (needed for ori but not sf or tf)
            case 'ori'
                %find prefered direction by converting vectors into complex numbers
                tuning_curve_fr_real = cosd(stim_cond_list)'.*tuning_curve_fr_evoked';
                tuning_curve_fr_imag = sind(stim_cond_list)'.*tuning_curve_fr_evoked';
                tuning_curve_fr_complex = complex(tuning_curve_fr_real, tuning_curve_fr_imag);
                preferred_stim_complex = sum(tuning_curve_fr_complex);
                pref_stim_angle = wrapTo2Pi(angle(preferred_stim_complex));
                [pref_angle_stim_cond_difference, pref_ind] = min(abs(stim_cond_list - rad2deg(pref_stim_angle)));
                pref_stim_angle = stim_cond_list(pref_ind);
                max_fr = tuning_curve_fr_evoked(pref_ind);
                %max_fr ind]= max(tuning_curve_fr_evoked);
                %pref_stim_angle = deg2rad(stim_cond_list(ind)); %%
                tuning_fr_evoked_norm = tuning_curve_fr_evoked/max_fr; %normalizing fr for all orientations by dividing by preferred direction
                norm_tuning_curve_high = [norm_tuning_curve_high; tuning_fr_evoked_norm']; %concatylating cluter orientation fr with all other units
                
                %calculate ratio better (pref - null/pref + null) directions
                null_angle = wrapTo360(pref_stim_angle + 180);
                [null_angle_stim_cond_difference, null_ind] = min(abs(stim_cond_list - null_angle));
                null_fr = tuning_curve_fr_evoked(null_ind);
                ratio = (max_fr-null_fr)/(max_fr+null_fr);
                preferred_over_null_high(n,:) = [cluster_high(n) ratio];
                
                %%%remove non direction tuned cells (where max fr is super low) for quantification of preferd direction
                if ratio > 0.3
                    preferred_stim_high(n,:)= [cluster_high(n) pref_stim_angle]; %save pref angle 
                else
                    fprintf(['excluded unit: ' num2str(cluster_high(n)) '\n']);
                end
        end
    end
    
    %compile data from multiple recordings
    cluster_ID_low = [cluster_ID_low; cluster_low];
    cluster_ID_high = [cluster_ID_high; cluster_high];
    temp_recording_low = table(cellstr(repmat(recordinglist{i}, length(cluster_low), 1)));
    temp_recording_high = table(cellstr(repmat(recordinglist{i}, length(cluster_high), 1)));
    recording_complile_low = [recording_complile_low; temp_recording_low];
    recording_complile_high = [recording_complile_high; temp_recording_high];
    preferred_stim_compile_low=[preferred_stim_compile_low; preferred_stim_low]; 
    preferred_stim_compile_high=[preferred_stim_compile_high; preferred_stim_high];
    num_cluster_low=num_cluster_low +length(cluster_low);
    num_cluster_high=num_cluster_high +length(cluster_high);
    tuning_curve_norm_compile_low=[tuning_curve_norm_compile_low; norm_tuning_curve_low];
    tuning_curve_norm_compile_high=[tuning_curve_norm_compile_high; norm_tuning_curve_high];
    preferred_over_null_compile_low = [preferred_over_null_compile_low; preferred_over_null_low];
    preferred_over_null_compile_high=[preferred_over_null_compile_high; preferred_over_null_high];
end

%% calculate values for figures

%store values in table
%rename var to make table headers cleaner
    recording = recording_complile_low.Var1;
    cluster = cluster_ID_low;
    pref_stim = preferred_stim_compile_low(:,2);  
    tuning_curve=tuning_curve_norm_compile_low;
    DSI = preferred_over_null_compile_low(:,2);
   summary_led = table(recording, cluster, tuning_curve, pref_stim, DSI); %create table
    recording = recording_complile_high.Var1;
    cluster = cluster_ID_high; %first column is cluster ID
    pref_stim= preferred_stim_compile_high(:,2);
    tuning_curve= tuning_curve_norm_compile_high;
    DSI = preferred_over_null_compile_high(:,2);
    summary_non_led = table(recording, cluster, tuning_curve, pref_stim, DSI);

    
%bar graph of distribution
    
% calculating tuning distribution 
percent_pref_stim_low=[stim_cond_list zeros(length(stim_cond_list),1)];
percent_pref_stim_high=[stim_cond_list zeros(length(stim_cond_list),1)];

%remove nan %%%NEED TO REMOVE INF
preferred_stim_compile_low= [preferred_stim_compile_low(~isnan(preferred_stim_compile_low(:,1)), 1) preferred_stim_compile_low(~isnan(preferred_stim_compile_low(:,2)),2)];
preferred_stim_compile_high= [preferred_stim_compile_high(~isnan(preferred_stim_compile_high(:,1)), 1) preferred_stim_compile_high(~isnan(preferred_stim_compile_high(:,2)),2)];

%number of units with preferred tuning and remove nan that remain from initialization since only units with preferred direction are included here 
num_pref_cluster_low = length(preferred_stim_compile_low);
num_pref_cluster_high = length(preferred_stim_compile_high);

% calculate percentage of prefered sf or ori
for n=1:length(stim_cond_list)
    percent_pref_stim_low(n,2)=sum(preferred_stim_compile_low(:,2)== stim_cond_list(n))/num_pref_cluster_low*100;%the number of clusters that preferred each orientation
    percent_pref_stim_high(n,2)=sum(preferred_stim_compile_high(:,2)== stim_cond_list(n))/num_pref_cluster_high*100;
end

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
        
        theta = stim_cond_list'/180*pi;
        theta = [theta theta(1)];
        
%         subplot(1,2,1)%plot distribution
%         %by default clusters with jitter<2ms is blue, the clusters with jitter>2ms is red
%         tuning_distribution_low=[tuning_distribution_low; tuning_distribution_low(1,:)];
%         tuning_distribution_high=[tuning_distribution_high;  tuning_distribution_high(1,:)];
        
%         %the axes of the polar plot is set by the first plotted curve
%         if max(tuning_distribution_high(:,2))>max(tuning_distribution_low(:,2))
%             polar(theta, tuning_distribution_high(:,2)','r')
%             hold on
%             polar(theta, tuning_distribution_low(:,2)','b')
%         else
%             polar(theta, tuning_distribution_low(:,2)','b')
%             hold on
%             polar(theta, tuning_distribution_high(:,2)','r')
%         end
        
%         subplot(1,2,1)%plot distribution
%         polarhistogram((tuning_distribution_low(:,2)), 24, 'FaceColor', 'r' )
%         hold on
%         polarhistogram(deg2rad(max_tuning_compile_low(:,2)), 24, 'FaceColor', 'b') 
%         title('preferred tuning distribution (% units)')
        C = categorical({'0', '30', '60', '90',  '120', '150', '180', '210', '240', '270', '300', '360'}); %add zero so that closes cirle
        subplot(1,3,1)
%         histogram(preferred_stim_compile_high(:,2),'NumBins', 12, 'BinLimits', [0 360], 'FaceColor', 'r');
%         hold on
%         histogram(preferred_stim_compile_low(:,2),'NumBins', 12, 'BinLimits', [0 360], 'FaceColor', 'b');
%          xlabel('Degrees', 'FontSize', 25)
%         ylabel('# units', 'FontSize', 25)
%         title('preferred orientation')
        hold on
        bar(stim_cond_list',[percent_pref_stim_low(:,2) percent_pref_stim_high(:,2)])
        xlabel('Preferred direction (degree)', 'FontSize', 15)
        ylabel('% units', 'FontSize', 20)
        title('% of perferred direction', 'FontSize', 15)
        
        %for grant:
%         figure;
%         bar([percent_pref_stim_low(:,2))

        subplot(1,3,2);%plot average tuning
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
        
        %(pref - null/ pref + null) ratio; 
        preferred_over_null_compile_high_mean=nanmean(preferred_over_null_compile_high(:,2));
        preferred_over_null_compile_high_std=nanstd(preferred_over_null_compile_high(:,2));
        preferred_over_null_compile_high_stem=preferred_over_null_compile_high_std/sqrt(length(preferred_over_null_compile_high(:,2)));
        
        preferred_over_null_compile_low_mean=nanmean(preferred_over_null_compile_low(:,2));
        preferred_over_null_compile_low_std=nanstd(preferred_over_null_compile_low(:,2));
        preferred_over_null_compile_low_stem=preferred_over_null_compile_low_std/sqrt(length(preferred_over_null_compile_low(:,2)));

        [H, P] =ttest2(preferred_over_null_compile_high(:,2), preferred_over_null_compile_low(:,2));
       
        %plot direction selectivity ratio 
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
            text(-0.4, max(preferred_over_null_compile_high_mean + preferred_over_null_compile_high_stem*2), '----------- * -----------')
        else
            text(-0.4, max(preferred_over_null_compile_high_mean + preferred_over_null_compile_high_stem*2), '----------- n.s.-----------')
        end

       save(['E:\extarcellular\Andreanne\data_analysis\photo-tag\population tuning\summary_led_' tuning datestr(now, 'dd-mmm-yyyy') user '.mat'], 'summary_led', 'summary_non_led', ...
            'percent_pref_stim_low', 'percent_pref_stim_high', 'stim_cond_list', 'ave_tuning_norm_compile_low', 'ave_tuning_norm_compile_high', 'jitter_cutoff');
        
end
end

