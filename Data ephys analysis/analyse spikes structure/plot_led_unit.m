function photo_tagging_param = plot_led_unit(spikes, assign, tuning, flag_plot)
%% This function will plot the psth, led histograms and tuning sf psth for the given unit
%INPUT
%spikes structure 
%assign -> integer (can be multiple values e.g. [12 14 46])
%tuning = 'sf'; or 'none'

%OUTPUT
%photo_tagging_param stucture containing:
    %spiketimes window near LED (default 15ms)
    %extract data from psth
    %mean latency 
    %jitter (SD)
    %success rate / reliability
    %First vs all spiketimes ratio
    %LED frequency vs background frquency (ratio)

%% parameters
bin_LED = 1; %0.5 %ms
bin_psth = 20;%20; %assumes 20 ms when calculating led-over-background firing rate ratio
tagging_analysis_window = 20; %15 ms %To calculate jitter and mean
bin_limit = [0 200]; %in ms; for led spiketimes histogram

%led over background timing parameters
num_trials = length(spikes.vs_params(:,1)); %number of visual stim trials
num_LED_pulses = spikes.vs_params(1,16); %number of led pulses per trial
led_period = spikes.vs_params(1,13);
led_delay = spikes.vs_params(1,11); %start of led train after each trigger for vis stim
led_end = led_delay + led_period*num_LED_pulses;
background_fr_range_delta = 300; %ms
mean_background_fr_range = [(led_delay - background_fr_range_delta) led_delay] ; %ms

%visual driven activity timing parameters
vis_stim_start = spikes.vs_params(1,9); %ms; trig time
vis_stim_dur = spikes.vs_params(1,7); %ms; vis stim duration
vis_stim_end = vis_stim_start + vis_stim_dur;
trial_dur = vis_stim_start + vis_stim_dur + spikes.vs_params(1,6) + spikes.vs_params(1,8); %add tbef and taft

%filter by unit
filt_spikes_unit = filtspikes(spikes, 0, 'assigns', assign);

%determine subplot parameters (n, m and p) ie; n x m matrix and p is the
%coordinate of the figure; this depends if tuning is added 
if tuning == 'sf'   
    subplot_param_n_m = [3 5]; % 3 x 5 grid
    subplot_p_psth = [1 2 3 4 5];
    subplot_p_led_spiketimes_1 = 6;
    subplot_p_led_spiketimes_all = 7;
    
elseif tuning == 'ori'
    
    pass
   
else
    subplot_param_n_m = [2,2];
    subplot_p_psth = [1:subplot_param_n_m(2)];
    subplot_p_led_spiketimes_1 = 3;
    subplot_p_led_spiketimes_all = 4;
end



%% PLOT
figure;

%plot psth (from trigger to trigger - includes vis stim and LED)
subplot(subplot_param_n_m(1), subplot_param_n_m(2), subplot_p_psth);
title('Average PSTH');
filt_unit_psth = psth_bl(filt_spikes_unit, bin_psth);


%plot histogram of first led spiketimes
subplot(subplot_param_n_m(1), subplot_param_n_m(2), subplot_p_led_spiketimes_1);
[length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, assign, bin_LED, bin_limit, 1, flag_plot);
title('1st LED st')

%plot histogram of all led spikes
subplot(subplot_param_n_m(1), subplot_param_n_m(2), subplot_p_led_spiketimes_all);
[length_led_all, led_spiketimes_all, led_hist_all] = make_first_led_spiketimes(spikes, assign, bin_LED, bin_limit, 0, flag_plot);
title('all LED st')

%plot tuning 
if tuning == 'sf'
    
    %plot tuning curve
    subplot(subplot_param_n_m(1), subplot_param_n_m(2), subplot_p_led_spiketimes_all + 2);
    [stim_cond_list, tuning_fr] = compute_tuning(spikes, assign, 'sf', 1);
    title('SF tuning curve');
    
    axis_max = stim_cond_list; %initalize param
    
    for i=(1 : length(stim_cond_list))
        
        filt_spikes_sf = filtspikes(spikes, 0, 'assigns', assign, 'stimcond', single(stim_cond_list(i)));
        subplot(subplot_param_n_m(1), subplot_param_n_m(2),(10 + i));
        h = psth_bl(filt_spikes_sf);
        title(['SF :' num2str(stim_cond_list(i))]);
        
        %calculate max value for each sf
        axis_max(i) = max(h.YData);
    end
    
    %replot everything with new axis limit
    y_min = 0;
    y_max = ceil(max(axis_max));
    
    for i=(1 : length(stim_cond_list))
        
        filt_spikes_sf = filtspikes(spikes, 0, 'assigns', assign, 'stimcond', single(stim_cond_list(i)));
        subplot(subplot_param_n_m(1), subplot_param_n_m(2),(10 + i));
        h = psth_bl(filt_spikes_sf);
        title(['SF :' num2str(stim_cond_list(i))]);
        
        %set a constant y-axis limit 
        ylim([y_min y_max]);
 
    end

end

%general title
suptitle(['Unit ' num2str(assign)])

%% CALUCULATE 'OPTO-TAG' PARAMETERS
photo_tagging_param.unit = assign;
photo_tagging_param.window = tagging_analysis_window; %sec

%remove first led spiketimes over analysis window (typically 15 ms) 
photo_tagging_param.led_spiketimes = led_spiketimes_1(led_spiketimes_1 < tagging_analysis_window/1000); 

%mean latency %%%SHOULD BE PEAK? WHERE BIGGEST Y VALUE ; ONLY WORKS IF
%SMOOTH PSTH CURVE
photo_tagging_param.mean_latency = mean(photo_tagging_param.led_spiketimes);

%jitter (SD)
photo_tagging_param.jitter = std(photo_tagging_param.led_spiketimes);

%% success rate / reliability %%%IS THIS CORRECT???
    %over whole recording period
photo_tagging_param.reliability = length(photo_tagging_param.led_spiketimes)/(num_trials*num_LED_pulses); 
    %during range of unit activity -> for unstable recording where units
    %THAT do not last whole recording periods 
    num_trials_with_unit_activity = length(unique(filt_spikes_unit.trials));  
    photo_tagging_param.percentage_active_trials = num_trials_with_unit_activity / length(unique(spikes.trials))*100; %trials with spikes over total number of trials

%% First vs all led_spiketimes ratio
number_led_spiketimes_1 = length(photo_tagging_param.led_spiketimes);
number_led_spiketimes_all = length(led_spiketimes_all(led_spiketimes_all < (tagging_analysis_window/1000)));
photo_tagging_param.first_over_all = number_led_spiketimes_1/number_led_spiketimes_all;

%% LED frequency vs background frquency (ratio)
trial_firing_rates = filt_unit_psth.YData; %20 ms (or value of bin_psth) binned firing rate for each trial 0 = trig time (includes vis stimand led events)

%background fr (range: between vis stim and first led pulse in a trial)
%%%%NOT USING ANYMORE
mean_background_fr_index = ((mean_background_fr_range(1)/bin_psth) + 1):((mean_background_fr_range(2)/bin_psth));
photo_tagging_param.background_fr = trial_firing_rates(mean_background_fr_index);
photo_tagging_param.mean_background_fr = mean(photo_tagging_param.background_fr);

%led fr 
led_firing_rates_index = ((led_delay/bin_psth) + 1): (led_period/bin_psth):((led_end/bin_psth)); %first 20ms firing rate after led pulse
photo_tagging_param.led_pulse_fr = trial_firing_rates(led_firing_rates_index)';
photo_tagging_param.mean_led_fr = mean(photo_tagging_param.led_pulse_fr); %average led firing rate (20ms bin mean)

%background led fr -> take the average fr of three 20ms bins before led
background_fr_led_pulse= zeros(length(led_firing_rates_index), 3);
background_fr_led_pulse(:,1) =  trial_firing_rates(led_firing_rates_index - 1); 
background_fr_led_pulse(:,2) =  trial_firing_rates(led_firing_rates_index - 2);
background_fr_led_pulse(:,3) =  trial_firing_rates(led_firing_rates_index - 3);
photo_tagging_param.background_fr_led_pulse = mean(background_fr_led_pulse, 2); %average of 3*20ms firing rate before led pulse

photo_tagging_param.background_fr_led_pulse_mean = mean(photo_tagging_param.background_fr_led_pulse);

%paired t-test between fr before led and after led
% [photo_tagging_param.reject_null, photo_tagging_param.p] = ttest(photo_tagging_param.led_pulse_fr, photo_tagging_param.background_fr_led_pulse);

%ratio
photo_tagging_param.led_over_background = mean(photo_tagging_param.led_pulse_fr) / mean(photo_tagging_param.background_fr_led_pulse);
%photo_tagging_param.led_over_background = photo_tagging_param.mean_led_fr / photo_tagging_param.mean_background_fr;

%% NEW LED vs BACKGROUND TTEST

photo_tagging_param.background_fr_trial = zeros(num_trials, 1);
photo_tagging_param.led_2sec_fr_trial = zeros(num_trials, 1);

for t=1:num_trials
    filtspikes_unit_trial = filtspikes(spikes, 0, 'assigns', assign, 'trials', t);
    
    photo_tagging_param.background_fr_trial(t) = computeFR(filtspikes_unit_trial, mean_background_fr_range/1000);
    photo_tagging_param.led_2sec_fr_trial(t) = computeFR(filtspikes_unit_trial, [led_delay led_end]/1000);
end

[photo_tagging_param.led_driven_ttest_h, photo_tagging_param.led_driven_ttest_p] = ttest(photo_tagging_param.led_2sec_fr_trial, photo_tagging_param.background_fr_trial);

%% visual driven activity
vis_stim_fr_index = ((vis_stim_start/bin_psth) + 1): ((vis_stim_end/bin_psth)); %first 20ms firing rate after led pulse
photo_tagging_param.mean_vis_stim_fr = mean(trial_firing_rates(vis_stim_fr_index)); %average led firing rate (20ms bin mean)
    %ratio
photo_tagging_param.visual_driven_activity = photo_tagging_param.mean_vis_stim_fr / photo_tagging_param.mean_background_fr;