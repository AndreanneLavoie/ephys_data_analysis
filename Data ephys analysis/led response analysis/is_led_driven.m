function spikes = is_led_driven(spikes, alpha_val)
%calculate wether units show any led driven activity using paired t-test
%for all trials (sum up all led pulses within trial)

%INPUT
%{1} spikes construct
%{2} alpha (t-test cutoff criteria) - i.e. 0.05 - 0.001

%OUTPUT
%{1} updated spikes contruct with new field : spikes.is_led_driven; 
%with 1st column = good units; 2nd = 0 or 1 logic value for is_led_driven

%assign = 81;
bin = 0.001; %sec 
bin_limit = [0.002 0.012]; %sec
led_window_delta = bin_limit(2) - bin_limit(1);
plot_flag = 0;
led_baseline_delta = 0.3; %sec
num_trials = length(spikes.vs_params(:, 16));
num_led_pulses = spikes.vs_params(1, 16);

% 1 - calculate var for average baseline firing rate right before led stim for each trial
led_delay = spikes.vs_params(1,11)/1000; %sec
baseline_start = led_delay - led_baseline_delta; %sec %baseline is 500ms before led on (give time for network to calm down form vis stim)
baseline_range = [baseline_start led_delay];

% 2 - identify good units
spikes.is_led_driven = spikes.labels(spikes.labels(:,2) == 2);
num_units = length(spikes.is_led_driven);

% 3 - each units, a)calculate baseline fr for each trial
for i = 1:num_units 
    filtered_spikes = filtspikes(spikes, 0, 'assigns', spikes.is_led_driven(i, 1));
    
    % compute baseline fr for each trial (third var output)
    [baseline_fr_avg, baseline_fr_sem, baseline_fr] = computeFR(filtered_spikes, baseline_range); %spikes/sec; third var is baseline fr for each trial
    
    % compute led fr for each trial
    led_fr = nan(size(baseline_fr));
    for j = 1:num_trials
        temp_spikes =  filtspikes(filtered_spikes, 0, 'trials', j); %this spikes contruct is technically filt by assign and trial 
        
        % compute histogram bin values for each trial (i.e. make hist bins of 
        % all spiketimes for all 10 led pulses) for eachassign  
        [hist_fr_bins, edges_temp, bin_temp]=histcounts(temp_spikes.led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);
        led_fr(j) = sum(hist_fr_bins)/(led_window_delta*num_led_pulses);
    end
    
    % run paired t-test between baseline and led fr (for each trial)
    [h,p] = ttest(baseline_fr,led_fr, 'Alpha', alpha_val);
    spikes.is_led_driven(i, 2) = h;
    %change nan (no spikes during led pulse window == non led-driven unit
    spikes.is_led_driven(isnan(spikes.is_led_driven(:,2)),2) = 0;
    

end


end


