function hist_width = calc_jitters(spikes, assign)
%calculate the jitter using histogram
%assign = 81;
bin = 1; %ms 
bin_limit = [0 25]; %ms
plot_flag = 0;
led_baseline_delta = 0.3; %sec
num_trials = length(spikes.vs_params(:, 16));
num_led_pulses = spikes.vs_params(1, 16);

% 1 - calculate average baseline firing rate right before led stim for each trial
led_delay = spikes.vs_params(1,11)/1000; %sec
baseline_start = led_delay - led_baseline_delta; %sec %baseline is 500ms before led on (give time for network to calm down form vis stim)
baseline_range = [baseline_start led_delay];

filtered_spikes = filtspikes(spikes, 0, 'assigns', assign);
led_baseline_fr = computeFR(filtered_spikes, baseline_range); %spikes/sec

% 2 - convert number to bin size (ie average number of spontaneous spikes / bin)
ave_num_spont_spikes_per_bin = led_baseline_fr*(bin/1000); %convert to spikes/bin (sec cancels out)
%multiply the ave number of spikes per bin by the total number of trials*numberof led pulses
sum_num_spont_spikes_per_bin = ave_num_spont_spikes_per_bin*num_trials*num_led_pulses;

% compute histogram bin values for led -> USE TOTAL NUMBER OF SPIKES 
[length_led, led_spiketimes] = make_first_led_spiketimes(spikes, assign, bin, bin_limit, 1, plot_flag); 
led_spiketimes = led_spiketimes*1000; %convert to ms
[hist_fr_bins, edges_temp, bin_temp]=histcounts(led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);

% 3 - substract average spontaneous spikes number per each bin of led histogram
hist_fr_bins_sub = hist_fr_bins - sum_num_spont_spikes_per_bin;



%check
%figure;
%cen = edges_temp(1:end-1) + mean(diff(edges_temp))/2;
%bar(cen, hist_fr_bins_sub, 1)

%remove neg values from hist_fr_bins_sub
hist_fr_bins_sub = hist_fr_bins_sub(hist_fr_bins_sub > 0);

if isempty(hist_fr_bins_sub)
    hist_width = nan;
else
    % 4 - calculate number of spikes * binwidth for each bin and sum them all up
    % 5 - divide by the peak (max number of spikes in any bin) 
    peak_bin = max(hist_fr_bins_sub);
    hist_width = sum(hist_fr_bins_sub*bin)/peak_bin;
end

end


